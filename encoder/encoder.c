/*****************************************************************************
 * encoder.c: top-level encoder functions
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at licensing@x264.com.
 *****************************************************************************/

#include "common/common.h"

#include "set.h"
#include "analyse.h"
#include "ratecontrol.h"
#include "macroblock.h"
#include "me.h"

#if HAVE_VISUALIZE
#include "common/visualize.h"
#endif

//#define DEBUG_MB_TYPE

#define bs_write_ue bs_write_ue_big

static int x264_encoder_frame_end( x264_t *h,x264_nal_t **pp_nal, int *pi_nal,
                                   x264_picture_t *pic_out );

/****************************************************************************
 *
 ******************************* x264 libs **********************************
 *
 ****************************************************************************/
static double x264_psnr( double sqe, double size )
{
    double mse = sqe / (PIXEL_MAX*PIXEL_MAX * size);
    if( mse <= 0.0000000001 ) /* Max 100dB */
        return 100;

    return -10.0 * log10( mse );
}

static double x264_ssim( double ssim )
{
    double inv_ssim = 1 - ssim;
    if( inv_ssim <= 0.0000000001 ) /* Max 100dB */
        return 100;

    return -10.0 * log10( inv_ssim );
}

static void x264_frame_dump( x264_t *h )
{
    FILE *f = x264_fopen( h->param.psz_dump_yuv, "r+b" );
    if( !f )
        return;

    /* Write the frame in display order */
    int frame_size = FRAME_SIZE( h->param.i_height * h->param.i_width * sizeof(pixel) );
    fseek( f, (uint64_t)h->fdec->i_frame * frame_size, SEEK_SET );
    for( int p = 0; p < (1); p++ )
        for( int y = 0; y < h->param.i_height; y++ )
            fwrite( &h->fdec->plane[p][y*h->fdec->i_stride[p]], sizeof(pixel), h->param.i_width, f );
    {
        int cw = h->param.i_width>>1;
        int ch = h->param.i_height>>CHROMA_V_SHIFT;
        pixel *planeu = x264_malloc( (cw*ch*2+32)*sizeof(pixel) );
        pixel *planev = planeu + cw*ch + 16;
        h->mc.plane_copy_deinterleave( planeu, cw, planev, cw, h->fdec->plane[1], h->fdec->i_stride[1], cw, ch );
        fwrite( planeu, 1, cw*ch*sizeof(pixel), f );
        fwrite( planev, 1, cw*ch*sizeof(pixel), f );
        x264_free( planeu );
    }
    fclose( f );
}

/* Fill "default" values */
static void x264_slice_header_init( x264_t *h, x264_slice_header_t *sh,
                                    x264_sps_t *sps, x264_pps_t *pps,
                                    int i_idr_pic_id, int i_frame, int i_qp )
{
    x264_param_t *param = &h->param;

    /* First we fill all fields */
    sh->sps = sps;
    sh->pps = pps;

    sh->i_first_mb  = 0;
    sh->i_last_mb   = h->mb.i_mb_count - 1;
    sh->i_pps_id    = pps->i_id;

    sh->i_frame_num = i_frame;

    sh->b_mbaff = 0;
    sh->b_field_pic = 0;    /* no field support for now */
    sh->b_bottom_field = 0; /* not yet used */

    sh->i_idr_pic_id = i_idr_pic_id;

    /* poc stuff, fixed later */
    sh->i_poc = 0;
    sh->i_delta_poc_bottom = 0;
    sh->i_delta_poc[0] = 0;
    sh->i_delta_poc[1] = 0;

    sh->i_redundant_pic_cnt = 0;

    sh->b_num_ref_idx_override = 0;
    sh->i_num_ref_idx_l0_active = 1; //use p ref frame
    sh->i_num_ref_idx_l1_active = 0; //do not use b ref frame

    sh->b_ref_pic_list_reordering[0] = h->b_ref_reorder;
    sh->b_ref_pic_list_reordering[1] = 0;

    /* If the ref list isn't in the default order, construct reordering header */
    {
        if( sh->b_ref_pic_list_reordering[0] )
        {
            int pred_frame_num = i_frame;
            for( int i = 0; i < h->i_ref; i++ )
            {
                int diff = h->fref[i]->i_frame_num - pred_frame_num;
                sh->ref_pic_list_order[0][i].idc = ( diff > 0 );
                sh->ref_pic_list_order[0][i].arg = (abs(diff) - 1) & ((1 << sps->i_log2_max_frame_num) - 1);
                pred_frame_num = h->fref[i]->i_frame_num;
            }
        }
    }

    sh->i_qp = SPEC_QP(i_qp);
    sh->i_qp_delta = sh->i_qp - pps->i_pic_init_qp;
    sh->b_sp_for_swidth = 0;
    sh->i_qs_delta = 0;

    int deblock_thresh = i_qp + 2 * X264_MIN(param->i_deblocking_filter_alphac0, param->i_deblocking_filter_beta);
    /* If effective qp <= 15, deblocking would have no effect anyway */
    if( param->b_deblocking_filter && (h->mb.b_variable_qp || 15 < deblock_thresh ) )
        sh->i_disable_deblocking_filter_idc = 0;
    else
        sh->i_disable_deblocking_filter_idc = 1;
    sh->i_alpha_c0_offset = param->i_deblocking_filter_alphac0 << 1;
    sh->i_beta_offset = param->i_deblocking_filter_beta << 1;
}

static void x264_slice_header_write( bs_t *s, x264_slice_header_t *sh, int i_nal_ref_idc )
{
    bs_write_ue( s, sh->i_first_mb );

    bs_write_ue( s, sh->i_type + 5 );   /* same type things */
    bs_write_ue( s, sh->i_pps_id );
    bs_write( s, sh->sps->i_log2_max_frame_num, sh->i_frame_num & ((1<<sh->sps->i_log2_max_frame_num)-1) );

    if( sh->i_idr_pic_id >= 0 ) /* NAL IDR */
        bs_write_ue( s, sh->i_idr_pic_id );

    if( sh->sps->i_poc_type == 0 )
    {
        bs_write( s, sh->sps->i_log2_max_poc_lsb, sh->i_poc & ((1<<sh->sps->i_log2_max_poc_lsb)-1) );
        if( sh->pps->b_pic_order && !sh->b_field_pic )
            bs_write_se( s, sh->i_delta_poc_bottom );
    }

    if( sh->pps->b_redundant_pic_cnt )
        bs_write_ue( s, sh->i_redundant_pic_cnt );

    if( sh->i_type == SLICE_TYPE_P )
    {
        bs_write1( s, sh->b_num_ref_idx_override );
        if( sh->b_num_ref_idx_override )
        {
            bs_write_ue( s, sh->i_num_ref_idx_l0_active - 1 );
        }
    }

    /* ref pic list reordering */
    if( sh->i_type != SLICE_TYPE_I )
    {
        bs_write1( s, sh->b_ref_pic_list_reordering[0] );
        if( sh->b_ref_pic_list_reordering[0] )
        {
            for( int i = 0; i < sh->i_num_ref_idx_l0_active; i++ )
            {
                bs_write_ue( s, sh->ref_pic_list_order[0][i].idc );
                bs_write_ue( s, sh->ref_pic_list_order[0][i].arg );
            }
            bs_write_ue( s, 3 );
        }
    }

    if( i_nal_ref_idc != 0 )
    {
        if( sh->i_idr_pic_id >= 0 )
        {
            bs_write1( s, 0 );  /* no output of prior pics flag */
            bs_write1( s, 0 );  /* long term reference flag */
        }
        else
        {
            bs_write1( s, sh->i_mmco_command_count > 0 ); /* adaptive_ref_pic_marking_mode_flag */
            if( sh->i_mmco_command_count > 0 )
            {
                for( int i = 0; i < sh->i_mmco_command_count; i++ )
                {
                    bs_write_ue( s, 1 ); /* mark short term ref as unused */
                    bs_write_ue( s, sh->mmco[i].i_difference_of_pic_nums - 1 );
                }
                bs_write_ue( s, 0 ); /* end command list */
            }
        }
    }

    bs_write_se( s, sh->i_qp_delta );      /* slice qp delta */

    if( sh->pps->b_deblocking_filter_control )
    {
        bs_write_ue( s, sh->i_disable_deblocking_filter_idc );
        if( sh->i_disable_deblocking_filter_idc != 1 )
        {
            bs_write_se( s, sh->i_alpha_c0_offset >> 1 );
            bs_write_se( s, sh->i_beta_offset >> 1 );
        }
    }
}

/* If we are within a reasonable distance of the end of the memory allocated for the bitstream, */
/* reallocate, adding an arbitrary amount of space. */
static int x264_bitstream_check_buffer_internal( x264_t *h, int size, int i_nal )
{
    if( (h->out.bs.p_end - h->out.bs.p < size) )
    {
        int buf_size = h->out.i_bitstream + size;
        uint8_t *buf = x264_malloc( buf_size );
        if( !buf )
            return -1;
        int aligned_size = h->out.i_bitstream & ~15;
        h->mc.memcpy_aligned( buf, h->out.p_bitstream, aligned_size );
        memcpy( buf + aligned_size, h->out.p_bitstream + aligned_size, h->out.i_bitstream - aligned_size );

        intptr_t delta = buf - h->out.p_bitstream;

        h->out.bs.p_start += delta;
        h->out.bs.p += delta;
        h->out.bs.p_end = buf + buf_size;

        for( int i = 0; i <= i_nal; i++ )
            h->out.nal[i].p_payload += delta;

        x264_free( h->out.p_bitstream );
        h->out.p_bitstream = buf;
        h->out.i_bitstream = buf_size;
    }
    return 0;
}

static int x264_bitstream_check_buffer( x264_t *h )
{
    int max_row_size = (2500) * h->mb.i_mb_width;
    return x264_bitstream_check_buffer_internal( h, max_row_size, h->out.i_nal );
}

static int x264_bitstream_check_buffer_filler( x264_t *h, int filler )
{
    filler += 32; // add padding for safety
    return x264_bitstream_check_buffer_internal( h, filler, -1 );
}

/****************************************************************************
 *
 ****************************************************************************
 ****************************** External API*********************************
 ****************************************************************************
 *
 ****************************************************************************/
//b_open - when open is 1, when reconfig is 0
static int x264_validate_parameters( x264_t *h, int b_open )
{
    if( h->param.i_width <= 0 || h->param.i_height <= 0 )
    {
        x264_log( h, X264_LOG_ERROR, "invalid width x height (%dx%d)\n",
                  h->param.i_width, h->param.i_height );
        return -1;
    }

    int i_csp = h->param.i_csp & X264_CSP_MASK;

    if( i_csp <= X264_CSP_NONE || i_csp >= X264_CSP_MAX )
    {
        x264_log( h, X264_LOG_ERROR, "invalid CSP (only I420/YV12/NV12 supported)\n" );
        return -1;
    }

    if( i_csp < X264_CSP_NV12 && h->param.i_width % 2 )
    {
        x264_log( h, X264_LOG_ERROR, "width not divisible by 2 (%dx%d)\n",
                  h->param.i_width, h->param.i_height );
        return -1;
    }

    if( i_csp < X264_CSP_NV12 && h->param.i_height % 2 )
    {
        x264_log( h, X264_LOG_ERROR, "height not divisible by 2 (%dx%d)\n",
                  h->param.i_width, h->param.i_height );
        return -1;
    }

    if( (h->param.crop_rect.i_left + h->param.crop_rect.i_right ) >= h->param.i_width ||
        (h->param.crop_rect.i_top  + h->param.crop_rect.i_bottom) >= h->param.i_height )
    {
        x264_log( h, X264_LOG_ERROR, "invalid crop-rect %u,%u,%u,%u\n", h->param.crop_rect.i_left,
                  h->param.crop_rect.i_top, h->param.crop_rect.i_right,  h->param.crop_rect.i_bottom );
        return -1;
    }

    h->param.i_keyint_max = x264_clip3( h->param.i_keyint_max, 1, X264_KEYINT_MAX_INFINITE );
    if( h->param.i_keyint_max == 1 )
    {
        h->param.b_intra_refresh = 0;
        h->param.i_frame_reference = 1;
        h->param.i_dpb_size = 1;
    }

    /* Detect default ffmpeg settings and terminate with an error. */
    if( b_open )
    {
        int score = 0;
        score += h->param.analyse.i_me_range == 0;
        score += h->param.rc.i_qp_step == 3;
        score += h->param.i_keyint_max == 12;
        score += h->param.rc.i_qp_min == 2;
        score += h->param.rc.i_qp_max == 31;
        score += h->param.rc.f_qcompress == 0.5;
        score += fabs(h->param.rc.f_ip_factor - 1.25) < 0.01;
        score += fabs(h->param.rc.f_pb_factor - 1.25) < 0.01;
        score += h->param.analyse.inter == 0 && h->param.analyse.i_subpel_refine == 8;
        if( score >= 5 )
        {
            x264_log( h, X264_LOG_ERROR, "broken ffmpeg default settings detected\n" );
            x264_log( h, X264_LOG_ERROR, "use an encoding preset (e.g. -vpre medium)\n" );
            x264_log( h, X264_LOG_ERROR, "preset usage: -vpre <speed> -vpre <profile>\n" );
            x264_log( h, X264_LOG_ERROR, "speed presets are listed in x264 --help\n" );
            x264_log( h, X264_LOG_ERROR, "profile is optional; x264 defaults to high\n" );
            return -1;
        }
    }

    if( h->param.rc.i_rc_method < 0 || h->param.rc.i_rc_method > 2 )
    {
        x264_log( h, X264_LOG_ERROR, "no ratecontrol method specified\n" );
        return -1;
    }

    h->param.rc.f_rf_constant = x264_clip3f( h->param.rc.f_rf_constant, 0, 51 );
    h->param.rc.f_rf_constant_max = x264_clip3f( h->param.rc.f_rf_constant_max, 0, 51 );
    h->param.rc.i_qp_constant = x264_clip3( h->param.rc.i_qp_constant, 0, QP_MAX );
    h->param.analyse.i_subpel_refine = x264_clip3( h->param.analyse.i_subpel_refine, 0, 11 );
    h->param.rc.f_ip_factor = X264_MAX( h->param.rc.f_ip_factor, 0.01f );
    h->param.rc.f_pb_factor = X264_MAX( h->param.rc.f_pb_factor, 0.01f );
    if( h->param.rc.i_rc_method == X264_RC_CRF )
    {
        h->param.rc.i_qp_constant = h->param.rc.f_rf_constant ;
        h->param.rc.i_bitrate = 0;
    }
    if( b_open && (h->param.rc.i_rc_method == X264_RC_CQP || h->param.rc.i_rc_method == X264_RC_CRF)
        && h->param.rc.i_qp_constant == 0 )
    {
        h->mb.b_lossless = 1;
        h->param.rc.i_rc_method = X264_RC_CQP;
        h->param.rc.f_ip_factor = 1;
        h->param.rc.f_pb_factor = 1;
        h->param.analyse.b_psnr = 0;
        h->param.analyse.b_ssim = 0;
        h->param.analyse.i_chroma_qp_offset = 0;
        h->param.analyse.b_fast_pskip = 0;
        h->param.analyse.i_noise_reduction = 0;
        h->param.analyse.b_psy = 0;
    }
    if( h->param.rc.i_rc_method == X264_RC_CQP )
    {
        float qp_p = h->param.rc.i_qp_constant;
        float qp_i = qp_p - 6*log2f( h->param.rc.f_ip_factor );
        float qp_b = qp_p + 6*log2f( h->param.rc.f_pb_factor );
        h->param.rc.i_qp_min = x264_clip3( (int)(X264_MIN3( qp_p, qp_i, qp_b )), 0, QP_MAX );
        h->param.rc.i_qp_max = x264_clip3( (int)(X264_MAX3( qp_p, qp_i, qp_b ) + .999), 0, QP_MAX );
        h->param.rc.i_aq_mode = 0;
        h->param.rc.i_bitrate = 0;
    }
    h->param.rc.i_qp_max = x264_clip3( h->param.rc.i_qp_max, 0, QP_MAX );
    h->param.rc.i_qp_min = x264_clip3( h->param.rc.i_qp_min, 0, h->param.rc.i_qp_max );
    h->param.rc.i_qp_step = x264_clip3( h->param.rc.i_qp_step, 2, QP_MAX );
    h->param.rc.i_bitrate = x264_clip3( h->param.rc.i_bitrate, 0, 2000000 );
    if( h->param.rc.i_rc_method == X264_RC_ABR && !h->param.rc.i_bitrate )
    {
        x264_log( h, X264_LOG_ERROR, "bitrate not specified\n" );
        return -1;
    }
    h->param.rc.i_vbv_buffer_size = x264_clip3( h->param.rc.i_vbv_buffer_size, 0, 2000000 );
    h->param.rc.i_vbv_max_bitrate = x264_clip3( h->param.rc.i_vbv_max_bitrate, 0, 2000000 );
    h->param.rc.f_vbv_buffer_init = x264_clip3f( h->param.rc.f_vbv_buffer_init, 0, 2000000 );
    if( h->param.rc.i_vbv_buffer_size )
    {
        if( h->param.rc.i_rc_method == X264_RC_CQP )
        {
            x264_log( h, X264_LOG_WARNING, "VBV is incompatible with constant QP, ignored.\n" );
            h->param.rc.i_vbv_max_bitrate = 0;
            h->param.rc.i_vbv_buffer_size = 0;
        }
        else if( h->param.rc.i_vbv_max_bitrate == 0 )
        {
            if( h->param.rc.i_rc_method == X264_RC_ABR )
            {
                x264_log( h, X264_LOG_WARNING, "VBV maxrate unspecified, assuming CBR\n" );
                h->param.rc.i_vbv_max_bitrate = h->param.rc.i_bitrate;
            }
            else
            {
                x264_log( h, X264_LOG_WARNING, "VBV bufsize set but maxrate unspecified, ignored\n" );
                h->param.rc.i_vbv_buffer_size = 0;
            }
        }
        else if( h->param.rc.i_vbv_max_bitrate < h->param.rc.i_bitrate &&
                 h->param.rc.i_rc_method == X264_RC_ABR )
        {
            x264_log( h, X264_LOG_WARNING, "max bitrate less than average bitrate, assuming CBR\n" );
            h->param.rc.i_bitrate = h->param.rc.i_vbv_max_bitrate;
        }
    }
    else if( h->param.rc.i_vbv_max_bitrate )
    {
        x264_log( h, X264_LOG_WARNING, "VBV maxrate specified, but no bufsize, ignored\n" );
        h->param.rc.i_vbv_max_bitrate = 0;
    }

    h->param.i_slice_max_size = X264_MAX( h->param.i_slice_max_size, 0 );
    h->param.i_slice_max_mbs = X264_MAX( h->param.i_slice_max_mbs, 0 );
    h->param.i_slice_min_mbs = X264_MAX( h->param.i_slice_min_mbs, 0 );
    if( h->param.i_slice_max_mbs )
        h->param.i_slice_min_mbs = X264_MIN( h->param.i_slice_min_mbs, h->param.i_slice_max_mbs/2 );
    else if( !h->param.i_slice_max_size )
        h->param.i_slice_min_mbs = 0;
    int mb_width = (h->param.i_width+15)/16;
    if( h->param.i_slice_min_mbs > mb_width )
    {
        x264_log( h, X264_LOG_WARNING, "slice-min-mbs > row mb size (%d) not implemented\n", mb_width );
        h->param.i_slice_min_mbs = mb_width;
    }

    int max_slices = (h->param.i_height+(16-1))/16;
    {
        h->param.i_slice_count = x264_clip3( h->param.i_slice_count, 0, max_slices );
        if( h->param.i_slice_max_mbs || h->param.i_slice_max_size )
            h->param.i_slice_count = 0;
    }
    if( h->param.i_slice_count_max > 0 )
        h->param.i_slice_count_max = X264_MAX( h->param.i_slice_count, h->param.i_slice_count_max );

    h->param.i_frame_reference = 1;
    h->param.i_dpb_size = x264_clip3( h->param.i_dpb_size, 1, X264_REF_MAX );
    if( h->param.i_scenecut_threshold < 0 )
        h->param.i_scenecut_threshold = 0;
    h->param.analyse.i_direct_mv_pred = x264_clip3( h->param.analyse.i_direct_mv_pred, X264_DIRECT_PRED_NONE, X264_DIRECT_PRED_AUTO );
    if( !h->param.analyse.i_subpel_refine && h->param.analyse.i_direct_mv_pred > X264_DIRECT_PRED_SPATIAL )
    {
        x264_log( h, X264_LOG_WARNING, "subme=0 + direct=temporal is not supported\n" );
        h->param.analyse.i_direct_mv_pred = X264_DIRECT_PRED_SPATIAL;
    }
    {
        h->param.analyse.i_direct_mv_pred = 0;
    }
    if( h->param.b_intra_refresh && (h->param.i_frame_reference > 1 || h->param.i_dpb_size > 1) )
    {
        x264_log( h, X264_LOG_WARNING, "ref > 1 + intra-refresh is not supported\n" );
        h->param.i_frame_reference = 1;
        h->param.i_dpb_size = 1;
    }
    if( !h->param.i_fps_num || !h->param.i_fps_den )
    {
        h->param.i_fps_num = 25;
        h->param.i_fps_den = 1;
    }
    float fps = (float) h->param.i_fps_num / h->param.i_fps_den;
    if( h->param.i_keyint_min == X264_KEYINT_MIN_AUTO )
        h->param.i_keyint_min = X264_MIN( h->param.i_keyint_max / 10, fps );
    h->param.i_keyint_min = x264_clip3( h->param.i_keyint_min, 1, h->param.i_keyint_max/2+1 );

    h->param.rc.f_qcompress = x264_clip3f( h->param.rc.f_qcompress, 0.0, 1.0 );

    h->param.i_deblocking_filter_alphac0 = x264_clip3( h->param.i_deblocking_filter_alphac0, -6, 6 );
    h->param.i_deblocking_filter_beta    = x264_clip3( h->param.i_deblocking_filter_beta, -6, 6 );
    h->param.analyse.i_luma_deadzone[0] = x264_clip3( h->param.analyse.i_luma_deadzone[0], 0, 32 );
    h->param.analyse.i_luma_deadzone[1] = x264_clip3( h->param.analyse.i_luma_deadzone[1], 0, 32 );

    if( h->param.analyse.i_me_method < X264_ME_DIA ||
        h->param.analyse.i_me_method > X264_ME_DIA )
        h->param.analyse.i_me_method = X264_ME_DIA;
    h->param.analyse.i_me_range = x264_clip3( h->param.analyse.i_me_range, 4, 1024 );
    if( h->param.analyse.i_me_range > 16 && h->param.analyse.i_me_method <= X264_ME_DIA )
        h->param.analyse.i_me_range = 16;
    h->param.analyse.inter &= X264_ANALYSE_PSUB16x16|X264_ANALYSE_PSUB8x8|X264_ANALYSE_I4x4;
    h->param.analyse.intra &= X264_ANALYSE_I4x4;
    if( !(h->param.analyse.inter & X264_ANALYSE_PSUB16x16) )
        h->param.analyse.inter &= ~X264_ANALYSE_PSUB8x8;

    h->param.rc.i_aq_mode = x264_clip3( h->param.rc.i_aq_mode, 0, 2 );
    h->param.rc.f_aq_strength = x264_clip3f( h->param.rc.f_aq_strength, 0, 3 );
    if( h->param.rc.f_aq_strength == 0 )
        h->param.rc.i_aq_mode = 0;

    if( h->param.i_log_level < X264_LOG_INFO )
    {
        h->param.analyse.b_psnr = 0;
        h->param.analyse.b_ssim = 0;
    }
    /* Warn users trying to measure PSNR/SSIM with psy opts on. */
    if( b_open && (h->param.analyse.b_psnr || h->param.analyse.b_ssim) )
    {
        char *s = NULL;

        if( h->param.analyse.b_psy )
        {
            s = h->param.analyse.b_psnr ? "psnr" : "ssim";
            x264_log( h, X264_LOG_WARNING, "--%s used with psy on: results will be invalid!\n", s );
        }
        else if( !h->param.rc.i_aq_mode && h->param.analyse.b_ssim )
        {
            x264_log( h, X264_LOG_WARNING, "--ssim used with AQ off: results will be invalid!\n" );
            s = "ssim";
        }
        else if(  h->param.rc.i_aq_mode && h->param.analyse.b_psnr )
        {
            x264_log( h, X264_LOG_WARNING, "--psnr used with AQ on: results will be invalid!\n" );
            s = "psnr";
        }
        if( s )
            x264_log( h, X264_LOG_WARNING, "--tune %s should be used if attempting to benchmark %s!\n", s, s );
    }

    if( !h->param.analyse.b_psy )
    {
        h->param.analyse.f_psy_rd = 0;
    }
    h->param.analyse.f_psy_rd = x264_clip3f( h->param.analyse.f_psy_rd, 0, 10 );
    h->mb.i_psy_rd = h->param.analyse.i_subpel_refine >= 6 ? FIX8( h->param.analyse.f_psy_rd ) : 0;
    h->param.analyse.i_chroma_qp_offset = x264_clip3(h->param.analyse.i_chroma_qp_offset, -32, 32);
    /* Psy RDO increases overall quantizers to improve the quality of luma--this indirectly hurts chroma quality */
    /* so we lower the chroma QP offset to compensate */
    if( b_open && h->mb.i_psy_rd)
        h->param.analyse.i_chroma_qp_offset -= h->param.analyse.f_psy_rd < 0.25 ? 1 : 2;
    h->param.analyse.i_chroma_qp_offset = x264_clip3(h->param.analyse.i_chroma_qp_offset, -12, 12);
    h->param.analyse.i_noise_reduction = x264_clip3( h->param.analyse.i_noise_reduction, 0, 1<<16 );

    {
        const x264_level_t *l = x264_levels;
        if( h->param.i_level_idc < 0 )
        {
            int maxrate_bak = h->param.rc.i_vbv_max_bitrate;
            if( h->param.rc.i_rc_method == X264_RC_ABR && h->param.rc.i_vbv_buffer_size <= 0 )
                h->param.rc.i_vbv_max_bitrate = h->param.rc.i_bitrate * 2;
            x264_sps_init( h->sps, h->param.i_sps_id, &h->param );
            do h->param.i_level_idc = l->level_idc;
                while( l[1].level_idc && x264_validate_levels( h, 0 ) && l++ );
            h->param.rc.i_vbv_max_bitrate = maxrate_bak;
        }
        else
        {
            while( l->level_idc && l->level_idc != h->param.i_level_idc )
                l++;
            if( l->level_idc == 0 )
            {
                x264_log( h, X264_LOG_ERROR, "invalid level_idc: %d\n", h->param.i_level_idc );
                return -1;
            }
        }
        if( h->param.analyse.i_mv_range <= 0 )
            h->param.analyse.i_mv_range = l->mv_range;
        else
            h->param.analyse.i_mv_range = x264_clip3(h->param.analyse.i_mv_range, 32, 512 );
    }

    if( h->param.rc.f_rate_tolerance < 0 )
        h->param.rc.f_rate_tolerance = 0;

    h->param.i_sps_id &= 31;

    h->param.i_nal_hrd = x264_clip3( h->param.i_nal_hrd, X264_NAL_HRD_NONE, X264_NAL_HRD_CBR );

    if( h->param.i_nal_hrd && !h->param.rc.i_vbv_buffer_size )
    {
        x264_log( h, X264_LOG_WARNING, "NAL HRD parameters require VBV parameters\n" );
        h->param.i_nal_hrd = X264_NAL_HRD_NONE;
    }

    if( h->param.i_nal_hrd == X264_NAL_HRD_CBR &&
       (h->param.rc.i_bitrate != h->param.rc.i_vbv_max_bitrate || !h->param.rc.i_vbv_max_bitrate) )
    {
        x264_log( h, X264_LOG_WARNING, "CBR HRD requires constant bitrate\n" );
        h->param.i_nal_hrd = X264_NAL_HRD_VBR;
    }

    /* ensure the booleans are 0 or 1 so they can be used in math */
#define BOOLIFY(x) h->param.x = !!h->param.x
    BOOLIFY( b_constrained_intra );
    BOOLIFY( b_deblocking_filter );
    BOOLIFY( b_intra_refresh );
    BOOLIFY( b_visualize );
    BOOLIFY( b_aud );
    BOOLIFY( b_repeat_headers );
    BOOLIFY( b_annexb );
    BOOLIFY( b_full_recon );
    BOOLIFY( analyse.b_chroma_me );
    BOOLIFY( analyse.b_fast_pskip );
    BOOLIFY( analyse.b_dct_decimate );
    BOOLIFY( analyse.b_psy );
    BOOLIFY( analyse.b_psnr );
    BOOLIFY( analyse.b_ssim );
#undef BOOLIFY

    return 0;
}

static void mbcmp_init( x264_t *h )
{
    int satd = !h->mb.b_lossless && h->param.analyse.i_subpel_refine > 1;
    memcpy( h->pixf.mbcmp, satd ? h->pixf.satd : h->pixf.sad_aligned, sizeof(h->pixf.mbcmp) );
    memcpy( h->pixf.mbcmp_unaligned, satd ? h->pixf.satd : h->pixf.sad, sizeof(h->pixf.mbcmp_unaligned) );
    h->pixf.intra_mbcmp_x3_16x16 = satd ? h->pixf.intra_satd_x3_16x16 : h->pixf.intra_sad_x3_16x16;
    h->pixf.intra_mbcmp_x3_8x16c = satd ? h->pixf.intra_satd_x3_8x16c : h->pixf.intra_sad_x3_8x16c;
    h->pixf.intra_mbcmp_x3_8x8c  = satd ? h->pixf.intra_satd_x3_8x8c  : h->pixf.intra_sad_x3_8x8c;
    h->pixf.intra_mbcmp_x3_4x4 = satd ? h->pixf.intra_satd_x3_4x4 : h->pixf.intra_sad_x3_4x4;
    memcpy( h->pixf.fpelcmp,  h->pixf.sad, sizeof(h->pixf.fpelcmp) );
    memcpy( h->pixf.fpelcmp_x4, h->pixf.sad_x4, sizeof(h->pixf.fpelcmp_x4) );
}

static void chroma_dsp_init( x264_t *h )
{
    memcpy( h->luma2chroma_pixel, x264_luma2chroma_pixel[CHROMA_FORMAT], sizeof(h->luma2chroma_pixel) );

    memcpy( h->predict_chroma, h->predict_8x8c, sizeof(h->predict_chroma) );
    h->loopf.deblock_chroma[0] = h->loopf.deblock_h_chroma_420;
    h->loopf.deblock_chroma_intra[0] = h->loopf.deblock_h_chroma_420_intra;
    h->loopf.deblock_chroma_mbaff = h->loopf.deblock_chroma_420_mbaff;
    h->loopf.deblock_chroma_intra_mbaff = h->loopf.deblock_chroma_420_intra_mbaff;
    h->pixf.intra_mbcmp_x3_chroma = h->pixf.intra_mbcmp_x3_8x8c;
    h->quantf.coeff_last[DCT_CHROMA_DC] = h->quantf.coeff_last4;
    h->quantf.coeff_level_run[DCT_CHROMA_DC] = h->quantf.coeff_level_run4;

}

static void x264_set_aspect_ratio( x264_t *h, x264_param_t *param, int initial )
{
    /* VUI */
    if( param->vui.i_sar_width > 0 && param->vui.i_sar_height > 0 )
    {
        uint32_t i_w = param->vui.i_sar_width;
        uint32_t i_h = param->vui.i_sar_height;
        uint32_t old_w = h->param.vui.i_sar_width;
        uint32_t old_h = h->param.vui.i_sar_height;

        x264_reduce_fraction( &i_w, &i_h );

        while( i_w > 65535 || i_h > 65535 )
        {
            i_w /= 2;
            i_h /= 2;
        }

        x264_reduce_fraction( &i_w, &i_h );

        if( i_w != old_w || i_h != old_h || initial )
        {
            h->param.vui.i_sar_width = 0;
            h->param.vui.i_sar_height = 0;
            if( i_w == 0 || i_h == 0 )
                x264_log( h, X264_LOG_WARNING, "cannot create valid sample aspect ratio\n" );
            else
            {
                x264_log( h, initial?X264_LOG_INFO:X264_LOG_DEBUG, "using SAR=%d/%d\n", i_w, i_h );
                h->param.vui.i_sar_width = i_w;
                h->param.vui.i_sar_height = i_h;
            }
            x264_sps_init( h->sps, h->param.i_sps_id, &h->param );
        }
    }
}

/****************************************************************************
 * x264_encoder_open:
 ****************************************************************************/
x264_t *x264_encoder_open( x264_param_t *param )
{
    x264_t *h;
    char buf[1000], *p;
    int qp;

    CHECKED_MALLOCZERO( h, sizeof(x264_t) );

    /* Create a copy of param */
    memcpy( &h->param, param, sizeof(x264_param_t) );

    if( param->param_free )
        param->param_free( param );

    if( x264_validate_parameters( h, 1 ) < 0 )
        goto fail;

    x264_reduce_fraction( &h->param.i_fps_num, &h->param.i_fps_den );

    /* Init x264_t */
    h->i_frame = -1;
    h->i_frame_num = 0;
    h->i_idr_pic_id = 0;

    x264_sps_init( h->sps, h->param.i_sps_id, &h->param );
    x264_pps_init( h->pps, h->param.i_sps_id, &h->param, h->sps );

    x264_set_aspect_ratio( h, &h->param, 1 );

    x264_validate_levels( h, 1 );

    h->chroma_qp_table = i_chroma_qp_table + 12 + h->pps->i_chroma_qp_index_offset;

    if( x264_cqm_init( h ) < 0 )
        goto fail;

    h->mb.i_mb_width = h->sps->i_mb_width;
    h->mb.i_mb_height = h->sps->i_mb_height;
    h->mb.i_mb_count = h->mb.i_mb_width * h->mb.i_mb_height;

    h->mb.chroma_h_shift = 1;
    h->mb.chroma_v_shift = 1;

    /* Init frames. */
    h->frames.i_delay = 0;

    h->frames.i_max_ref0 = 1;
    h->frames.i_max_dpb  = h->sps->vui.i_max_dec_frame_buffering;
    h->frames.b_have_lowres = ( h->param.rc.i_rc_method == X264_RC_ABR
          || h->param.rc.i_rc_method == X264_RC_CRF
          || h->param.i_scenecut_threshold );
    h->frames.b_have_sub8x8_esa = !!(h->param.analyse.inter & X264_ANALYSE_PSUB8x8);

    h->frames.i_last_idr =
    h->frames.i_last_keyframe = - h->param.i_keyint_max;
    h->frames.i_input    = 0;

    CHECKED_MALLOCZERO( h->frames.unused[0], (3) * sizeof(x264_frame_t *) );
    /* Allocate room for max refs plus a few extra just in case. */
    CHECKED_MALLOCZERO( h->frames.unused[1], (1 + X264_REF_MAX + 4) * sizeof(x264_frame_t *) );

    h->i_ref = 0;
    h->i_cpb_delay = h->i_coded_fields = h->i_disp_fields = 0;
    h->i_disp_fields_last_frame = -1;
    x264_rdo_init();

    /* init CPU functions */
    x264_predict_16x16_init(h->predict_16x16 );
    x264_predict_8x8c_init(  h->predict_8x8c );
    x264_predict_4x4_init(  h->predict_4x4 );
    x264_pixel_init(  &h->pixf );
    x264_dct_init(  &h->dctf );
    x264_zigzag_init(  &h->zigzagf_progressive, &h->zigzagf_interlaced );
    memcpy( &h->zigzagf, &h->zigzagf_progressive, sizeof(h->zigzagf) );
    x264_mc_init( &h->mc );
    x264_quant_init( h,  &h->quantf );
    x264_deblock_init(  &h->loopf );
    x264_bitstream_init(  &h->bsf );
    x264_cavlc_init( h );

    mbcmp_init( h );
    chroma_dsp_init( h );

    p = buf + sprintf( buf, "using cpu capabilities:" );
    p += sprintf( p, " none!" );
    x264_log( h, X264_LOG_INFO, "%s\n", buf );

    float *logs = x264_analyse_prepare_costs( h );
    if( !logs )
        goto fail;
    for( qp = X264_MIN( h->param.rc.i_qp_min, QP_MAX_SPEC ); qp <= h->param.rc.i_qp_max; qp++ )
        if( x264_analyse_init_costs( h, logs, qp ) )
            goto fail;
    if( x264_analyse_init_costs( h, logs, X264_LOOKAHEAD_QP ) )
        goto fail;
    x264_free( logs );

    static const uint16_t cost_mv_correct[7] = { 24, 47, 95, 189, 379, 757, 1515 };
    /* Checks for known miscompilation issues. */
    if( h->cost_mv[X264_LOOKAHEAD_QP][2013] != cost_mv_correct[BIT_DEPTH-8] )
    {
        x264_log( h, X264_LOG_ERROR, "MV cost test failed: x264 has been miscompiled!\n" );
        goto fail;
    }

    /* Must be volatile or else GCC will optimize it out. */
    volatile int temp = 392;
    if( x264_clz( temp ) != 23 )
    {
        x264_log( h, X264_LOG_ERROR, "CLZ test failed: x264 has been miscompiled!\n" );
#if ARCH_X86 || ARCH_X86_64
        x264_log( h, X264_LOG_ERROR, "Are you attempting to run an SSE4a/LZCNT-targeted build on a CPU that\n" );
        x264_log( h, X264_LOG_ERROR, "doesn't support it?\n" );
#endif
        goto fail;
    }

    h->out.i_nal = 0;
    h->out.i_bitstream = X264_MAX( 1000000, h->param.i_width * h->param.i_height * 4
        * ( h->param.rc.i_rc_method == X264_RC_ABR ? pow( 0.95, h->param.rc.i_qp_min )
          : pow( 0.95, h->param.rc.i_qp_constant ) * X264_MAX( 1, h->param.rc.f_ip_factor )));

    h->nal_buffer_size = h->out.i_bitstream * 3/2 + 4 + 64; /* +4 for startcode, +64 for nal_escape assembly padding */
    CHECKED_MALLOC( h->nal_buffer, h->nal_buffer_size );

    {
        int init_nal_count = h->param.i_slice_count + 3;

        h->fdec = x264_frame_pop_unused( h, 1 );
        if( !h->fdec )
            goto fail;

        CHECKED_MALLOC( h->out.p_bitstream, h->out.i_bitstream );
        /* Start each thread with room for init_nal_count NAL units; it'll realloc later if needed. */
        CHECKED_MALLOC( h->out.nal, init_nal_count*sizeof(x264_nal_t) );
        h->out.i_nals_allocated = init_nal_count;

        if( x264_macroblock_cache_allocate( h ) < 0 )
            goto fail;
    }


    if( x264_lookahead_init( h ) )
        goto fail;

   if( x264_macroblock_thread_allocate( h ) < 0 )
        goto fail;

    if( x264_ratecontrol_new( h ) < 0 )
        goto fail;

    if( h->param.i_nal_hrd )
    {
        x264_log( h, X264_LOG_DEBUG, "HRD bitrate: %i bits/sec\n", h->sps->vui.hrd.i_bit_rate_unscaled );
        x264_log( h, X264_LOG_DEBUG, "CPB size: %i bits\n", h->sps->vui.hrd.i_cpb_size_unscaled );
    }

    if( h->param.psz_dump_yuv )
    {
        /* create or truncate the reconstructed video file */
        FILE *f = x264_fopen( h->param.psz_dump_yuv, "w" );
        if( !f )
        {
            x264_log( h, X264_LOG_ERROR, "dump_yuv: can't write to %s\n", h->param.psz_dump_yuv );
            goto fail;
        }
        else if( !x264_is_regular_file( f ) )
        {
            x264_log( h, X264_LOG_ERROR, "dump_yuv: incompatible with non-regular file %s\n", h->param.psz_dump_yuv );
            goto fail;
        }
        fclose( f );
    }

    const char *profile = "Constrained Baseline";
    char level[4];
    snprintf( level, sizeof(level), "%d.%d", h->sps->i_level_idc/10, h->sps->i_level_idc%10 );
    if( h->sps->i_level_idc == 9 || ( h->sps->i_level_idc == 11 && h->sps->b_constraint_set3 ) )
        strcpy( level, "1b" );

    x264_log( h, X264_LOG_INFO, "profile %s, level %s\n", profile, level );

    return h;
fail:
    x264_free( h );
    return NULL;
}

/****************************************************************************
 * x264_encoder_reconfig:
 ****************************************************************************/
int x264_encoder_reconfig( x264_t *h, x264_param_t *param )
{
    int rc_reconfig = 0;
    x264_set_aspect_ratio( h, param, 0 );
#define COPY(var) h->param.var = param->var
    COPY( i_frame_reference ); // but never uses more refs than initially specified
    if( h->param.i_scenecut_threshold )
        COPY( i_scenecut_threshold ); // can't turn it on or off, only vary the threshold
    COPY( b_deblocking_filter );
    COPY( i_deblocking_filter_alphac0 );
    COPY( i_deblocking_filter_beta );
    COPY( analyse.inter );
    COPY( analyse.intra );
    COPY( analyse.i_direct_mv_pred );
    /* Scratch buffer prevents me_range from being increased for esa/tesa */
    COPY( analyse.i_me_range );
    COPY( analyse.i_noise_reduction );
    /* We can't switch out of subme=0 during encoding. */
    if( h->param.analyse.i_subpel_refine )
        COPY( analyse.i_subpel_refine );
    COPY( analyse.b_chroma_me );
    COPY( analyse.b_dct_decimate );
    COPY( analyse.b_fast_pskip );
    COPY( analyse.f_psy_rd );
    COPY( crop_rect );
    COPY( i_slice_max_size );
    COPY( i_slice_max_mbs );
    COPY( i_slice_min_mbs );
    COPY( i_slice_count );
    COPY( i_slice_count_max );

    /* VBV can't be turned on if it wasn't on to begin with */
    if( h->param.rc.i_vbv_max_bitrate > 0 && h->param.rc.i_vbv_buffer_size > 0 &&
          param->rc.i_vbv_max_bitrate > 0 &&   param->rc.i_vbv_buffer_size > 0 )
    {
        rc_reconfig |= h->param.rc.i_vbv_max_bitrate != param->rc.i_vbv_max_bitrate;
        rc_reconfig |= h->param.rc.i_vbv_buffer_size != param->rc.i_vbv_buffer_size;
        rc_reconfig |= h->param.rc.i_bitrate != param->rc.i_bitrate;
        COPY( rc.i_vbv_max_bitrate );
        COPY( rc.i_vbv_buffer_size );
        COPY( rc.i_bitrate );
    }
    rc_reconfig |= h->param.rc.f_rf_constant != param->rc.f_rf_constant;
    rc_reconfig |= h->param.rc.f_rf_constant_max != param->rc.f_rf_constant_max;
    COPY( rc.f_rf_constant );
    COPY( rc.f_rf_constant_max );
#undef COPY

    mbcmp_init( h );

    int ret = x264_validate_parameters( h, 0 );

    /* Supported reconfiguration options (1-pass only):
     * vbv-maxrate
     * vbv-bufsize
     * crf
     * bitrate (CBR only) */
    if( !ret && rc_reconfig )
        x264_ratecontrol_init_reconfigurable( h, 0 );

    return ret;
}

/****************************************************************************
 * x264_encoder_parameters:
 ****************************************************************************/
void x264_encoder_parameters( x264_t *h, x264_param_t *param )
{
    memcpy( param, &h->param, sizeof(x264_param_t) );
}

/* internal usage */
static void x264_nal_start( x264_t *h, int i_type, int i_ref_idc )
{
    x264_nal_t *nal = &h->out.nal[h->out.i_nal];

    nal->i_ref_idc        = i_ref_idc;
    nal->i_type           = i_type;
    nal->b_long_startcode = 1;

    nal->i_payload= 0;
    nal->p_payload= &h->out.p_bitstream[bs_pos( &h->out.bs ) / 8];
    nal->i_padding= 0;
}

/* if number of allocated nals is not enough, re-allocate a larger one. */
static int x264_nal_check_buffer( x264_t *h )
{
    if( h->out.i_nal >= h->out.i_nals_allocated )
    {
        x264_nal_t *new_out = x264_malloc( sizeof(x264_nal_t) * (h->out.i_nals_allocated*2) );
        if( !new_out )
            return -1;
        memcpy( new_out, h->out.nal, sizeof(x264_nal_t) * (h->out.i_nals_allocated) );
        x264_free( h->out.nal );
        h->out.nal = new_out;
        h->out.i_nals_allocated *= 2;
    }
    return 0;
}

static int x264_nal_end( x264_t *h )
{
    x264_nal_t *nal = &h->out.nal[h->out.i_nal];
    uint8_t *end = &h->out.p_bitstream[bs_pos( &h->out.bs ) / 8];
    nal->i_payload = end - nal->p_payload;
    /* Assembly implementation of nal_escape reads past the end of the input.
     * While undefined padding wouldn't actually affect the output, it makes valgrind unhappy. */
    memset( end, 0xff, 64 );
    if( h->param.nalu_process )
        h->param.nalu_process( h, nal, h->fenc->opaque );
    h->out.i_nal++;

    return x264_nal_check_buffer( h );
}

static int x264_check_encapsulated_buffer( x264_t *h, x264_t *h0, int start,
                                           int previous_nal_size, int necessary_size )
{
    if( h0->nal_buffer_size < necessary_size )
    {
        necessary_size *= 2;
        uint8_t *buf = x264_malloc( necessary_size );
        if( !buf )
            return -1;
        if( previous_nal_size )
            memcpy( buf, h0->nal_buffer, previous_nal_size );

        intptr_t delta = buf - h0->nal_buffer;
        for( int i = 0; i < start; i++ )
            h->out.nal[i].p_payload += delta;

        x264_free( h0->nal_buffer );
        h0->nal_buffer = buf;
        h0->nal_buffer_size = necessary_size;
    }

    return 0;
}

static int x264_encoder_encapsulate_nals( x264_t *h, int start )
{
    x264_t *h0 = h;
    int nal_size = 0, previous_nal_size = 0;

    if( h->param.nalu_process )
    {
        for( int i = start; i < h->out.i_nal; i++ )
            nal_size += h->out.nal[i].i_payload;
        return nal_size;
    }

    for( int i = 0; i < start; i++ )
        previous_nal_size += h->out.nal[i].i_payload;

    for( int i = start; i < h->out.i_nal; i++ )
        nal_size += h->out.nal[i].i_payload;

    /* Worst-case NAL unit escaping: reallocate the buffer if it's too small. */
    int necessary_size = previous_nal_size + nal_size * 3/2 + h->out.i_nal * 4 + 4 + 64;
    for( int i = start; i < h->out.i_nal; i++ )
        necessary_size += h->out.nal[i].i_padding;
    if( x264_check_encapsulated_buffer( h, h0, start, previous_nal_size, necessary_size ) )
        return -1;

    uint8_t *nal_buffer = h0->nal_buffer + previous_nal_size;

    for( int i = start; i < h->out.i_nal; i++ )
    {
        h->out.nal[i].b_long_startcode = !i || h->out.nal[i].i_type == NAL_SPS || h->out.nal[i].i_type == NAL_PPS ;
        x264_nal_encode( h, nal_buffer, &h->out.nal[i] );
        nal_buffer += h->out.nal[i].i_payload;
    }



    return nal_buffer - (h0->nal_buffer + previous_nal_size);
}

/****************************************************************************
 * x264_encoder_headers:
 ****************************************************************************/
int x264_encoder_headers( x264_t *h, x264_nal_t **pp_nal, int *pi_nal )
{
    int frame_size = 0;
    /* init bitstream context */
    h->out.i_nal = 0;
    bs_init( &h->out.bs, h->out.p_bitstream, h->out.i_bitstream );

    /* Write SEI, SPS and PPS. */

    /* generate sequence parameters */
    x264_nal_start( h, NAL_SPS, NAL_PRIORITY_HIGHEST );
    x264_sps_write( &h->out.bs, h->sps );
    if( x264_nal_end( h ) )
        return -1;

    /* generate picture parameters */
    x264_nal_start( h, NAL_PPS, NAL_PRIORITY_HIGHEST );
    x264_pps_write( &h->out.bs, h->sps, h->pps );
    if( x264_nal_end( h ) )
        return -1;

    /* identify ourselves */
    x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
    if( x264_sei_version_write( h, &h->out.bs ) )
        return -1;
    if( x264_nal_end( h ) )
        return -1;

    frame_size = x264_encoder_encapsulate_nals( h, 0 );
    if( frame_size < 0 )
        return -1;

    /* now set output*/
    *pi_nal = h->out.i_nal;
    *pp_nal = &h->out.nal[0];
    h->out.i_nal = 0;

    return frame_size;
}

/* Check to see whether we have chosen a reference list ordering different
 * from the standard's default. */
static  void x264_reference_check_reorder( x264_t *h )
{
    for( int i = 0; i < h->i_ref - 1; i++ )
    {
        int framenum_diff = h->fref[i+1]->i_frame_num - h->fref[i]->i_frame_num;
        int poc_diff = h->fref[i+1]->i_poc - h->fref[i]->i_poc;
        /* P and B-frames use different default orders. */
        if( h->sh.i_type == SLICE_TYPE_P ? framenum_diff > 0 : poc_diff > 0 )
        {
            h->b_ref_reorder = 1;
            return;
        }
    }
}

static  int x264_reference_distance( x264_t *h, x264_frame_t *frame )
{
   return abs(h->fenc->i_frame - frame->i_frame);
}

static  void x264_reference_build_list( x264_t *h, int i_poc )
{
    int b_ok;

    /* build ref list 0/1 */
    h->mb.pic.i_fref = h->i_ref = 0;
    if( h->sh.i_type == SLICE_TYPE_I )
        return;

    for( int i = 0; h->frames.reference[i]; i++ )
    {
        if( h->frames.reference[i]->i_poc < i_poc )
            h->fref[h->i_ref++] = h->frames.reference[i];
    }

    /* Order reference lists by distance from the current frame. */
    {
        h->fref_nearest = h->fref[0];
        do
        {
            b_ok = 1;
            for( int i = 0; i < h->i_ref - 1; i++ )
            {
                if( h->fref[i+1]->i_poc > h->fref_nearest->i_poc )
                    h->fref_nearest = h->fref[i+1];
                if( x264_reference_distance( h, h->fref[i] ) > x264_reference_distance( h, h->fref[i+1] ) )
                {
                    XCHG( x264_frame_t*, h->fref[i], h->fref[i+1] );
                    b_ok = 0;
                    break;
                }
            }
        } while( !b_ok );
    }

    if( h->sh.i_mmco_remove_from_end )
        for( int i = h->i_ref-1; i >= h->i_ref - h->sh.i_mmco_remove_from_end; i-- )
        {
            int diff = h->i_frame_num - h->fref[i]->i_frame_num;
            h->sh.mmco[h->sh.i_mmco_command_count].i_poc = h->fref[i]->i_poc;
            h->sh.mmco[h->sh.i_mmco_command_count++].i_difference_of_pic_nums = diff;
        }

    x264_reference_check_reorder( h );

    h->i_ref = X264_MIN( h->i_ref, h->frames.i_max_ref0 );
    h->i_ref = X264_MIN( h->i_ref, h->param.i_frame_reference ); // if reconfig() has lowered the limit

    /* add duplicates */
    if( h->fenc->i_type == X264_TYPE_P )
    {
        int idx = -1;
        h->mb.ref_blind_dupe = idx;
    }

    h->mb.pic.i_fref = h->i_ref;
}

static void x264_fdec_filter_row( x264_t *h, int mb_y, int pass )
{
    /* mb_y is the mb to be encoded next, not the mb to be filtered here */
    int b_hpel = h->fdec->b_kept_as_ref;
    int b_deblock = h->sh.i_disable_deblocking_filter_idc != 1;
    int b_end = mb_y == h->i_threadslice_end;
    int b_measure_quality = 1;
    int min_y = mb_y - (1);
    int b_start = min_y == h->i_threadslice_start;
    /* Even in interlaced mode, deblocking never modifies more than 4 pixels
     * above each MB, as bS=4 doesn't happen for the top of interlaced mbpairs. */
    int minpix_y = min_y*16 - 4 * !b_start;
    int maxpix_y = mb_y*16 - 4 * !b_end;
    b_deblock &= b_hpel || h->param.b_full_recon || h->param.psz_dump_yuv;
    if( mb_y & 0 )
        return;
    if( min_y < h->i_threadslice_start )
        return;

    if( b_deblock )
        for( int y = min_y; y < mb_y; y += (1) )
            x264_frame_deblock_row( h, y );

    /* FIXME: Prediction requires different borders for interlaced/progressive mc,
     * but the actual image data is equivalent. For now, maintain this
     * consistency by copying deblocked pixels between planes. */

    if( h->fdec->b_kept_as_ref )
        x264_frame_expand_border( h, h->fdec, min_y );
    if( b_hpel )
    {
        int end = mb_y == h->mb.i_mb_height;
        /* Can't do hpel until the previous slice is done encoding. */
        if( h->param.analyse.i_subpel_refine )
        {
            x264_frame_filter( h, h->fdec, min_y, end );
            x264_frame_expand_border_filtered( h, h->fdec, min_y, end );
        }
    }

    if( b_measure_quality )
    {
        maxpix_y = X264_MIN( maxpix_y, h->param.i_height );
        if( h->param.analyse.b_psnr )
        {
            for( int p = 0; p < (1); p++ )
                h->stat.frame.i_ssd[p] += x264_pixel_ssd_wxh( &h->pixf,
                    h->fdec->plane[p] + minpix_y * h->fdec->i_stride[p], h->fdec->i_stride[p],
                    h->fenc->plane[p] + minpix_y * h->fenc->i_stride[p], h->fenc->i_stride[p],
                    h->param.i_width, maxpix_y-minpix_y );
            {
                uint64_t ssd_u, ssd_v;
                int v_shift = CHROMA_V_SHIFT;
                x264_pixel_ssd_nv12( &h->pixf,
                    h->fdec->plane[1] + (minpix_y>>v_shift) * h->fdec->i_stride[1], h->fdec->i_stride[1],
                    h->fenc->plane[1] + (minpix_y>>v_shift) * h->fenc->i_stride[1], h->fenc->i_stride[1],
                    h->param.i_width>>1, (maxpix_y-minpix_y)>>v_shift, &ssd_u, &ssd_v );
                h->stat.frame.i_ssd[1] += ssd_u;
                h->stat.frame.i_ssd[2] += ssd_v;
            }
        }

        if( h->param.analyse.b_ssim )
        {
            int ssim_cnt;

            /* offset by 2 pixels to avoid alignment of ssim blocks with dct blocks,
             * and overlap by 4 */
            minpix_y += b_start ? 2 : -6;
            h->stat.frame.f_ssim +=
                x264_pixel_ssim_wxh( &h->pixf,
                    h->fdec->plane[0] + 2+minpix_y*h->fdec->i_stride[0], h->fdec->i_stride[0],
                    h->fenc->plane[0] + 2+minpix_y*h->fenc->i_stride[0], h->fenc->i_stride[0],
                    h->param.i_width-2, maxpix_y-minpix_y, h->scratch_buffer, &ssim_cnt );
            h->stat.frame.i_ssim_cnt += ssim_cnt;
        }
    }
}

static  int x264_reference_update( x264_t *h )
{
    if( !h->fdec->b_kept_as_ref )
    {
        return 0;
    }

    /* apply mmco from previous frame. */
    for( int i = 0; i < h->sh.i_mmco_command_count; i++ )
        for( int j = 0; h->frames.reference[j]; j++ )
            if( h->frames.reference[j]->i_poc == h->sh.mmco[i].i_poc )
                x264_frame_push_unused( h, x264_frame_shift( &h->frames.reference[j] ) );

    /* move frame in the buffer */
    x264_frame_push( h->frames.reference, h->fdec );
    if( h->frames.reference[h->sps->i_num_ref_frames] )
        x264_frame_push_unused( h, x264_frame_shift( h->frames.reference ) );
    h->fdec = x264_frame_pop_unused( h, 1 );
    if( !h->fdec )
        return -1;
    return 0;
}

static  void x264_reference_reset( x264_t *h )
{
    while( h->frames.reference[0] )
        x264_frame_push_unused( h, x264_frame_pop( h->frames.reference ) );
    h->fdec->i_poc =
    h->fenc->i_poc = 0;
}


static  void x264_slice_init( x264_t *h, int i_nal_type, int i_global_qp )
{
    /* ------------------------ Create slice header  ----------------------- */
    if( i_nal_type == NAL_SLICE_IDR )
    {
        x264_slice_header_init( h, &h->sh, h->sps, h->pps, h->i_idr_pic_id, h->i_frame_num, i_global_qp );

        /* alternate id */
        h->i_idr_pic_id ^= 1;
    }
    else
    {
        x264_slice_header_init( h, &h->sh, h->sps, h->pps, -1, h->i_frame_num, i_global_qp );

        h->sh.i_num_ref_idx_l0_active = h->i_ref <= 0 ? 1 : h->i_ref;
        h->sh.i_num_ref_idx_l1_active = 0;
        if( h->sh.i_num_ref_idx_l0_active != h->pps->i_num_ref_idx_l0_default_active )
        {
            h->sh.b_num_ref_idx_override = 1;
        }
    }

    h->fdec->i_frame_num = h->sh.i_frame_num;

    if( h->sps->i_poc_type == 0 )
    {
        h->sh.i_poc = h->fdec->i_poc;
        h->sh.i_delta_poc_bottom = 0;
        h->fdec->i_delta_poc[0] = h->sh.i_delta_poc_bottom == -1;
        h->fdec->i_delta_poc[1] = h->sh.i_delta_poc_bottom ==  1;
    }
    else
    {
        /* Nothing to do ? */
    }

    x264_macroblock_slice_init( h );
}

typedef struct
{
    int skip;
    uint8_t cabac_prevbyte;
    bs_t bs;
    x264_frame_stat_t stat;
    int last_qp;
    int last_dqp;
    int field_decoding_flag;
} x264_bs_bak_t;

static  void x264_bitstream_backup( x264_t *h, x264_bs_bak_t *bak, int i_skip, int full )
{
    if( full )
    {
        bak->stat = h->stat.frame;
        bak->last_qp = h->mb.i_last_qp;
        bak->last_dqp = h->mb.i_last_dqp;
        bak->field_decoding_flag = h->mb.field_decoding_flag;
    }
    else
    {
        bak->stat.i_mv_bits = h->stat.frame.i_mv_bits;
        bak->stat.i_tex_bits = h->stat.frame.i_tex_bits;
    }
    /* In the per-MB backup, we don't need the contexts because flushing the CABAC
     * encoder has no context dependency and in this case, a slice is ended (and
     * thus the content of all contexts are thrown away). */
    {
        bak->bs = h->out.bs;
        bak->skip = i_skip;
    }
}

static  void x264_bitstream_restore( x264_t *h, x264_bs_bak_t *bak, int *skip, int full )
{
    if( full )
    {
        h->stat.frame = bak->stat;
        h->mb.i_last_qp = bak->last_qp;
        h->mb.i_last_dqp = bak->last_dqp;
        h->mb.field_decoding_flag = bak->field_decoding_flag;
    }
    else
    {
        h->stat.frame.i_mv_bits = bak->stat.i_mv_bits;
        h->stat.frame.i_tex_bits = bak->stat.i_tex_bits;
    }
    {
        h->out.bs = bak->bs;
        *skip = bak->skip;
    }
}

static int x264_slice_write( x264_t *h )
{
    int i_skip;
    int mb_xy, i_mb_x, i_mb_y;
    /* NALUs other than the first use a 3-byte startcode.
     * Add one extra byte for the rbsp, and one more for the final CABAC putbyte.
     * Then add an extra 5 bytes just in case, to account for random NAL escapes and
     * other inaccuracies. */
    int overhead_guess = (NALU_OVERHEAD - (h->param.b_annexb && h->out.i_nal)) + 1 + 5;
    int slice_max_size = h->param.i_slice_max_size > 0 ? (h->param.i_slice_max_size-overhead_guess)*8 : 0;
    int back_up_bitstream = slice_max_size || (h->sps->i_profile_idc < PROFILE_HIGH);
    int starting_bits = bs_pos(&h->out.bs);
    int b_deblock = h->sh.i_disable_deblocking_filter_idc != 1;
    int b_hpel = h->fdec->b_kept_as_ref;
    int orig_last_mb = h->sh.i_last_mb;
    int thread_last_mb = h->i_threadslice_end * h->mb.i_mb_width - 1;
    uint8_t *last_emu_check;
#define BS_BAK_SLICE_MAX_SIZE 0
#define BS_BAK_SLICE_MIN_MBS  1
#define BS_BAK_ROW_VBV        2
    x264_bs_bak_t bs_bak[3];
    b_deblock &= b_hpel || h->param.b_full_recon || h->param.psz_dump_yuv;
    bs_realign( &h->out.bs );

    /* Slice */
    x264_nal_start( h, h->i_nal_type, h->i_nal_ref_idc );
    h->out.nal[h->out.i_nal].i_first_mb = h->sh.i_first_mb;

    /* Slice header */
    x264_macroblock_thread_init( h );

    /* Set the QP equal to the first QP in the slice for more accurate CABAC initialization. */
    h->mb.i_mb_xy = h->sh.i_first_mb;
    h->sh.i_qp = x264_ratecontrol_mb_qp( h );
    h->sh.i_qp = SPEC_QP( h->sh.i_qp );
    h->sh.i_qp_delta = h->sh.i_qp - h->pps->i_pic_init_qp;

    x264_slice_header_write( &h->out.bs, &h->sh, h->i_nal_ref_idc );

    last_emu_check = h->out.bs.p;
    h->mb.i_last_qp = h->sh.i_qp;
    h->mb.i_last_dqp = 0;
    h->mb.field_decoding_flag = 0;

    i_mb_y = h->sh.i_first_mb / h->mb.i_mb_width;
    i_mb_x = h->sh.i_first_mb % h->mb.i_mb_width;
    i_skip = 0;

    while( 1 )
    {
        mb_xy = i_mb_x + i_mb_y * h->mb.i_mb_width;
        int mb_spos = bs_pos(&h->out.bs) ;

        if( i_mb_x == 0 )
        {
            if( x264_bitstream_check_buffer( h ) )
                return -1;
            if( h->param.rc.i_vbv_buffer_size )
                x264_bitstream_backup( h, &bs_bak[BS_BAK_ROW_VBV], i_skip, 1 );
            if( !h->mb.b_reencode_mb )
                x264_fdec_filter_row( h, i_mb_y, 0 );
        }

        if( !(i_mb_y & 0) && back_up_bitstream )
        {
            x264_bitstream_backup( h, &bs_bak[BS_BAK_SLICE_MAX_SIZE], i_skip, 0 );
            if( slice_max_size && (thread_last_mb+1-mb_xy) == h->param.i_slice_min_mbs )
                x264_bitstream_backup( h, &bs_bak[BS_BAK_SLICE_MIN_MBS], i_skip, 0 );
        }

        /* load cache */
        x264_macroblock_cache_load_progressive( h, i_mb_x, i_mb_y );

        x264_macroblock_analyse( h );

        /* encode this macroblock -> be careful it can change the mb type to P_SKIP if needed */
reencode:
        x264_macroblock_encode( h );

        {
            if( IS_SKIP( h->mb.i_type ) )
                i_skip++;
            else
            {
                if( h->sh.i_type != SLICE_TYPE_I )
                {
                    bs_write_ue( &h->out.bs, i_skip );  /* skip run */
                    i_skip = 0;
                }
                x264_macroblock_write_cavlc( h );
                /* If there was a CAVLC level code overflow, try again at a higher QP. */
                if( h->mb.b_overflow )
                {
                    h->mb.i_chroma_qp = h->chroma_qp_table[++h->mb.i_qp];
                    h->mb.i_skip_intra = 0;
                    h->mb.b_skip_mc = 0;
                    h->mb.b_overflow = 0;
                    x264_bitstream_restore( h, &bs_bak[BS_BAK_SLICE_MAX_SIZE], &i_skip, 0 );
                    goto reencode;
                }
            }
        }

        int total_bits = bs_pos(&h->out.bs);
        int mb_size = total_bits - mb_spos;

        if( slice_max_size )
        {
            /* Count the skip run, just in case. */
            total_bits += bs_size_ue_big( i_skip );
            /* Check for escape bytes. */
            uint8_t *end = h->out.bs.p;
            for( ; last_emu_check < end - 2; last_emu_check++ )
                if( last_emu_check[0] == 0 && last_emu_check[1] == 0 && last_emu_check[2] <= 3 )
                {
                    slice_max_size -= 8;
                    last_emu_check++;
                }
            /* We'll just re-encode this last macroblock if we go over the max slice size. */
            if( total_bits - starting_bits > slice_max_size && !h->mb.b_reencode_mb )
            {
                if( !x264_frame_new_slice( h, h->fdec ) )
                {
                    /* Handle the most obnoxious slice-min-mbs edge case: we need to end the slice
                     * because it's gone over the maximum size, but doing so would violate slice-min-mbs.
                     * If possible, roll back to the last checkpoint and try again.
                     * We could try raising QP, but that would break in the case where a slice spans multiple
                     * rows, which the re-encoding infrastructure can't currently handle. */
                    if( mb_xy <= thread_last_mb && (thread_last_mb+1-mb_xy) < h->param.i_slice_min_mbs )
                    {
                        if( thread_last_mb-h->param.i_slice_min_mbs < h->sh.i_first_mb+h->param.i_slice_min_mbs )
                        {
                            x264_log( h, X264_LOG_WARNING, "slice-max-size violated (frame %d, cause: slice-min-mbs)\n", h->i_frame );
                            slice_max_size = 0;
                            goto cont;
                        }
                        x264_bitstream_restore( h, &bs_bak[BS_BAK_SLICE_MIN_MBS], &i_skip, 0 );
                        h->mb.b_reencode_mb = 1;
                        h->sh.i_last_mb = thread_last_mb-h->param.i_slice_min_mbs;
                        break;
                    }
                    if( mb_xy != h->sh.i_first_mb )
                    {
                        x264_bitstream_restore( h, &bs_bak[BS_BAK_SLICE_MAX_SIZE], &i_skip, 0 );
                        h->mb.b_reencode_mb = 1;
                        h->sh.i_last_mb = mb_xy-1;
                        break;
                    }
                    else
                        h->sh.i_last_mb = mb_xy;
                }
                else
                    slice_max_size = 0;
            }
        }
cont:
        h->mb.b_reencode_mb = 0;

#if HAVE_VISUALIZE
        if( h->param.b_visualize )
            x264_visualize_mb( h );
#endif

        /* save cache */
        x264_macroblock_cache_save( h );

        if( x264_ratecontrol_mb( h, mb_size ) < 0 )
        {
            x264_bitstream_restore( h, &bs_bak[BS_BAK_ROW_VBV], &i_skip, 1 );
            h->mb.b_reencode_mb = 1;
            i_mb_x = 0;
            h->mb.i_mb_prev_xy = i_mb_y * h->mb.i_mb_stride - 1;
            h->sh.i_last_mb = orig_last_mb;
            continue;
        }

        /* accumulate mb stats */
        h->stat.frame.i_mb_count[h->mb.i_type]++;

        int b_intra = IS_INTRA( h->mb.i_type );
        int b_skip = IS_SKIP( h->mb.i_type );
        if( h->param.i_log_level >= X264_LOG_INFO )
        {
            if( !b_intra && !b_skip )
            {
                if( h->mb.i_partition != D_8x8 )
                        h->stat.frame.i_mb_partition[h->mb.i_partition] += 4;
                    else
                        for( int i = 0; i < 4; i++ )
                            h->stat.frame.i_mb_partition[h->mb.i_sub_partition[i]] ++;
            }
        }

        if( h->param.i_log_level >= X264_LOG_INFO )
        {
            if( h->mb.i_cbp_luma | h->mb.i_cbp_chroma )
            {
                {
                    int cbpsum = (h->mb.i_cbp_luma&1) + ((h->mb.i_cbp_luma>>1)&1)
                               + ((h->mb.i_cbp_luma>>2)&1) + (h->mb.i_cbp_luma>>3);
                    h->stat.frame.i_mb_cbp[!b_intra + 0] += cbpsum;
                    h->stat.frame.i_mb_cbp[!b_intra + 2] += !!h->mb.i_cbp_chroma;
                    h->stat.frame.i_mb_cbp[!b_intra + 4] += h->mb.i_cbp_chroma >> 1;
                }
            }
            if( b_intra && h->mb.i_type != I_PCM )
            {
                if( h->mb.i_type == I_16x16 )
                    h->stat.frame.i_mb_pred_mode[0][h->mb.i_intra16x16_pred_mode]++;
                else //if( h->mb.i_type == I_4x4 )
                    for( int i = 0; i < 16; i++ )
                        h->stat.frame.i_mb_pred_mode[2][h->mb.cache.intra4x4_pred_mode[x264_scan8[i]]]++;
                h->stat.frame.i_mb_pred_mode[3][x264_mb_chroma_pred_mode_fix[h->mb.i_chroma_pred_mode]]++;
            }
        }

        /* calculate deblock strength values (actual deblocking is done per-row along with hpel) */
        if( b_deblock )
            x264_macroblock_deblock_strength( h );

        if( mb_xy == h->sh.i_last_mb )
            break;

        i_mb_x++;
        if( i_mb_x == h->mb.i_mb_width )
        {
            i_mb_y++;
            i_mb_x = 0;
        }
    }
    if( h->sh.i_last_mb < h->sh.i_first_mb )
        return 0;

    h->out.nal[h->out.i_nal].i_last_mb = h->sh.i_last_mb;

    {
        if( i_skip > 0 )
            bs_write_ue( &h->out.bs, i_skip );  /* last skip run */
        /* rbsp_slice_trailing_bits */
        bs_rbsp_trailing( &h->out.bs );
        bs_flush( &h->out.bs );
    }
    if( x264_nal_end( h ) )
        return -1;

    if( h->sh.i_last_mb == (h->i_threadslice_end * h->mb.i_mb_width - 1) )
    {
        h->stat.frame.i_misc_bits = bs_pos( &h->out.bs )
                                  + (h->out.i_nal*NALU_OVERHEAD * 8)
                                  - h->stat.frame.i_tex_bits
                                  - h->stat.frame.i_mv_bits;
        x264_fdec_filter_row( h, h->i_threadslice_end, 0 );

    }

    return 0;
}

static void *x264_slices_write( x264_t *h )
{
    int i_slice_num = 0;
    //last_thread_mb means the last mb of current frame
    //h->sh.i_last_mb means current slice las mb
    int last_thread_mb = h->sh.i_last_mb;

#if HAVE_VISUALIZE
    if( h->param.b_visualize )
        if( x264_visualize_init( h ) )
            goto fail;
#endif

    /* init stats */
    memset( &h->stat.frame, 0, sizeof(h->stat.frame) );
    h->mb.b_reencode_mb = 0;
    while( h->sh.i_first_mb <= last_thread_mb )
    {
        //if is the first slice or new slice success
        if( !i_slice_num || !x264_frame_new_slice( h, h->fdec ) )
        {
            if( h->param.i_slice_max_mbs )
            {
                h->sh.i_last_mb = h->sh.i_first_mb + h->param.i_slice_max_mbs - 1;
                if( h->sh.i_last_mb < last_thread_mb && last_thread_mb - h->sh.i_last_mb < h->param.i_slice_min_mbs )
                    h->sh.i_last_mb = last_thread_mb - h->param.i_slice_min_mbs;
            }
            else if( h->param.i_slice_count )
            {
                h->sh.i_last_mb = (h->mb.i_mb_height * i_slice_num + h->param.i_slice_count/2) / h->param.i_slice_count * h->mb.i_mb_width - 1;
            }
            i_slice_num++;
        }
        h->sh.i_last_mb = X264_MIN( h->sh.i_last_mb, last_thread_mb );
        if( x264_slice_write( h ) )
            goto fail;
        h->sh.i_first_mb = h->sh.i_last_mb + 1;
    }

#if HAVE_VISUALIZE
    if( h->param.b_visualize )
    {
        x264_visualize_show( h );
        x264_visualize_close( h );
    }
#endif

    return (void *)0;

fail:
    return (void *)-1;
}

void x264_encoder_intra_refresh( x264_t *h )
{
    h->b_queued_intra_refresh = 1;
}

int x264_encoder_invalidate_reference( x264_t *h, int64_t pts )
{
    if( h->param.b_intra_refresh )
    {
        x264_log( h, X264_LOG_ERROR, "x264_encoder_invalidate_reference is not supported with intra refresh enabled\n" );
        return -1;
    }

    return 0;
}

/****************************************************************************
 * x264_encoder_encode:
 *  XXX: i_poc   : is the poc of the current given picture
 *       i_frame : is the number of the frame being coded
 *  ex:  type frame poc
 *       I      0   2*0
 *       P      1   2*3
 *       B      2   2*1
 *       B      3   2*2
 *       P      4   2*6
 *       B      5   2*4
 *       B      6   2*5
 ****************************************************************************/
int     x264_encoder_encode( x264_t *h,
                             x264_nal_t **pp_nal, int *pi_nal,
                             x264_picture_t *pic_in,
                             x264_picture_t *pic_out )
{
    int i_nal_type, i_nal_ref_idc, i_global_qp;
    int overhead = NALU_OVERHEAD;

    h->i_cpb_delay_pir_offset = h->i_cpb_delay_pir_offset_next;

    /* no data out */
    *pi_nal = 0;
    *pp_nal = NULL;

    /* ------------------- Setup new frame from picture -------------------- */

    /* 1: Copy the picture to a frame and move it to a buffer */
    x264_frame_t *fenc = x264_frame_pop_unused( h, 0 );
    if( !fenc )
        return -1;

    if( x264_frame_copy_picture( h, fenc, pic_in ) < 0 )
        return -1;

    if( h->param.i_width != 16 * h->mb.i_mb_width ||
            h->param.i_height != 16 * h->mb.i_mb_height )
        x264_frame_expand_border_mod16( h, fenc );

    fenc->i_frame = h->frames.i_input++;

    x264_adaptive_quant_frame( h, fenc );

    if( h->frames.b_have_lowres )
        x264_frame_init_lowres( h, fenc );

    /* 2: Place the frame into the queue for its slice type decision */
    h->lookahead->current_frame = fenc;

    assert( h->frames.i_input > h->frames.i_delay );


    h->i_frame++;
    /* 3: The picture is analyzed in the lookahead */
    if( !h->frames.current )
        x264_lookahead_get_frame( h );

    //now we has the slice type

    /* ------------------- Get frame to be encoded ------------------------- */
    /* 4: get picture to encode */
    h->fenc = h->frames.current;
    h->frames.current = 0;

    if( h->fenc->param )
    {
        x264_encoder_reconfig( h, h->fenc->param );
        if( h->fenc->param->param_free )
        {
            h->fenc->param->param_free( h->fenc->param );
            h->fenc->param = NULL;
        }
    }

    // ok to call this before encoding any frames, since the initial values of fdec have b_kept_as_ref=0
    if( x264_reference_update( h ) )
        return -1;
    h->fdec->i_lines_completed = -1;

    if( h->fenc->b_keyframe )
    {
        h->frames.i_last_keyframe = h->fenc->i_frame;
        if( h->fenc->i_type == X264_TYPE_IDR )
        {
            h->i_frame_num = 0;
            h->frames.i_last_idr = h->fenc->i_frame;
        }
    }
    h->sh.i_mmco_command_count =
    h->sh.i_mmco_remove_from_end = 0;
    h->b_ref_reorder = 0;
    h->fdec->i_poc =
    h->fenc->i_poc = 2 * ( h->fenc->i_frame - X264_MAX( h->frames.i_last_idr, 0 ) );

    /* ------------------- Setup frame context ----------------------------- */
    /* 5: Init data dependent of frame type */
    if( h->fenc->i_type == X264_TYPE_IDR )
    {
        /* reset ref pictures */
        i_nal_type    = NAL_SLICE_IDR;
        i_nal_ref_idc = NAL_PRIORITY_HIGHEST;
        h->sh.i_type = SLICE_TYPE_I;
        x264_reference_reset( h );
    }
    else if( h->fenc->i_type == X264_TYPE_I )
    {
        i_nal_type    = NAL_SLICE;
        i_nal_ref_idc = NAL_PRIORITY_HIGH; /* Not completely true but for now it is (as all I/P are kept as ref)*/
        h->sh.i_type = SLICE_TYPE_I;
    }
    else if( h->fenc->i_type == X264_TYPE_P )
    {
        i_nal_type    = NAL_SLICE;
        i_nal_ref_idc = NAL_PRIORITY_HIGH; /* Not completely true but for now it is (as all I/P are kept as ref)*/
        h->sh.i_type = SLICE_TYPE_P;
    }

    h->fdec->i_type = h->fenc->i_type;
    h->fdec->i_frame = h->fenc->i_frame;
    h->fenc->b_kept_as_ref =
    h->fdec->b_kept_as_ref = i_nal_ref_idc != NAL_PRIORITY_DISPOSABLE && h->param.i_keyint_max > 1;

    /* ------------------- Init                ----------------------------- */
    /* build ref list 0/1 */
    x264_reference_build_list( h, h->fdec->i_poc );

    /* ---------------------- Write the bitstream -------------------------- */
    /* Init bitstream context */
    {
        bs_init( &h->out.bs, h->out.p_bitstream, h->out.i_bitstream );
        h->out.i_nal = 0;
    }

    if( h->param.b_aud )
    {
        int pic_type;

        if( h->sh.i_type == SLICE_TYPE_I )
            pic_type = 0;
        else if( h->sh.i_type == SLICE_TYPE_P )
            pic_type = 1;
        else
            pic_type = 7;

        x264_nal_start( h, NAL_AUD, NAL_PRIORITY_DISPOSABLE );
        bs_write( &h->out.bs, 3, pic_type );
        bs_rbsp_trailing( &h->out.bs );
        if( x264_nal_end( h ) )
            return -1;
        overhead += h->out.nal[h->out.i_nal-1].i_payload + NALU_OVERHEAD;
    }

    h->i_nal_type = i_nal_type;
    h->i_nal_ref_idc = i_nal_ref_idc;

    if( h->param.b_intra_refresh )
    {
        if( IS_X264_TYPE_I( h->fenc->i_type ) )
        {
            h->fdec->i_frames_since_pir = 0;
            h->b_queued_intra_refresh = 0;
            /* PIR is currently only supported with ref == 1, so any intra frame effectively refreshes
             * the whole frame and counts as an intra refresh. */
            h->fdec->f_pir_position = h->mb.i_mb_width;
        }
        else if( h->fenc->i_type == X264_TYPE_P )
        {
            int pocdiff = (h->fdec->i_poc - h->fref[0]->i_poc)/2;
            float increment = X264_MAX( ((float)h->mb.i_mb_width-1) / h->param.i_keyint_max, 1 );
            h->fdec->f_pir_position = h->fref[0]->f_pir_position;
            h->fdec->i_frames_since_pir = h->fref[0]->i_frames_since_pir + pocdiff;
            if( h->fdec->i_frames_since_pir >= h->param.i_keyint_max ||
                (h->b_queued_intra_refresh && h->fdec->f_pir_position + 0.5 >= h->mb.i_mb_width) )
            {
                h->fdec->f_pir_position = 0;
                h->fdec->i_frames_since_pir = 0;
                h->b_queued_intra_refresh = 0;
                h->fenc->b_keyframe = 1;
            }
            h->fdec->i_pir_start_col = h->fdec->f_pir_position+0.5;
            h->fdec->f_pir_position += increment * pocdiff;
            h->fdec->i_pir_end_col = h->fdec->f_pir_position+0.5;
            /* If our intra refresh has reached the right side of the frame, we're done. */
            if( h->fdec->i_pir_end_col >= h->mb.i_mb_width - 1 )
            {
                h->fdec->f_pir_position = h->mb.i_mb_width;
                h->fdec->i_pir_end_col = h->mb.i_mb_width - 1;
            }
        }
    }

    if( h->fenc->b_keyframe )
    {
        /* Write SPS and PPS */
        if( h->param.b_repeat_headers )
        {
            /* generate sequence parameters */
            x264_nal_start( h, NAL_SPS, NAL_PRIORITY_HIGHEST );
            x264_sps_write( &h->out.bs, h->sps );
            if( x264_nal_end( h ) )
                return -1;
            overhead += h->out.nal[h->out.i_nal-1].i_payload + h->out.nal[h->out.i_nal-1].i_padding + NALU_OVERHEAD;

            /* generate picture parameters */
            x264_nal_start( h, NAL_PPS, NAL_PRIORITY_HIGHEST );
            x264_pps_write( &h->out.bs, h->sps, h->pps );
            if( x264_nal_end( h ) )
                return -1;
            overhead += h->out.nal[h->out.i_nal-1].i_payload + h->out.nal[h->out.i_nal-1].i_padding + NALU_OVERHEAD;
        }

        /* when frame threading is used, buffering period sei is written in x264_encoder_frame_end */
        if(h->sps->vui.b_nal_hrd_parameters_present )
        {
            x264_hrd_fullness( h );
            x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
            x264_sei_buffering_period_write( h, &h->out.bs );
            if( x264_nal_end( h ) )
               return -1;
            overhead += h->out.nal[h->out.i_nal-1].i_payload + SEI_OVERHEAD;
        }
    }

    /* write extra sei */
    for( int i = 0; i < h->fenc->extra_sei.num_payloads; i++ )
    {
        x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
        x264_sei_write( &h->out.bs, h->fenc->extra_sei.payloads[i].payload, h->fenc->extra_sei.payloads[i].payload_size,
                        h->fenc->extra_sei.payloads[i].payload_type );
        if( x264_nal_end( h ) )
            return -1;
        overhead += h->out.nal[h->out.i_nal-1].i_payload + SEI_OVERHEAD;
        if( h->fenc->extra_sei.sei_free )
        {
            h->fenc->extra_sei.sei_free( h->fenc->extra_sei.payloads[i].payload );
            h->fenc->extra_sei.payloads[i].payload = NULL;
        }
    }

    if( h->fenc->extra_sei.sei_free )
    {
        h->fenc->extra_sei.sei_free( h->fenc->extra_sei.payloads );
        h->fenc->extra_sei.payloads = NULL;
        h->fenc->extra_sei.sei_free = NULL;
    }

    if( h->fenc->b_keyframe )
    {
        /* Avid's decoder strictly wants two SEIs for AVC-Intra so we can't insert the x264 SEI */
        if( h->param.b_repeat_headers && h->fenc->i_frame == 0 )
        {
            /* identify ourself */
            x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
            if( x264_sei_version_write( h, &h->out.bs ) )
                return -1;
            if( x264_nal_end( h ) )
                return -1;
            overhead += h->out.nal[h->out.i_nal-1].i_payload + SEI_OVERHEAD;
        }

        if( h->fenc->i_type != X264_TYPE_IDR )
        {
            int time_to_recovery = X264_MIN( h->mb.i_mb_width - 1, h->param.i_keyint_max ) - 1;
            x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
            x264_sei_recovery_point_write( h, &h->out.bs, time_to_recovery );
            if( x264_nal_end( h ) )
                return -1;
            overhead += h->out.nal[h->out.i_nal-1].i_payload + SEI_OVERHEAD;
        }

    }

    /* generate sei pic timing */
    if( h->sps->vui.b_nal_hrd_parameters_present )
    {
        x264_nal_start( h, NAL_SEI, NAL_PRIORITY_DISPOSABLE );
        x264_sei_pic_timing_write( h, &h->out.bs );
        if( x264_nal_end( h ) )
            return -1;
        overhead += h->out.nal[h->out.i_nal-1].i_payload + SEI_OVERHEAD;
    }

    if( h->fenc->b_keyframe && h->param.b_intra_refresh )
        h->i_cpb_delay_pir_offset_next = h->fenc->i_cpb_delay;

    /* Init the rate control */
    /* FIXME: Include slice header bit cost. */
    x264_ratecontrol_start( h, h->fenc->i_qpplus1, overhead*8 );
    i_global_qp = x264_ratecontrol_qp( h );

    pic_out->i_qpplus1 =
    h->fdec->i_qpplus1 = i_global_qp + 1;

    if( h->i_ref )
        h->fdec->i_poc_l0ref0 = h->fref[0]->i_poc;

    /* ------------------------ Create slice header  ----------------------- */
    x264_slice_init( h, i_nal_type, i_global_qp );

    if( i_nal_ref_idc != NAL_PRIORITY_DISPOSABLE )
        h->i_frame_num++;

    /* Write frame */
    h->i_threadslice_start = 0;
    h->i_threadslice_end = h->mb.i_mb_height;

    if( (intptr_t)x264_slices_write( h ) )
        return -1;

    return x264_encoder_frame_end( h, pp_nal, pi_nal, pic_out );
}

static int x264_encoder_frame_end( x264_t *h,
                                   x264_nal_t **pp_nal, int *pi_nal,
                                   x264_picture_t *pic_out )
{
    char psz_message[80];

    if( !h->out.i_nal )
    {
        pic_out->i_type = X264_TYPE_AUTO;
        return 0;
    }



    int frame_size = x264_encoder_encapsulate_nals( h, 0 );
    if( frame_size < 0 )
        return -1;

    /* Set output picture properties */
    pic_out->i_type = h->fenc->i_type;

    pic_out->b_keyframe = h->fenc->b_keyframe;

    pic_out->opaque = h->fenc->opaque;

    pic_out->img.i_csp = h->fdec->i_csp;
    pic_out->img.i_plane = h->fdec->i_plane;
    for( int i = 0; i < pic_out->img.i_plane; i++ )
    {
        pic_out->img.i_stride[i] = h->fdec->i_stride[i] * sizeof(pixel);
        pic_out->img.plane[i] = (uint8_t*)h->fdec->plane[i];
    }

    x264_frame_push_unused( h, h->fenc );

    /* ---------------------- Update encoder state ------------------------- */

    /* update rc */
    int filler = 0;
    if( x264_ratecontrol_end( h, frame_size * 8, &filler ) < 0 )
        return -1;

    pic_out->hrd_timing = h->fenc->hrd_timing;
    pic_out->prop.f_crf_avg = h->fdec->f_crf_avg;

    {
        while( filler > 0 )
        {
            int f, overhead;
            overhead = (FILLER_OVERHEAD - h->param.b_annexb);
            if( h->param.i_slice_max_size && filler > h->param.i_slice_max_size )
            {
                int next_size = filler - h->param.i_slice_max_size;
                int overflow = X264_MAX( overhead - next_size, 0 );
                f = h->param.i_slice_max_size - overhead - overflow;
            }
            else
                f = X264_MAX( 0, filler - overhead );

            if( x264_bitstream_check_buffer_filler( h, f ) )
                return -1;
            x264_nal_start( h, NAL_FILLER, NAL_PRIORITY_DISPOSABLE );
            x264_filler_write( h, &h->out.bs, f );
            if( x264_nal_end( h ) )
                return -1;
            int total_size = x264_encoder_encapsulate_nals( h, h->out.i_nal-1 );
            if( total_size < 0 )
                return -1;
            frame_size += total_size;
            filler -= total_size;
        }
    }

    /* End bitstream, set output  */
    *pi_nal = h->out.i_nal;
    *pp_nal = h->out.nal;

    h->out.i_nal = 0;

    x264_noise_reduction_update( h );

    /* ---------------------- Compute/Print statistics --------------------- */

    /* Slice stat */
    h->stat.i_frame_count[h->sh.i_type]++;
    h->stat.i_frame_size[h->sh.i_type] += frame_size;
    h->stat.f_frame_qp[h->sh.i_type] += h->fdec->f_qp_avg_aq;

    for( int i = 0; i < X264_MBTYPE_MAX; i++ )
        h->stat.i_mb_count[h->sh.i_type][i] += h->stat.frame.i_mb_count[i];
    for( int i = 0; i < X264_PARTTYPE_MAX; i++ )
        h->stat.i_mb_partition[h->sh.i_type][i] += h->stat.frame.i_mb_partition[i];
    for( int i = 0; i < 6; i++ )
        h->stat.i_mb_cbp[i] += h->stat.frame.i_mb_cbp[i];
    for( int i = 0; i < 4; i++ )
        for( int j = 0; j < 13; j++ )
            h->stat.i_mb_pred_mode[i][j] += h->stat.frame.i_mb_pred_mode[i][j];
    if( h->sh.i_type != SLICE_TYPE_I )
            for( int i = 0; i < X264_REF_MAX*2; i++ )
                h->stat.i_mb_count_ref[h->sh.i_type][i] += h->stat.frame.i_mb_count_ref[i];
    for( int i = 0; i < 3; i++ )
        h->stat.i_mb_field[i] += h->stat.frame.i_mb_field[i];

    psz_message[0] = '\0';
    double dur = h->fenc->f_duration;
    h->stat.f_frame_duration[h->sh.i_type] += dur;
    if( h->param.analyse.b_psnr )
    {
        int64_t ssd[3] =
        {
            h->stat.frame.i_ssd[0],
            h->stat.frame.i_ssd[1],
            h->stat.frame.i_ssd[2],
        };
        int luma_size = h->param.i_width * h->param.i_height;
        int chroma_size = CHROMA_SIZE( luma_size );
        pic_out->prop.f_psnr[0] = x264_psnr( ssd[0], luma_size );
        pic_out->prop.f_psnr[1] = x264_psnr( ssd[1], chroma_size );
        pic_out->prop.f_psnr[2] = x264_psnr( ssd[2], chroma_size );
        pic_out->prop.f_psnr_avg = x264_psnr( ssd[0] + ssd[1] + ssd[2], luma_size + chroma_size*2 );

        h->stat.f_ssd_global[h->sh.i_type]   += dur * (ssd[0] + ssd[1] + ssd[2]);
        h->stat.f_psnr_average[h->sh.i_type] += dur * pic_out->prop.f_psnr_avg;
        h->stat.f_psnr_mean_y[h->sh.i_type]  += dur * pic_out->prop.f_psnr[0];
        h->stat.f_psnr_mean_u[h->sh.i_type]  += dur * pic_out->prop.f_psnr[1];
        h->stat.f_psnr_mean_v[h->sh.i_type]  += dur * pic_out->prop.f_psnr[2];

        snprintf( psz_message, 80, " PSNR Y:%5.2f U:%5.2f V:%5.2f", pic_out->prop.f_psnr[0],
                                                                    pic_out->prop.f_psnr[1],
                                                                    pic_out->prop.f_psnr[2] );
    }

    if( h->param.analyse.b_ssim )
    {
        pic_out->prop.f_ssim = h->stat.frame.f_ssim / h->stat.frame.i_ssim_cnt;
        h->stat.f_ssim_mean_y[h->sh.i_type] += pic_out->prop.f_ssim * dur;
        snprintf( psz_message + strlen(psz_message), 80 - strlen(psz_message),
                  " SSIM Y:%.5f", pic_out->prop.f_ssim );
    }
    psz_message[79] = '\0';

    x264_log( h, X264_LOG_DEBUG,
                  "frame=%4d QP=%.2f NAL=%d Slice:%c Poc:%-3d I:%-4d P:%-4d SKIP:%-4d size=%d bytes%s\n",
              h->i_frame,
              h->fdec->f_qp_avg_aq,
              h->i_nal_ref_idc,
              h->sh.i_type == SLICE_TYPE_I ? 'I' : (h->sh.i_type == SLICE_TYPE_P ? 'P' : 'B' ),
              h->fdec->i_poc,
              h->stat.frame.i_mb_count_i,
              h->stat.frame.i_mb_count_p,
              h->stat.frame.i_mb_count_skip,
              frame_size,
              psz_message );

#ifdef DEBUG_MB_TYPE
{
    static const char mb_chars[] = { 'i', 'i', 'I', 'C', 'P', '8', 'S',
        'D', '<', 'X', 'B', 'X', '>', 'B', 'B', 'B', 'B', '8', 'S' };
    for( int mb_xy = 0; mb_xy < h->mb.i_mb_width * h->mb.i_mb_height; mb_xy++ )
    {
        if( h->mb.type[mb_xy] < X264_MBTYPE_MAX && h->mb.type[mb_xy] >= 0 )
            fprintf( stderr, "%c ", mb_chars[ h->mb.type[mb_xy] ] );
        else
            fprintf( stderr, "? " );

        if( (mb_xy+1) % h->mb.i_mb_width == 0 )
            fprintf( stderr, "\n" );
    }
}
#endif

    /* Remove duplicates, must be done near the end as breaks h->fref0 array
     * by freeing some of its pointers. */
    for( int i = 0; i < h->i_ref; i++ )
        if( h->fref[i] && h->fref[i]->b_duplicate )
        {
            x264_frame_push_blank_unused( h, h->fref[i] );
            h->fref[i] = 0;
        }

    if( h->param.psz_dump_yuv )
        x264_frame_dump( h );


    return frame_size;
}

static void x264_print_intra( int64_t *i_mb_count, double i_count, int b_print_pcm, char *intra )
{
    intra += sprintf( intra, "I16..4%s: %4.1f%% %4.1f%% %4.1f%%",
        b_print_pcm ? "..PCM" : "",
        i_mb_count[I_16x16]/ i_count,
        i_mb_count[I_8x8]  / i_count,
        i_mb_count[I_4x4]  / i_count );
    if( b_print_pcm )
        sprintf( intra, " %4.1f%%", i_mb_count[I_PCM]  / i_count );
}

/****************************************************************************
 * x264_encoder_close:
 ****************************************************************************/
void    x264_encoder_close  ( x264_t *h )
{
    int64_t i_yuv_size = FRAME_SIZE( h->param.i_width * h->param.i_height );
    int64_t i_mb_count_size[2][7] = {{0}};
    char buf[200];
    int b_print_pcm = h->stat.i_mb_count[SLICE_TYPE_I][I_PCM]
                   || h->stat.i_mb_count[SLICE_TYPE_P][I_PCM];

    x264_lookahead_delete( h );

    h->i_frame++;

    /* Slices used and PSNR */
    for( int i = 0; i < 2; i++ )
    {
        static const uint8_t slice_order[] = { SLICE_TYPE_I, SLICE_TYPE_P};
        int i_slice = slice_order[i];

        if( h->stat.i_frame_count[i_slice] > 0 )
        {
            int i_count = h->stat.i_frame_count[i_slice];
            double dur =  h->stat.f_frame_duration[i_slice];
            if( h->param.analyse.b_psnr )
            {
                x264_log( h, X264_LOG_INFO,
                          "frame %c:%-5d Avg QP:%5.2f  size:%6.0f  PSNR Mean Y:%5.2f U:%5.2f V:%5.2f Avg:%5.2f Global:%5.2f\n",
                          slice_type_to_char[i_slice],
                          i_count,
                          h->stat.f_frame_qp[i_slice] / i_count,
                          (double)h->stat.i_frame_size[i_slice] / i_count,
                          h->stat.f_psnr_mean_y[i_slice] / dur, h->stat.f_psnr_mean_u[i_slice] / dur, h->stat.f_psnr_mean_v[i_slice] / dur,
                          h->stat.f_psnr_average[i_slice] / dur,
                          x264_psnr( h->stat.f_ssd_global[i_slice], dur * i_yuv_size ) );
            }
            else
            {
                x264_log( h, X264_LOG_INFO,
                          "frame %c:%-5d Avg QP:%5.2f  size:%6.0f\n",
                          slice_type_to_char[i_slice],
                          i_count,
                          h->stat.f_frame_qp[i_slice] / i_count,
                          (double)h->stat.i_frame_size[i_slice] / i_count );
            }
        }
    }

    for( int i_type = 0; i_type < 2; i_type++ )
        for( int i = 0; i < X264_PARTTYPE_MAX; i++ )
        {
            i_mb_count_size[i_type][x264_mb_partition_pixel_table[i]] += h->stat.i_mb_partition[i_type][i];
        }

    /* MB types used */
    if( h->stat.i_frame_count[SLICE_TYPE_I] > 0 )
    {
        int64_t *i_mb_count = h->stat.i_mb_count[SLICE_TYPE_I];
        double i_count = h->stat.i_frame_count[SLICE_TYPE_I] * h->mb.i_mb_count / 100.0;
        x264_print_intra( i_mb_count, i_count, b_print_pcm, buf );
        x264_log( h, X264_LOG_INFO, "mb I  %s\n", buf );
    }
    if( h->stat.i_frame_count[SLICE_TYPE_P] > 0 )
    {
        int64_t *i_mb_count = h->stat.i_mb_count[SLICE_TYPE_P];
        double i_count = h->stat.i_frame_count[SLICE_TYPE_P] * h->mb.i_mb_count / 100.0;
        int64_t *i_mb_size = i_mb_count_size[SLICE_TYPE_P];
        x264_print_intra( i_mb_count, i_count, b_print_pcm, buf );
        x264_log( h, X264_LOG_INFO,
                  "mb P  %s  P16..4: %4.1f%% %4.1f%% %4.1f%% %4.1f%% %4.1f%%    skip:%4.1f%%\n",
                  buf,
                  i_mb_size[PIXEL_16x16] / (i_count*4),
                  (i_mb_size[PIXEL_16x8] + i_mb_size[PIXEL_8x16]) / (i_count*4),
                  i_mb_size[PIXEL_8x8] / (i_count*4),
                  (i_mb_size[PIXEL_8x4] + i_mb_size[PIXEL_4x8]) / (i_count*4),
                  i_mb_size[PIXEL_4x4] / (i_count*4),
                  i_mb_count[P_SKIP] / i_count );
    }

    x264_ratecontrol_summary( h );

    if( h->stat.i_frame_count[SLICE_TYPE_I] + h->stat.i_frame_count[SLICE_TYPE_P]  > 0 )
    {
#define SUM3(p) (p[SLICE_TYPE_I] + p[SLICE_TYPE_P] )
#define SUM3b(p,o) (p[SLICE_TYPE_I][o] + p[SLICE_TYPE_P][o] )
        int64_t i_i8x8 = SUM3b( h->stat.i_mb_count, I_8x8 );
        int64_t i_intra = i_i8x8 + SUM3b( h->stat.i_mb_count, I_4x4 )
                                 + SUM3b( h->stat.i_mb_count, I_16x16 );
        int64_t i_all_intra = i_intra + SUM3b( h->stat.i_mb_count, I_PCM);
        const int i_count = h->stat.i_frame_count[SLICE_TYPE_I] +
                            h->stat.i_frame_count[SLICE_TYPE_P] ;
        int64_t i_mb_count = (int64_t)i_count * h->mb.i_mb_count;
        const double duration = h->stat.f_frame_duration[SLICE_TYPE_I] +
                                h->stat.f_frame_duration[SLICE_TYPE_P];
        float f_bitrate = SUM3(h->stat.i_frame_size) / duration / 125;

        buf[0] = 0;
        int csize = 1;
        if( i_mb_count != i_all_intra )
            sprintf( buf, " inter: %.1f%% %.1f%% %.1f%%",
                     h->stat.i_mb_cbp[1] * 100.0 / ((i_mb_count - i_all_intra)*4),
                     h->stat.i_mb_cbp[3] * 100.0 / ((i_mb_count - i_all_intra)*csize),
                     h->stat.i_mb_cbp[5] * 100.0 / ((i_mb_count - i_all_intra)*csize) );
        x264_log( h, X264_LOG_INFO, "coded y,%s,%s intra: %.1f%% %.1f%% %.1f%%%s\n",
                  "uvDC", "uvAC",
                  h->stat.i_mb_cbp[0] * 100.0 / (i_all_intra*4),
                  h->stat.i_mb_cbp[2] * 100.0 / (i_all_intra*csize),
                  h->stat.i_mb_cbp[4] * 100.0 / (i_all_intra*csize), buf );

        int64_t fixed_pred_modes[4][9] = {{0}};
        int64_t sum_pred_modes[4] = {0};
        for( int i = 0; i <= I_PRED_16x16_DC_128; i++ )
        {
            fixed_pred_modes[0][x264_mb_pred_mode16x16_fix[i]] += h->stat.i_mb_pred_mode[0][i];
            sum_pred_modes[0] += h->stat.i_mb_pred_mode[0][i];
        }
        if( sum_pred_modes[0] )
            x264_log( h, X264_LOG_INFO, "i16 v,h,dc,p: %2.0f%% %2.0f%% %2.0f%% %2.0f%%\n",
                      fixed_pred_modes[0][0] * 100.0 / sum_pred_modes[0],
                      fixed_pred_modes[0][1] * 100.0 / sum_pred_modes[0],
                      fixed_pred_modes[0][2] * 100.0 / sum_pred_modes[0],
                      fixed_pred_modes[0][3] * 100.0 / sum_pred_modes[0] );
        for( int i = 1; i <= 2; i++ )
        {
            if( sum_pred_modes[i] )
                x264_log( h, X264_LOG_INFO, "i%d v,h,dc,ddl,ddr,vr,hd,vl,hu: %2.0f%% %2.0f%% %2.0f%% %2.0f%% %2.0f%% %2.0f%% %2.0f%% %2.0f%% %2.0f%%\n", (3-i)*4,
                          fixed_pred_modes[i][0] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][1] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][2] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][3] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][4] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][5] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][6] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][7] * 100.0 / sum_pred_modes[i],
                          fixed_pred_modes[i][8] * 100.0 / sum_pred_modes[i] );
        }
        for( int i = 0; i <= I_PRED_CHROMA_DC_128; i++ )
        {
            fixed_pred_modes[3][x264_mb_chroma_pred_mode_fix[i]] += h->stat.i_mb_pred_mode[3][i];
            sum_pred_modes[3] += h->stat.i_mb_pred_mode[3][i];
        }
        if( sum_pred_modes[3])
            x264_log( h, X264_LOG_INFO, "i8c dc,h,v,p: %2.0f%% %2.0f%% %2.0f%% %2.0f%%\n",
                      fixed_pred_modes[3][0] * 100.0 / sum_pred_modes[3],
                      fixed_pred_modes[3][1] * 100.0 / sum_pred_modes[3],
                      fixed_pred_modes[3][2] * 100.0 / sum_pred_modes[3],
                      fixed_pred_modes[3][3] * 100.0 / sum_pred_modes[3] );

            for( int i_slice = 0; i_slice < 2; i_slice++ )
            {
                char *p = buf;
                int64_t i_den = 0;
                int i_max = 0;
                for( int i = 0; i < X264_REF_MAX*2; i++ )
                    if( h->stat.i_mb_count_ref[i_slice][i] )
                    {
                        i_den += h->stat.i_mb_count_ref[i_slice][i];
                        i_max = i;
                    }
                if( i_max == 0 )
                    continue;
                for( int i = 0; i <= i_max; i++ )
                    p += sprintf( p, " %4.1f%%", 100. * h->stat.i_mb_count_ref[i_slice][i] / i_den );
                x264_log( h, X264_LOG_INFO, "ref %c L%d:%s\n", "PB"[i_slice], 0, buf );
            }

        if( h->param.analyse.b_ssim )
        {
            float ssim = SUM3( h->stat.f_ssim_mean_y ) / duration;
            x264_log( h, X264_LOG_INFO, "SSIM Mean Y:%.7f (%6.3fdb)\n", ssim, x264_ssim( ssim ) );
        }
        if( h->param.analyse.b_psnr )
        {
            x264_log( h, X264_LOG_INFO,
                      "PSNR Mean Y:%6.3f U:%6.3f V:%6.3f Avg:%6.3f Global:%6.3f kb/s:%.2f\n",
                      SUM3( h->stat.f_psnr_mean_y ) / duration,
                      SUM3( h->stat.f_psnr_mean_u ) / duration,
                      SUM3( h->stat.f_psnr_mean_v ) / duration,
                      SUM3( h->stat.f_psnr_average ) / duration,
                      x264_psnr( SUM3( h->stat.f_ssd_global ), duration * i_yuv_size ),
                      f_bitrate );
        }
        else
            x264_log( h, X264_LOG_INFO, "kb/s:%.2f\n", f_bitrate );
    }

    /* rc */
    x264_ratecontrol_delete( h );

    x264_cqm_delete( h );
    x264_free( h->nal_buffer );
    x264_analyse_free_costs( h );

    /* frames */
    x264_frame_delete_list( h->frames.unused[0] );
    x264_frame_delete_list( h->frames.unused[1] );
    if(h->frames.current){
        x264_frame_delete( h->frames.current );
    }
    x264_frame_delete_list( h->frames.blank_unused );

    for( int j = 0; j < h->i_ref; j++ )
        if( h->fref[j] && h->fref[j]->b_duplicate )
            x264_frame_delete( h->fref[j] );

    {
        x264_frame_t **frame;

        {
            for( frame = h->frames.reference; *frame; frame++ )
            {
                assert( (*frame)->i_reference_count > 0 );
                (*frame)->i_reference_count--;
                if( (*frame)->i_reference_count == 0 )
                    x264_frame_delete( *frame );
            }
            frame = &h->fdec;
            if( *frame )
            {
                assert( (*frame)->i_reference_count > 0 );
                (*frame)->i_reference_count--;
                if( (*frame)->i_reference_count == 0 )
                    x264_frame_delete( *frame );
            }
            x264_macroblock_cache_free( h );
        }
        x264_macroblock_thread_free( h );
        x264_free( h->out.p_bitstream );
        x264_free( h->out.nal );
    }
}

int x264_encoder_delayed_frames( x264_t *h )
{
    int delayed_frames = 0;

    delayed_frames += h->frames.current?1:0;

    delayed_frames += h->lookahead->current_frame?1:0;

    assert(delayed_frames == 0);
    return delayed_frames;
}

int x264_encoder_maximum_delayed_frames( x264_t *h )
{
    return h->frames.i_delay;
}
