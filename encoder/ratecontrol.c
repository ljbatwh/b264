/*****************************************************************************
 * ratecontrol.c: ratecontrol
 *****************************************************************************
 * Copyright (C) 2005-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
 *          Michael Niedermayer <michaelni@gmx.at>
 *          Gabriel Bouvigne <gabriel.bouvigne@joost.com>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *          M錸s Rullg錼d <mru@mru.ath.cx>
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

#define _ISOC99_SOURCE
#undef NDEBUG // always check asserts, the speed effect is far too small to disable them

#include "common/common.h"
#include "ratecontrol.h"
#include "me.h"

typedef struct
{
    int pict_type;
    int frame_type;
    int kept_as_ref;
    double qscale;
    int mv_bits;
    int tex_bits;
    int misc_bits;
    uint64_t expected_bits; /*total expected bits up to the current frame (current one excluded)*/
    double expected_vbv;
    double new_qscale;
    int new_qp;
    int i_count;
    int p_count;
    int s_count;
    float blurred_complexity;
    char direct_mode;
    int refcount[16];
    int refs;
    int64_t i_duration;
    int64_t i_cpb_duration;
} ratecontrol_entry_t;

typedef struct
{
    float coeff_min;
    float coeff;
    float count;
    float decay;
    float offset;
} predictor_t;

struct x264_ratecontrol_t
{
    /* constants */
    int b_abr;
    int b_vbv;
    int b_vbv_min_rate;
    double fps;
    double bitrate;
    double rate_tolerance;
    double qcompress;
    int nmb;                    /* number of macroblocks in a frame */
    int qp_constant[3];

    /* current frame */
    ratecontrol_entry_t *rce;
    int qp;                     /* qp for current frame */
    float qpm;                  /* qp for current macroblock: precise float for AQ */
    float qpa_rc;               /* average of macroblocks' qp before aq */
    float qpa_rc_prev;
    int   qpa_aq;               /* average of macroblocks' qp after aq */
    int   qpa_aq_prev;
    float qp_novbv;             /* QP for the current frame if 1-pass VBV was disabled. */

    /* VBV stuff */
    double buffer_size;
    int64_t buffer_fill_final;
    double buffer_fill;         /* planned buffer, if all in-progress frames hit their bit budget */
    double buffer_rate;         /* # of bits added to buffer_fill after each frame */
    double vbv_max_rate;        /* # of bits added to buffer_fill per second */
    predictor_t *pred;          /* predict frame size from satd */
    int single_frame_vbv;
    double rate_factor_max_increment; /* Don't allow RF above (CRF + this value). */

    /* ABR stuff */
    int    last_satd;
    double last_rceq;
    double cplxr_sum;           /* sum of bits*qscale/rceq */
    int64_t filler_bits_sum;    /* sum in bits of finished frames' filler data */
    double wanted_bits_window;  /* target bitrate * window */
    double cbr_decay;
    double short_term_cplxsum;
    double short_term_cplxcount;
    double rate_factor_constant;
    double ip_offset;
    double pb_offset;

    int num_entries;            /* number of ratecontrol_entry_ts */
    ratecontrol_entry_t *entry; /* FIXME: copy needed data and free this once init is done */
    double last_qscale;
    double last_qscale_for[3];  /* last qscale for a specific pict type, used for max_diff & ipb factor stuff */
    int last_non_b_pict_type;
    double accum_p_qp;          /* for determining I-frame quant */
    double accum_p_norm;
    double last_accum_p_norm;
    double lmin[3];             /* min qscale by frame type */
    double lmax[3];
    double lstep;               /* max change (multiply) in qscale per frame */

    /* MBRC stuff */
    float frame_size_estimated; /* Access to this variable must be atomic: double is
                                 * not atomic on all arches we care about */
    double frame_size_maximum;  /* Maximum frame size due to MinCR */
    double frame_size_planned;
    double slice_size_planned;
    predictor_t (*row_pred)[2];
    predictor_t row_preds[3][2];

    /* hrd stuff */
    int initial_cpb_removal_delay;
    int initial_cpb_removal_delay_offset;
    double nrt_first_access_unit; /* nominal removal time */
    double previous_cpb_final_arrival_time;
    uint64_t hrd_multiply_denom;
};

static float rate_estimate_qscale( x264_t *h );
static int update_vbv( x264_t *h, int bits );
static void update_vbv_plan( x264_t *h, int overhead );
static float predict_size( predictor_t *p, float q, float var );
static void update_predictor( predictor_t *p, float q, float var, float bits );

#define CMP_OPT_FIRST_PASS( opt, param_val )\
{\
    if( ( p = strstr( opts, opt "=" ) ) && sscanf( p, opt "=%d" , &i ) && param_val != i )\
    {\
        x264_log( h, X264_LOG_ERROR, "different " opt " setting than first pass (%d vs %d)\n", param_val, i );\
        return -1;\
    }\
}

/* Terminology:
 * qp = h.264's quantizer
 * qscale = linearized quantizer = Lagrange multiplier
 */
static  float qp2qscale( float qp )
{
    return 0.85f * powf( 2.0f, ( qp - 12.0f ) / 6.0f );
}
static  float qscale2qp( float qscale )
{
    return 12.0f + 6.0f * log2f( qscale/0.85f );
}

static  uint32_t ac_energy_var( uint64_t sum_ssd, int shift, x264_frame_t *frame, int i, int b_store )
{
    uint32_t sum = sum_ssd;
    uint32_t ssd = sum_ssd >> 32;
    if( b_store )
    {
        frame->i_pixel_sum[i] += sum;
        frame->i_pixel_ssd[i] += ssd;
    }
    return ssd - ((uint64_t)sum * sum >> shift);
}

static  uint32_t ac_energy_plane( x264_t *h, int mb_x, int mb_y, x264_frame_t *frame, int i, int b_chroma, int b_field, int b_store )
{
    int height = b_chroma ? 16>>CHROMA_V_SHIFT : 16;
    int stride = frame->i_stride[i];
    int offset = b_field
        ? 16 * mb_x + height * (mb_y&~1) * stride + (mb_y&1) * stride
        : 16 * mb_x + height * mb_y * stride;
    stride <<= b_field;
    if( b_chroma )
    {
        ALIGNED_ARRAY_16( pixel, pix,[FENC_STRIDE*16] );
        int chromapix = h->luma2chroma_pixel[PIXEL_16x16];
        int shift = 7 - CHROMA_V_SHIFT;

        h->mc.load_deinterleave_chroma_fenc( pix, frame->plane[1] + offset, stride, height );
        return ac_energy_var( h->pixf.var[chromapix]( pix,               FENC_STRIDE ), shift, frame, 1, b_store )
             + ac_energy_var( h->pixf.var[chromapix]( pix+FENC_STRIDE/2, FENC_STRIDE ), shift, frame, 2, b_store );
    }
    else
        return ac_energy_var( h->pixf.var[PIXEL_16x16]( frame->plane[i] + offset, stride ), 8, frame, i, b_store );
}

// Find the total AC energy of the block in all planes.
static NOINLINE uint32_t x264_ac_energy_mb( x264_t *h, int mb_x, int mb_y, x264_frame_t *frame )
{
    /* This function contains annoying hacks because GCC has a habit of reordering emms
     * and putting it after floating point ops.  As a result, we put the emms at the end of the
     * function and make sure that its always called before the float math.  Noinline makes
     * sure no reordering goes on. */
    uint32_t var;
    {
        var  = ac_energy_plane( h, mb_x, mb_y, frame, 0, 0, 0, 1 );
        var += ac_energy_plane( h, mb_x, mb_y, frame, 1, 1, 0, 1 );
    }

    return var;
}

void x264_adaptive_quant_frame( x264_t *h, x264_frame_t *frame )
{
    /* constants chosen to result in approximately the same overall bitrate as without AQ.
     * FIXME: while they're written in 5 significant digits, they're only tuned to 2. */
    float strength;
    float avg_adj = 0.f;
    /* Initialize frame stats */
    for( int i = 0; i < 3; i++ )
    {
        frame->i_pixel_sum[i] = 0;
        frame->i_pixel_ssd[i] = 0;
    }

    /* Degenerate cases */
    if( h->param.rc.i_aq_mode == X264_AQ_NONE || h->param.rc.f_aq_strength == 0 )
    {
        /* Need to init it anyways for MB tree */
        if( h->param.rc.i_aq_mode && h->param.rc.f_aq_strength == 0 )
        {
            {
                memset( frame->f_qp_offset, 0, h->mb.i_mb_count * sizeof(float) );
                memset( frame->f_qp_offset_aq, 0, h->mb.i_mb_count * sizeof(float) );
                if( h->frames.b_have_lowres )
                    for( int mb_xy = 0; mb_xy < h->mb.i_mb_count; mb_xy++ )
                        frame->i_inv_qscale_factor[mb_xy] = 256;
            }
        }
        return;
    }
    /* Actual adaptive quantization */
    else
    {
        if( h->param.rc.i_aq_mode == X264_AQ_AUTOVARIANCE )
        {
            float bit_depth_correction = powf(1 << (BIT_DEPTH-8), 0.5f);
            float avg_adj_pow2 = 0.f;
            for( int mb_y = 0; mb_y < h->mb.i_mb_height; mb_y++ )
                for( int mb_x = 0; mb_x < h->mb.i_mb_width; mb_x++ )
                {
                    uint32_t energy = x264_ac_energy_mb( h, mb_x, mb_y, frame );
                    float qp_adj = powf( energy + 1, 0.125f );
                    frame->f_qp_offset[mb_x + mb_y*h->mb.i_mb_stride] = qp_adj;
                    avg_adj += qp_adj;
                    avg_adj_pow2 += qp_adj * qp_adj;
                }
            avg_adj /= h->mb.i_mb_count;
            avg_adj_pow2 /= h->mb.i_mb_count;
            strength = h->param.rc.f_aq_strength * avg_adj / bit_depth_correction;
            avg_adj = avg_adj - 0.5f * (avg_adj_pow2 - (14.f * bit_depth_correction)) / avg_adj;
        }
        else
            strength = h->param.rc.f_aq_strength * 1.0397f;

        for( int mb_y = 0; mb_y < h->mb.i_mb_height; mb_y++ )
            for( int mb_x = 0; mb_x < h->mb.i_mb_width; mb_x++ )
            {
                float qp_adj;
                int mb_xy = mb_x + mb_y*h->mb.i_mb_stride;
                if( h->param.rc.i_aq_mode == X264_AQ_AUTOVARIANCE )
                {
                    qp_adj = frame->f_qp_offset[mb_xy];
                    qp_adj = strength * (qp_adj - avg_adj);
                }
                else
                {
                    uint32_t energy = x264_ac_energy_mb( h, mb_x, mb_y, frame );
                    qp_adj = strength * (x264_log2( X264_MAX(energy, 1) ) - (14.427f + 2*(BIT_DEPTH-8)));
                }
                frame->f_qp_offset[mb_xy] =
                frame->f_qp_offset_aq[mb_xy] = qp_adj;
                if( h->frames.b_have_lowres )
                    frame->i_inv_qscale_factor[mb_xy] = x264_exp2fix8(qp_adj);
            }
    }

    /* Remove mean from SSD calculation */
    for( int i = 0; i < 3; i++ )
    {
        uint64_t ssd = frame->i_pixel_ssd[i];
        uint64_t sum = frame->i_pixel_sum[i];
        int width  = 16*h->mb.i_mb_width  >> (i && CHROMA_H_SHIFT);
        int height = 16*h->mb.i_mb_height >> (i && CHROMA_V_SHIFT);
        frame->i_pixel_ssd[i] = ssd - (sum * sum + width * height / 2) / (width * height);
    }
}


int x264_reference_build_list_optimal( x264_t *h )
{
    ratecontrol_entry_t *rce = h->rc->rce;
    x264_frame_t *frames[16];
    int refcount[16];

    if( rce->refs != h->i_ref )
        return -1;

    memcpy( frames, h->fref, sizeof(frames) );
    memcpy( refcount, rce->refcount, sizeof(refcount) );

    /* For now don't reorder ref 0; it seems to lower quality
       in most cases due to skips. */
    for( int ref = 1; ref < h->i_ref; ref++ )
    {
        int max = -1;
        int bestref = 1;

        for( int i = 1; i < h->i_ref; i++ )
            /* Favor lower POC as a tiebreaker. */
            COPY2_IF_GT( max, refcount[i], bestref, i );

        /* FIXME: If there are duplicates from frames other than ref0 then it is possible
         * that the optimal ordering doesnt place every duplicate. */

        refcount[bestref] = -1;
        h->fref[ref] = frames[bestref];
    }

    return 0;
}

void x264_ratecontrol_init_reconfigurable( x264_t *h, int b_init )
{
    x264_ratecontrol_t *rc = h->rc;

    if( h->param.rc.i_rc_method == X264_RC_CRF )
    {
        /* Arbitrary rescaling to make CRF somewhat similar to QP.
         * Try to compensate for MB-tree's effects as well. */
        double base_cplx = h->mb.i_mb_count * (80);
        rc->rate_factor_constant = pow( base_cplx, 1 - rc->qcompress )
                                 / qp2qscale( h->param.rc.f_rf_constant  );
    }

    if( h->param.rc.i_vbv_max_bitrate > 0 && h->param.rc.i_vbv_buffer_size > 0 )
    {
        /* We don't support changing the ABR bitrate right now,
           so if the stream starts as CBR, keep it CBR. */
        if( rc->b_vbv_min_rate )
            h->param.rc.i_vbv_max_bitrate = h->param.rc.i_bitrate;

        if( h->param.rc.i_vbv_buffer_size < (int)(h->param.rc.i_vbv_max_bitrate / rc->fps) )
        {
            h->param.rc.i_vbv_buffer_size = h->param.rc.i_vbv_max_bitrate / rc->fps;
            x264_log( h, X264_LOG_WARNING, "VBV buffer size cannot be smaller than one frame, using %d kbit\n",
                      h->param.rc.i_vbv_buffer_size );
        }

        int kilobit_size = 1000;
        int vbv_buffer_size = h->param.rc.i_vbv_buffer_size * kilobit_size;
        int vbv_max_bitrate = h->param.rc.i_vbv_max_bitrate * kilobit_size;

        /* Init HRD */
        if( h->param.i_nal_hrd && b_init )
        {
            h->sps->vui.hrd.i_cpb_cnt = 1;
            h->sps->vui.hrd.b_cbr_hrd = h->param.i_nal_hrd == X264_NAL_HRD_CBR;
            h->sps->vui.hrd.i_time_offset_length = 0;

            #define BR_SHIFT  6
            #define CPB_SHIFT 4

            // normalize HRD size and rate to the value / scale notation
            h->sps->vui.hrd.i_bit_rate_scale = x264_clip3( x264_ctz( vbv_max_bitrate ) - BR_SHIFT, 0, 15 );
            h->sps->vui.hrd.i_bit_rate_value = vbv_max_bitrate >> ( h->sps->vui.hrd.i_bit_rate_scale + BR_SHIFT );
            h->sps->vui.hrd.i_bit_rate_unscaled = h->sps->vui.hrd.i_bit_rate_value << ( h->sps->vui.hrd.i_bit_rate_scale + BR_SHIFT );
            h->sps->vui.hrd.i_cpb_size_scale = x264_clip3( x264_ctz( vbv_buffer_size ) - CPB_SHIFT, 0, 15 );
            h->sps->vui.hrd.i_cpb_size_value = vbv_buffer_size >> ( h->sps->vui.hrd.i_cpb_size_scale + CPB_SHIFT );
            h->sps->vui.hrd.i_cpb_size_unscaled = h->sps->vui.hrd.i_cpb_size_value << ( h->sps->vui.hrd.i_cpb_size_scale + CPB_SHIFT );

            #undef CPB_SHIFT
            #undef BR_SHIFT

            // arbitrary
            #define MAX_DURATION 0.5

            int max_cpb_output_delay = X264_MIN( h->param.i_keyint_max * MAX_DURATION * h->sps->vui.i_time_scale / h->sps->vui.i_num_units_in_tick, INT_MAX );
            int max_dpb_output_delay = h->sps->vui.i_max_dec_frame_buffering * MAX_DURATION * h->sps->vui.i_time_scale / h->sps->vui.i_num_units_in_tick;
            int max_delay = (int)(90000.0 * (double)h->sps->vui.hrd.i_cpb_size_unscaled / h->sps->vui.hrd.i_bit_rate_unscaled + 0.5);

            h->sps->vui.hrd.i_initial_cpb_removal_delay_length = 2 + x264_clip3( 32 - x264_clz( max_delay ), 4, 22 );
            h->sps->vui.hrd.i_cpb_removal_delay_length = x264_clip3( 32 - x264_clz( max_cpb_output_delay ), 4, 31 );
            h->sps->vui.hrd.i_dpb_output_delay_length  = x264_clip3( 32 - x264_clz( max_dpb_output_delay ), 4, 31 );

            #undef MAX_DURATION

            vbv_buffer_size = h->sps->vui.hrd.i_cpb_size_unscaled;
            vbv_max_bitrate = h->sps->vui.hrd.i_bit_rate_unscaled;
        }
        else if( h->param.i_nal_hrd && !b_init )
        {
            x264_log( h, X264_LOG_WARNING, "VBV parameters cannot be changed when NAL HRD is in use\n" );
            return;
        }
        h->sps->vui.hrd.i_bit_rate_unscaled = vbv_max_bitrate;
        h->sps->vui.hrd.i_cpb_size_unscaled = vbv_buffer_size;

        if( rc->b_vbv_min_rate )
            rc->bitrate = (double)h->param.rc.i_bitrate * kilobit_size;
        rc->buffer_rate = vbv_max_bitrate / rc->fps;
        rc->vbv_max_rate = vbv_max_bitrate;
        rc->buffer_size = vbv_buffer_size;
        rc->single_frame_vbv = rc->buffer_rate * 1.1 > rc->buffer_size;
        rc->cbr_decay = 1.0 - rc->buffer_rate / rc->buffer_size
                      * 0.5 * X264_MAX(0, 1.5 - rc->buffer_rate * rc->fps / rc->bitrate);
        if( h->param.rc.i_rc_method == X264_RC_CRF && h->param.rc.f_rf_constant_max )
        {
            rc->rate_factor_max_increment = h->param.rc.f_rf_constant_max - h->param.rc.f_rf_constant;
            if( rc->rate_factor_max_increment <= 0 )
            {
                x264_log( h, X264_LOG_WARNING, "CRF max must be greater than CRF\n" );
                rc->rate_factor_max_increment = 0;
            }
        }
        if( b_init )
        {
            if( h->param.rc.f_vbv_buffer_init > 1. )
                h->param.rc.f_vbv_buffer_init = x264_clip3f( h->param.rc.f_vbv_buffer_init / h->param.rc.i_vbv_buffer_size, 0, 1 );
            h->param.rc.f_vbv_buffer_init = x264_clip3f( X264_MAX( h->param.rc.f_vbv_buffer_init, rc->buffer_rate / rc->buffer_size ), 0, 1);
            rc->buffer_fill_final = rc->buffer_size * h->param.rc.f_vbv_buffer_init * h->sps->vui.i_time_scale;
            rc->b_vbv = 1;
            rc->b_vbv_min_rate = h->param.rc.i_rc_method == X264_RC_ABR
                          && h->param.rc.i_vbv_max_bitrate <= h->param.rc.i_bitrate;
        }
    }
}

int x264_ratecontrol_new( x264_t *h )
{
    x264_ratecontrol_t *rc;



    CHECKED_MALLOCZERO( h->rc, sizeof(x264_ratecontrol_t) );
    rc = h->rc;

    rc->b_abr = h->param.rc.i_rc_method != X264_RC_CQP;

    /* FIXME: use integers */
    if( h->param.i_fps_num > 0 && h->param.i_fps_den > 0 )
        rc->fps = (float) h->param.i_fps_num / h->param.i_fps_den;
    else
        rc->fps = 25.0;

    rc->qcompress = h->param.rc.f_qcompress;

    rc->bitrate = h->param.rc.i_bitrate * (1000.);
    rc->rate_tolerance = h->param.rc.f_rate_tolerance;
    rc->nmb = h->mb.i_mb_count;
    rc->last_non_b_pict_type = -1;
    rc->cbr_decay = 1.0;

    x264_ratecontrol_init_reconfigurable( h, 1 );

    if( h->param.i_nal_hrd )
    {
        uint64_t denom = (uint64_t)h->sps->vui.hrd.i_bit_rate_unscaled * h->sps->vui.i_time_scale;
        uint64_t num = 180000;
        x264_reduce_fraction64( &num, &denom );
        rc->hrd_multiply_denom = 180000 / num;

        double bits_required = log2( 180000 / rc->hrd_multiply_denom )
                             + log2( h->sps->vui.i_time_scale )
                             + log2( h->sps->vui.hrd.i_cpb_size_unscaled );
        if( bits_required >= 63 )
        {
            x264_log( h, X264_LOG_ERROR, "HRD with very large timescale and bufsize not supported\n" );
            return -1;
        }
    }

    if( rc->rate_tolerance < 0.01 )
    {
        x264_log( h, X264_LOG_WARNING, "bitrate tolerance too small, using .01\n" );
        rc->rate_tolerance = 0.01;
    }

    h->mb.b_variable_qp = rc->b_vbv || h->param.rc.i_aq_mode;

    if( rc->b_abr )
    {
        /* FIXME ABR_INIT_QP is actually used only in CRF */
#define ABR_INIT_QP (( h->param.rc.i_rc_method == X264_RC_CRF ? h->param.rc.f_rf_constant : 24 ) )
        rc->accum_p_norm = .01;
        rc->accum_p_qp = ABR_INIT_QP * rc->accum_p_norm;
        /* estimated ratio that produces a reasonable QP for the first I-frame */
        rc->cplxr_sum = .01 * pow( 7.0e5, rc->qcompress ) * pow( h->mb.i_mb_count, 0.5 );
        rc->wanted_bits_window = 1.0 * rc->bitrate / rc->fps;
        rc->last_non_b_pict_type = SLICE_TYPE_I;
    }

    rc->ip_offset = 6.0 * log2f( h->param.rc.f_ip_factor );
    rc->pb_offset = 6.0 * log2f( h->param.rc.f_pb_factor );
    rc->qp_constant[SLICE_TYPE_P] = h->param.rc.i_qp_constant;
    rc->qp_constant[SLICE_TYPE_I] = x264_clip3( h->param.rc.i_qp_constant - rc->ip_offset + 0.5, 0, QP_MAX );
    h->mb.ip_offset = rc->ip_offset + 0.5;

    rc->lstep = pow( 2, h->param.rc.i_qp_step / 6.0 );
    rc->last_qscale = qp2qscale( 26 );
    int num_preds = 1;
    CHECKED_MALLOC( rc->pred, 5 * sizeof(predictor_t) * num_preds );
    for( int i = 0; i < 3; i++ )
    {
        rc->last_qscale_for[i] = qp2qscale( ABR_INIT_QP );
        rc->lmin[i] = qp2qscale( h->param.rc.i_qp_min );
        rc->lmax[i] = qp2qscale( h->param.rc.i_qp_max );
        for( int j = 0; j < num_preds; j++ )
        {
            rc->pred[i+j*5].coeff_min = 2.0 / 4;
            rc->pred[i+j*5].coeff = 2.0;
            rc->pred[i+j*5].count = 1.0;
            rc->pred[i+j*5].decay = 0.5;
            rc->pred[i+j*5].offset = 0.0;
        }
        for( int j = 0; j < 2; j++ )
        {
            rc->row_preds[i][j].coeff_min = .25 / 4;
            rc->row_preds[i][j].coeff = .25;
            rc->row_preds[i][j].count = 1.0;
            rc->row_preds[i][j].decay = 0.5;
            rc->row_preds[i][j].offset = 0.0;
        }
    }

    return 0;
fail:
    return -1;
}

void x264_ratecontrol_summary( x264_t *h )
{
    x264_ratecontrol_t *rc = h->rc;
    if( rc->b_abr && h->param.rc.i_rc_method == X264_RC_ABR && rc->cbr_decay > .9999 )
    {
        double base_cplx = h->mb.i_mb_count * (80);
        x264_log( h, X264_LOG_INFO, "final ratefactor: %.2f\n",
                  qscale2qp( pow( base_cplx, 1 - rc->qcompress )
                             * rc->cplxr_sum / rc->wanted_bits_window ) );
    }
}

void x264_ratecontrol_delete( x264_t *h )
{
    x264_ratecontrol_t *rc = h->rc;

    x264_free( rc->pred );
    x264_free( rc->entry );
    x264_free( rc );
}

static void accum_p_qp_update( x264_t *h, float qp )
{
    x264_ratecontrol_t *rc = h->rc;
    rc->accum_p_qp   *= .95;
    rc->accum_p_norm *= .95;
    rc->accum_p_norm += 1;
    if( h->sh.i_type == SLICE_TYPE_I )
        rc->accum_p_qp += qp + rc->ip_offset;
    else
        rc->accum_p_qp += qp;
}

/* Before encoding a frame, choose a QP for it */
void x264_ratecontrol_start( x264_t *h, int i_force_qp, int overhead )
{
    x264_ratecontrol_t *rc = h->rc;
    ratecontrol_entry_t *rce = NULL;
    float q;

    if( rc->b_vbv )
    {
        memset( h->fdec->i_row_bits, 0, h->mb.i_mb_height * sizeof(int) );
        memset( h->fdec->f_row_qp, 0, h->mb.i_mb_height * sizeof(float) );
        memset( h->fdec->f_row_qscale, 0, h->mb.i_mb_height * sizeof(float) );
        rc->row_pred = &rc->row_preds[h->sh.i_type];
        rc->buffer_rate = h->fenc->i_cpb_duration * rc->vbv_max_rate * h->sps->vui.i_num_units_in_tick / h->sps->vui.i_time_scale;
        update_vbv_plan( h, overhead );

        const x264_level_t *l = x264_levels;
        while( l->level_idc != 0 && l->level_idc != h->param.i_level_idc )
            l++;

        int mincr = l->mincr;

        /* Profiles above High don't require minCR, so just set the maximum to a large value. */
        {
            /* The spec has a bizarre special case for the first frame. */
            if( h->i_frame == 0 )
            {
                //384 * ( Max( PicSizeInMbs, fR * MaxMBPS ) + MaxMBPS * ( tr( 0 ) - tr,n( 0 ) ) ) / MinCR
                double fr = 1. / 172;
                int pic_size_in_mbs = h->mb.i_mb_width * h->mb.i_mb_height;
                rc->frame_size_maximum = 384 * BIT_DEPTH * X264_MAX( pic_size_in_mbs, fr*l->mbps ) / mincr;
            }
            else
            {
                //384 * MaxMBPS * ( tr( n ) - tr( n - 1 ) ) / MinCR
                rc->frame_size_maximum = 384 * BIT_DEPTH * ((double)h->fenc->i_cpb_duration * h->sps->vui.i_num_units_in_tick / h->sps->vui.i_time_scale) * l->mbps / mincr;
            }
        }
    }


    if( rc->b_abr )
    {
        q = qscale2qp( rate_estimate_qscale( h ) );
    }
    else /* CQP */
    {
        q = rc->qp_constant[ h->sh.i_type ];

    }
    if( i_force_qp != X264_QP_AUTO )
        q = i_force_qp - 1;

    q = x264_clip3f( q, h->param.rc.i_qp_min, h->param.rc.i_qp_max );

    rc->qpa_rc = rc->qpa_rc_prev =
    rc->qpa_aq = rc->qpa_aq_prev = 0;
    rc->qp = x264_clip3( q + 0.5f, 0, QP_MAX );
    h->fdec->f_qp_avg_rc =
    h->fdec->f_qp_avg_aq =
    rc->qpm = q;
    if( rce )
        rce->new_qp = rc->qp;

    accum_p_qp_update( h, rc->qpm );

    rc->last_non_b_pict_type = h->sh.i_type;
}

static float predict_row_size( x264_t *h, int y, float qscale )
{
    /* average between two predictors:
     * absolute SATD, and scaled bit cost of the colocated row in the previous frame */
    x264_ratecontrol_t *rc = h->rc;
    float pred_s = predict_size( rc->row_pred[0], qscale, h->fdec->i_row_satd[y] );
    if( h->sh.i_type == SLICE_TYPE_I || qscale >= h->fref[0]->f_row_qscale[y] )
    {
        if( h->sh.i_type == SLICE_TYPE_P
            && h->fref[0]->i_type == h->fdec->i_type
            && h->fref[0]->f_row_qscale[y] > 0
            && h->fref[0]->i_row_satd[y] > 0
            && (abs(h->fref[0]->i_row_satd[y] - h->fdec->i_row_satd[y]) < h->fdec->i_row_satd[y]/2))
        {
            float pred_t = h->fref[0]->i_row_bits[y] * h->fdec->i_row_satd[y] / h->fref[0]->i_row_satd[y]
                         * h->fref[0]->f_row_qscale[y] / qscale;
            return (pred_s + pred_t) * 0.5f;
        }
        return pred_s;
    }
    /* Our QP is lower than the reference! */
    else
    {
        float pred_intra = predict_size( rc->row_pred[1], qscale, h->fdec->i_row_satds[0][y] );
        /* Sum: better to overestimate than underestimate by using only one of the two predictors. */
        return pred_intra + pred_s;
    }
}

static int row_bits_so_far( x264_t *h, int y )
{
    int bits = 0;
    for( int i = h->i_threadslice_start; i <= y; i++ )
        bits += h->fdec->i_row_bits[i];
    return bits;
}

static float predict_row_size_sum( x264_t *h, int y, float qp )
{
    float qscale = qp2qscale( qp );
    float bits = row_bits_so_far( h, y );
    for( int i = y+1; i < h->i_threadslice_end; i++ )
        bits += predict_row_size( h, i, qscale );
    return bits;
}

/* TODO:
 *  eliminate all use of qp in row ratecontrol: make it entirely qscale-based.
 *  make this function stop being needlessly O(N^2)
 *  update more often than once per row? */
int x264_ratecontrol_mb( x264_t *h, int bits )
{
    x264_ratecontrol_t *rc = h->rc;
    const int y = h->mb.i_mb_y;

    h->fdec->i_row_bits[y] += bits;
    rc->qpa_aq += h->mb.i_qp;

    if( h->mb.i_mb_x != h->mb.i_mb_width - 1 )
        return 0;


    rc->qpa_rc += rc->qpm * h->mb.i_mb_width;

    if( !rc->b_vbv )
        return 0;

    float qscale = qp2qscale( rc->qpm );
    h->fdec->f_row_qp[y] = rc->qpm;
    h->fdec->f_row_qscale[y] = qscale;

    update_predictor( rc->row_pred[0], qscale, h->fdec->i_row_satd[y], h->fdec->i_row_bits[y] );
    if( h->sh.i_type == SLICE_TYPE_P && rc->qpm < h->fref[0]->f_row_qp[y] )
        update_predictor( rc->row_pred[1], qscale, h->fdec->i_row_satds[0][y], h->fdec->i_row_bits[y] );

    /* FIXME: We don't currently support the case where there's a slice
     * boundary in between. */
    int can_reencode_row = h->sh.i_first_mb <= ((h->mb.i_mb_y) * h->mb.i_mb_stride);

    /* tweak quality based on difference from predicted size */
    float prev_row_qp = h->fdec->f_row_qp[y];
    float qp_absolute_max = h->param.rc.i_qp_max;
    if( rc->rate_factor_max_increment )
        qp_absolute_max = X264_MIN( qp_absolute_max, rc->qp_novbv + rc->rate_factor_max_increment );
    float qp_max = X264_MIN( prev_row_qp + h->param.rc.i_qp_step, qp_absolute_max );
    float qp_min = X264_MAX( prev_row_qp - h->param.rc.i_qp_step, h->param.rc.i_qp_min );
    float step_size = 0.5f;
    float buffer_left_planned = rc->buffer_fill - rc->frame_size_planned;
    float slice_size_planned = rc->frame_size_planned;
    float max_frame_error = X264_MAX( 0.05f, 1.0f / h->mb.i_mb_height );
    float size_of_other_slices = 0;
    if( y < h->i_threadslice_end-1 )
    {
        /* More threads means we have to be more cautious in letting ratecontrol use up extra bits. */
        float rc_tol = buffer_left_planned / rc->rate_tolerance;
        float b1 = predict_row_size_sum( h, y, rc->qpm ) + size_of_other_slices;

        /* Don't increase the row QPs until a sufficent amount of the bits of the frame have been processed, in case a flat */
        /* area at the top of the frame was measured inaccurately. */
        if( row_bits_so_far( h, y ) < 0.05f * slice_size_planned )
            qp_max = qp_absolute_max = prev_row_qp;

        if( h->sh.i_type != SLICE_TYPE_I )
            rc_tol *= 0.5f;

        if( !rc->b_vbv_min_rate )
            qp_min = X264_MAX( qp_min, rc->qp_novbv );

        while( rc->qpm < qp_max
               && ((b1 > rc->frame_size_planned + rc_tol) ||
                   (rc->buffer_fill - b1 < buffer_left_planned * 0.5f) ||
                   (b1 > rc->frame_size_planned && rc->qpm < rc->qp_novbv)) )
        {
            rc->qpm += step_size;
            b1 = predict_row_size_sum( h, y, rc->qpm ) + size_of_other_slices;
        }

        while( rc->qpm > qp_min
               && (rc->qpm > h->fdec->f_row_qp[0] || rc->single_frame_vbv)
               && ((b1 < rc->frame_size_planned * 0.8f && rc->qpm <= prev_row_qp)
               || b1 < (rc->buffer_fill - rc->buffer_size + rc->buffer_rate) * 1.1f) )
        {
            rc->qpm -= step_size;
            b1 = predict_row_size_sum( h, y, rc->qpm ) + size_of_other_slices;
        }

        /* avoid VBV underflow or MinCR violation */
        while( (rc->qpm < qp_absolute_max)
               && ((rc->buffer_fill - b1 < rc->buffer_rate * max_frame_error) ||
                   (rc->frame_size_maximum - b1 < rc->frame_size_maximum * max_frame_error)))
        {
            rc->qpm += step_size;
            b1 = predict_row_size_sum( h, y, rc->qpm ) + size_of_other_slices;
        }

        h->rc->frame_size_estimated = b1 - size_of_other_slices;

        /* If the current row was large enough to cause a large QP jump, try re-encoding it. */
        if( rc->qpm > qp_max && prev_row_qp < qp_max && can_reencode_row )
        {
            /* Bump QP to halfway in between... close enough. */
            rc->qpm = x264_clip3f( (prev_row_qp + rc->qpm)*0.5f, prev_row_qp + 1.0f, qp_max );
            rc->qpa_rc = rc->qpa_rc_prev;
            rc->qpa_aq = rc->qpa_aq_prev;
            h->fdec->i_row_bits[y] = 0;
            return -1;
        }
    }
    else
    {
        h->rc->frame_size_estimated = predict_row_size_sum( h, y, rc->qpm );

        /* Last-ditch attempt: if the last row of the frame underflowed the VBV,
         * try again. */
        if( (h->rc->frame_size_estimated + size_of_other_slices) > (rc->buffer_fill - rc->buffer_rate * max_frame_error) &&
             rc->qpm < qp_max && can_reencode_row )
        {
            rc->qpm = qp_max;
            rc->qpa_rc = rc->qpa_rc_prev;
            rc->qpa_aq = rc->qpa_aq_prev;
            h->fdec->i_row_bits[y] = 0;
            return -1;
        }
    }

    rc->qpa_rc_prev = rc->qpa_rc;
    rc->qpa_aq_prev = rc->qpa_aq;

    return 0;
}

int x264_ratecontrol_qp( x264_t *h )
{
    return x264_clip3( h->rc->qpm + 0.5f, h->param.rc.i_qp_min, h->param.rc.i_qp_max );
}

int x264_ratecontrol_mb_qp( x264_t *h )
{

    float qp = h->rc->qpm;
    if( h->param.rc.i_aq_mode )
    {
         /* MB-tree currently doesn't adjust quantizers in unreferenced frames. */
        float qp_offset = h->fdec->b_kept_as_ref ? h->fenc->f_qp_offset[h->mb.i_mb_xy] : h->fenc->f_qp_offset_aq[h->mb.i_mb_xy];
        /* Scale AQ's effect towards zero in emergency mode. */
        if( qp > QP_MAX_SPEC )
            qp_offset *= (QP_MAX - qp) / (QP_MAX - QP_MAX_SPEC);
        qp += qp_offset;
    }
    return x264_clip3( qp + 0.5f, h->param.rc.i_qp_min, h->param.rc.i_qp_max );
}

/* In 2pass, force the same frame types as in the 1st pass */
int x264_ratecontrol_slice_type( x264_t *h, int frame_num )
{
    return X264_TYPE_AUTO;
}

/* After encoding one frame, save stats and update ratecontrol state */
int x264_ratecontrol_end( x264_t *h, int bits, int *filler )
{
    x264_ratecontrol_t *rc = h->rc;
    const int *mbs = h->stat.frame.i_mb_count;

    h->stat.frame.i_mb_count_skip = mbs[P_SKIP] ;
    h->stat.frame.i_mb_count_i = mbs[I_16x16] + mbs[I_4x4];
    h->stat.frame.i_mb_count_p = mbs[P_L0] + mbs[P_8x8];

    h->fdec->f_qp_avg_rc = rc->qpa_rc /= h->mb.i_mb_count;
    h->fdec->f_qp_avg_aq = (float)rc->qpa_aq / h->mb.i_mb_count;
    h->fdec->f_crf_avg = h->param.rc.f_rf_constant + h->fdec->f_qp_avg_rc - rc->qp_novbv;

    if( rc->b_abr )
    {
        rc->cplxr_sum += bits * qp2qscale( rc->qpa_rc ) / rc->last_rceq;

        rc->cplxr_sum *= rc->cbr_decay;
        rc->wanted_bits_window += h->fenc->f_duration * rc->bitrate;
        rc->wanted_bits_window *= rc->cbr_decay;
    }

    *filler = update_vbv( h, bits );
    rc->filler_bits_sum += *filler * 8;

    if( h->sps->vui.b_nal_hrd_parameters_present )
    {
        if( h->fenc->i_frame == 0 )
        {
            // access unit initialises the HRD
            h->fenc->hrd_timing.cpb_initial_arrival_time = 0;
            rc->initial_cpb_removal_delay = h->initial_cpb_removal_delay;
            rc->initial_cpb_removal_delay_offset = h->initial_cpb_removal_delay_offset;
            h->fenc->hrd_timing.cpb_removal_time = rc->nrt_first_access_unit = (double)rc->initial_cpb_removal_delay / 90000;
        }
        else
        {
            h->fenc->hrd_timing.cpb_removal_time = rc->nrt_first_access_unit + (double)(h->fenc->i_cpb_delay - h->i_cpb_delay_pir_offset) *
                                                   h->sps->vui.i_num_units_in_tick / h->sps->vui.i_time_scale;

            double cpb_earliest_arrival_time = h->fenc->hrd_timing.cpb_removal_time - (double)rc->initial_cpb_removal_delay / 90000;
            if( h->fenc->b_keyframe )
            {
                 rc->nrt_first_access_unit = h->fenc->hrd_timing.cpb_removal_time;
                 rc->initial_cpb_removal_delay = h->initial_cpb_removal_delay;
                 rc->initial_cpb_removal_delay_offset = h->initial_cpb_removal_delay_offset;
            }
            else
                 cpb_earliest_arrival_time -= (double)rc->initial_cpb_removal_delay_offset / 90000;

            if( h->sps->vui.hrd.b_cbr_hrd )
                h->fenc->hrd_timing.cpb_initial_arrival_time = rc->previous_cpb_final_arrival_time;
            else
                h->fenc->hrd_timing.cpb_initial_arrival_time = X264_MAX( rc->previous_cpb_final_arrival_time, cpb_earliest_arrival_time );
        }
        int filler_bits = *filler ? X264_MAX( (FILLER_OVERHEAD - h->param.b_annexb), *filler )*8 : 0;
        // Equation C-6
        h->fenc->hrd_timing.cpb_final_arrival_time = rc->previous_cpb_final_arrival_time = h->fenc->hrd_timing.cpb_initial_arrival_time +
                                                     (double)(bits + filler_bits) / h->sps->vui.hrd.i_bit_rate_unscaled;

        h->fenc->hrd_timing.dpb_output_time = (double)h->fenc->i_dpb_output_delay * h->sps->vui.i_num_units_in_tick / h->sps->vui.i_time_scale +
                                              h->fenc->hrd_timing.cpb_removal_time;
    }

    return 0;
}

/****************************************************************************
 * 2 pass functions
 ***************************************************************************/

/**
 * modify the bitrate curve from pass1 for one frame
 */
static double get_qscale(x264_t *h, ratecontrol_entry_t *rce, double rate_factor, int frame_num)
{
    x264_ratecontrol_t *rcc= h->rc;
    double q;
    q = pow( rce->blurred_complexity, 1 - rcc->qcompress );

    // avoid NaN's in the rc_eq
    if( !isfinite(q) || rce->tex_bits + rce->mv_bits == 0 )
        q = rcc->last_qscale_for[rce->pict_type];
    else
    {
        rcc->last_rceq = q;
        q /= rate_factor;
        rcc->last_qscale = q;
    }

    return q;
}

static float predict_size( predictor_t *p, float q, float var )
{
    return (p->coeff*var + p->offset) / (q*p->count);
}

static void update_predictor( predictor_t *p, float q, float var, float bits )
{
    float range = 1.5;
    if( var < 10 )
        return;
    float old_coeff = p->coeff / p->count;
    float new_coeff = X264_MAX( bits*q / var, p->coeff_min );
    float new_coeff_clipped = x264_clip3f( new_coeff, old_coeff/range, old_coeff*range );
    float new_offset = bits*q - new_coeff_clipped * var;
    if( new_offset >= 0 )
        new_coeff = new_coeff_clipped;
    else
        new_offset = 0;
    p->count  *= p->decay;
    p->coeff  *= p->decay;
    p->offset *= p->decay;
    p->count  ++;
    p->coeff  += new_coeff;
    p->offset += new_offset;
}

// update VBV after encoding a frame
static int update_vbv( x264_t *h, int bits )
{
    int filler = 0;
    int bitrate = h->sps->vui.hrd.i_bit_rate_unscaled;
    x264_ratecontrol_t *rcc = h->rc;
    x264_ratecontrol_t *rct = h->rc;
    uint64_t buffer_size = (uint64_t)h->sps->vui.hrd.i_cpb_size_unscaled * h->sps->vui.i_time_scale;

    if( rcc->last_satd >= h->mb.i_mb_count )
        update_predictor( &rct->pred[h->sh.i_type], qp2qscale( rcc->qpa_rc ), rcc->last_satd, bits );

    if( !rcc->b_vbv )
        return filler;

    rct->buffer_fill_final -= (uint64_t)bits * h->sps->vui.i_time_scale;

    if( rct->buffer_fill_final < 0 )
        x264_log( h, X264_LOG_WARNING, "VBV underflow (frame %d, %.0f bits)\n", h->i_frame, (double)rct->buffer_fill_final / h->sps->vui.i_time_scale );
    rct->buffer_fill_final = X264_MAX( rct->buffer_fill_final, 0 );

    rct->buffer_fill_final += (uint64_t)bitrate * h->sps->vui.i_num_units_in_tick * h->fenc->i_cpb_duration;

    if( h->sps->vui.hrd.b_cbr_hrd && rct->buffer_fill_final > buffer_size )
    {
        int64_t scale = (int64_t)h->sps->vui.i_time_scale * 8;
        filler = (rct->buffer_fill_final - buffer_size + scale - 1) / scale;
        bits = X264_MAX( (FILLER_OVERHEAD - h->param.b_annexb), filler ) * 8;
        rct->buffer_fill_final -= (uint64_t)bits * h->sps->vui.i_time_scale;
    }
    else
        rct->buffer_fill_final = X264_MIN( rct->buffer_fill_final, buffer_size );

    return filler;
}

void x264_hrd_fullness( x264_t *h )
{
    x264_ratecontrol_t *rct = h->rc;
    uint64_t denom = (uint64_t)h->sps->vui.hrd.i_bit_rate_unscaled * h->sps->vui.i_time_scale / rct->hrd_multiply_denom;
    uint64_t cpb_state = rct->buffer_fill_final;
    uint64_t cpb_size = (uint64_t)h->sps->vui.hrd.i_cpb_size_unscaled * h->sps->vui.i_time_scale;
    uint64_t multiply_factor = 180000 / rct->hrd_multiply_denom;

    if( rct->buffer_fill_final < 0 || rct->buffer_fill_final > cpb_size )
    {
         x264_log( h, X264_LOG_WARNING, "CPB %s: %.0lf bits in a %.0lf-bit buffer\n",
                   rct->buffer_fill_final < 0 ? "underflow" : "overflow", (float)rct->buffer_fill_final/denom, (float)cpb_size/denom );
    }

    h->initial_cpb_removal_delay = (multiply_factor * cpb_state + denom) / (2*denom);
    h->initial_cpb_removal_delay_offset = (multiply_factor * cpb_size + denom) / (2*denom) - h->initial_cpb_removal_delay;
}

// provisionally update VBV according to the planned size of all frames currently in progress
static void update_vbv_plan( x264_t *h, int overhead )
{
    x264_ratecontrol_t *rcc = h->rc;
    rcc->buffer_fill = h->rc->buffer_fill_final / h->sps->vui.i_time_scale;
    rcc->buffer_fill = X264_MIN( rcc->buffer_fill, rcc->buffer_size );
    rcc->buffer_fill -= overhead;
}

// apply VBV constraints and clip qscale to between lmin and lmax
static double clip_qscale( x264_t *h, int pict_type, double q )
{
    x264_ratecontrol_t *rcc = h->rc;
    double lmin = rcc->lmin[pict_type];
    double lmax = rcc->lmax[pict_type];
    if( rcc->rate_factor_max_increment )
        lmax = X264_MIN( lmax, qp2qscale( rcc->qp_novbv + rcc->rate_factor_max_increment ) );
    double q0 = q;

    /* B-frames are not directly subject to VBV,
     * since they are controlled by the P-frames' QPs. */

    if( rcc->b_vbv && rcc->last_satd > 0 )
    {
        /* Fallback to old purely-reactive algorithm: no lookahead. */
        {
            if( ( pict_type == SLICE_TYPE_P ||
                ( pict_type == SLICE_TYPE_I && rcc->last_non_b_pict_type == SLICE_TYPE_I ) ) &&
                rcc->buffer_fill/rcc->buffer_size < 0.5 )
            {
                q /= x264_clip3f( 2.0*rcc->buffer_fill/rcc->buffer_size, 0.5, 1.0 );
            }

            /* Now a hard threshold to make sure the frame fits in VBV.
             * This one is mostly for I-frames. */
            double bits = predict_size( &rcc->pred[h->sh.i_type], q, rcc->last_satd );
            double qf = 1.0;
            /* For small VBVs, allow the frame to use up the entire VBV. */
            double max_fill_factor = h->param.rc.i_vbv_buffer_size >= 5*h->param.rc.i_vbv_max_bitrate / rcc->fps ? 2 : 1;
            /* For single-frame VBVs, request that the frame use up the entire VBV. */
            double min_fill_factor = rcc->single_frame_vbv ? 1 : 2;

            if( bits > rcc->buffer_fill/max_fill_factor )
                qf = x264_clip3f( rcc->buffer_fill/(max_fill_factor*bits), 0.2, 1.0 );
            q /= qf;
            bits *= qf;
            if( bits < rcc->buffer_rate/min_fill_factor )
                q *= bits*min_fill_factor/rcc->buffer_rate;
            q = X264_MAX( q0, q );
        }

        /* Apply MinCR restrictions */
        double bits = predict_size( &rcc->pred[h->sh.i_type], q, rcc->last_satd );
        if( bits > rcc->frame_size_maximum )
            q *= bits / rcc->frame_size_maximum;
        bits = predict_size( &rcc->pred[h->sh.i_type], q, rcc->last_satd );

        /* Check B-frame complexity, and use up any bits that would
         * overflow before the next P-frame. */
        if( h->sh.i_type == SLICE_TYPE_P && !rcc->single_frame_vbv )
        {
            double pbbits = bits;
            double space;
            double minigop_cpb_duration;

            minigop_cpb_duration = h->fenc->f_planned_cpb_duration;
            space = rcc->buffer_fill + minigop_cpb_duration*rcc->vbv_max_rate - rcc->buffer_size;
            if( pbbits < space )
            {
                q *= X264_MAX( pbbits / space, bits / (0.5 * rcc->buffer_size) );
            }
            q = X264_MAX( q0/2, q );
        }

        if( !rcc->b_vbv_min_rate )
            q = X264_MAX( q0, q );
    }

    if( lmin==lmax )
        return lmin;
    else
        return x264_clip3f( q, lmin, lmax );
}

// update qscale for 1 frame based on actual bits used so far
static float rate_estimate_qscale( x264_t *h )
{
    float q;
    x264_ratecontrol_t *rcc = h->rc;
    ratecontrol_entry_t UNINIT(rce);
    int pict_type = h->sh.i_type;
    int64_t total_bits = 8*(h->stat.i_frame_size[SLICE_TYPE_I]
                          + h->stat.i_frame_size[SLICE_TYPE_P])
                       - rcc->filler_bits_sum;

    {
        double abr_buffer = 2 * rcc->rate_tolerance * rcc->bitrate;

        /* 1pass ABR */
        {
            /* Calculate the quantizer which would have produced the desired
             * average bitrate if it had been applied to all frames so far.
             * Then modulate that quant based on the current frame's complexity
             * relative to the average complexity so far (using the 2pass RCEQ).
             * Then bias the quant up or down if total size so far was far from
             * the target.
             * Result: Depending on the value of rate_tolerance, there is a
             * tradeoff between quality and bitrate precision. But at large
             * tolerances, the bit distribution approaches that of 2pass. */

            double wanted_bits, overflow = 1;

            rcc->last_satd = x264_rc_analyse_slice( h );
            rcc->short_term_cplxsum *= 0.5;
            rcc->short_term_cplxcount *= 0.5;
            rcc->short_term_cplxsum += rcc->last_satd / (CLIP_DURATION(h->fenc->f_duration) / BASE_FRAME_DURATION);
            rcc->short_term_cplxcount ++;

            rce.tex_bits = rcc->last_satd;
            rce.blurred_complexity = rcc->short_term_cplxsum / rcc->short_term_cplxcount;
            rce.mv_bits = 0;
            rce.p_count = rcc->nmb;
            rce.i_count = 0;
            rce.s_count = 0;
            rce.qscale = 1;
            rce.pict_type = pict_type;
            rce.i_duration = h->fenc->i_duration;

            if( h->param.rc.i_rc_method == X264_RC_CRF )
            {
                q = get_qscale( h, &rce, rcc->rate_factor_constant, h->fenc->i_frame );
            }
            else
            {
                q = get_qscale( h, &rce, rcc->wanted_bits_window / rcc->cplxr_sum, h->fenc->i_frame );

                /* ABR code can potentially be counterproductive in CBR, so just don't bother.
                 * Don't run it if the frame complexity is zero either. */
                if( !rcc->b_vbv_min_rate && rcc->last_satd )
                {
                    // FIXME is it simpler to keep track of wanted_bits in ratecontrol_end?
                    int i_frame_done = h->i_frame;
                    double time_done = i_frame_done / rcc->fps;
                    wanted_bits = time_done * rcc->bitrate;
                    if( wanted_bits > 0 )
                    {
                        abr_buffer *= X264_MAX( 1, sqrt( time_done ) );
                        overflow = x264_clip3f( 1.0 + (total_bits - wanted_bits) / abr_buffer, .5, 2 );
                        q *= overflow;
                    }
                }
            }

            if( pict_type == SLICE_TYPE_I && h->param.i_keyint_max > 1
                /* should test _next_ pict type, but that isn't decided yet */
                && rcc->last_non_b_pict_type != SLICE_TYPE_I )
            {
                q = qp2qscale( rcc->accum_p_qp / rcc->accum_p_norm );
                q /= fabs( h->param.rc.f_ip_factor );
            }
            else if( h->i_frame > 0 )
            {
                if( h->param.rc.i_rc_method != X264_RC_CRF )
                {
                    /* Asymmetric clipping, because symmetric would prevent
                     * overflow control in areas of rapidly oscillating complexity */
                    double lmin = rcc->last_qscale_for[pict_type] / rcc->lstep;
                    double lmax = rcc->last_qscale_for[pict_type] * rcc->lstep;
                    if( overflow > 1.1 && h->i_frame > 3 )
                        lmax *= rcc->lstep;
                    else if( overflow < 0.9 )
                        lmin /= rcc->lstep;

                    q = x264_clip3f(q, lmin, lmax);
                }
            }
            else if( h->param.rc.i_rc_method == X264_RC_CRF && rcc->qcompress != 1 )
            {
                q = qp2qscale( ABR_INIT_QP ) / fabs( h->param.rc.f_ip_factor );
            }
            rcc->qp_novbv = qscale2qp( q );

            //FIXME use get_diff_limited_q() ?
            q = clip_qscale( h, pict_type, q );
        }

        rcc->last_qscale_for[pict_type] =
        rcc->last_qscale = q;

        if( h->fenc->i_frame == 0 )
            rcc->last_qscale_for[SLICE_TYPE_P] = q * fabs( h->param.rc.f_ip_factor );

        rcc->frame_size_planned = predict_size( &rcc->pred[h->sh.i_type], q, rcc->last_satd );

        /* Always use up the whole VBV in this case. */
        if( rcc->single_frame_vbv )
            rcc->frame_size_planned = rcc->buffer_rate;
        /* Limit planned size by MinCR */
        if( rcc->b_vbv )
            rcc->frame_size_planned = X264_MIN( rcc->frame_size_planned, rcc->frame_size_maximum );
        h->rc->frame_size_estimated = rcc->frame_size_planned;
        return q;
    }
}



