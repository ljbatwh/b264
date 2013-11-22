/*****************************************************************************
 * set: header writing
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
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

#define bs_write_ue bs_write_ue_big

void x264_sei_write( bs_t *s, uint8_t *payload, int payload_size, int payload_type )
{
    int i;

    bs_realign( s );

    for( i = 0; i <= payload_type-255; i += 255 )
        bs_write( s, 8, 255 );
    bs_write( s, 8, payload_type-i );

    for( i = 0; i <= payload_size-255; i += 255 )
        bs_write( s, 8, 255 );
    bs_write( s, 8, payload_size-i );

    for( i = 0; i < payload_size; i++ )
        bs_write( s, 8, payload[i] );

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sps_init( x264_sps_t *sps, int i_id, x264_param_t *param )
{
    sps->i_id = i_id;
    sps->i_mb_width = ( param->i_width + 15 ) / 16;
    sps->i_mb_height= ( param->i_height + 15 ) / 16;
    sps->i_chroma_format_idc = CHROMA_420;

    sps->b_qpprime_y_zero_transform_bypass = param->rc.i_rc_method == X264_RC_CQP && param->rc.i_qp_constant == 0;
    sps->i_profile_idc  = PROFILE_BASELINE;

    sps->b_constraint_set0  = sps->i_profile_idc == PROFILE_BASELINE;
    /* x264 doesn't support the features that are in Baseline and not in Main,
     * namely arbitrary_slice_order and slice_groups. */
    sps->b_constraint_set1  = sps->i_profile_idc <= PROFILE_MAIN;
    /* Never set constraint_set2, it is not necessary and not used in real world. */
    sps->b_constraint_set2  = 0;
    sps->b_constraint_set3  = 0;

    sps->i_level_idc = param->i_level_idc;
    if( param->i_level_idc == 9 && ( sps->i_profile_idc == PROFILE_BASELINE || sps->i_profile_idc == PROFILE_MAIN ) )
    {
        sps->b_constraint_set3 = 1; /* level 1b with Baseline or Main profile is signalled via constraint_set3 */
        sps->i_level_idc      = 11;
    }
    /* Intra profiles */
    if( param->i_keyint_max == 1 && sps->i_profile_idc > PROFILE_HIGH )
        sps->b_constraint_set3 = 1;

    sps->vui.i_num_reorder_frames = 0;
    /* extra slot with pyramid so that we don't have to override the
     * order of forgetting old pictures */
    sps->vui.i_max_dec_frame_buffering =
    sps->i_num_ref_frames = X264_MIN(X264_REF_MAX, X264_MAX(1 + sps->vui.i_num_reorder_frames,param->i_dpb_size));
    if( param->i_keyint_max == 1 )
    {
        sps->i_num_ref_frames = 0;
        sps->vui.i_max_dec_frame_buffering = 0;
    }

    /* number of refs + current frame */
    int max_frame_num = sps->vui.i_max_dec_frame_buffering + 1;
    /* Intra refresh cannot write a recovery time greater than max frame num-1 */
    if( param->b_intra_refresh )
    {
        int time_to_recovery = X264_MIN( sps->i_mb_width - 1, param->i_keyint_max ) - 1;
        max_frame_num = X264_MAX( max_frame_num, time_to_recovery+1 );
    }

    sps->i_log2_max_frame_num = 4;
    while( (1 << sps->i_log2_max_frame_num) <= max_frame_num )
        sps->i_log2_max_frame_num++;

    sps->i_poc_type = 2;
    if( sps->i_poc_type == 0 )
    {
        int max_delta_poc = (2) * 2;
        sps->i_log2_max_poc_lsb = 4;
        while( (1 << sps->i_log2_max_poc_lsb) <= max_delta_poc * 2 )
            sps->i_log2_max_poc_lsb++;
    }

    sps->b_vui = 1;

    sps->b_gaps_in_frame_num_value_allowed = 0;
    sps->b_frame_mbs_only = 1;
    if( !sps->b_frame_mbs_only )
        sps->i_mb_height = ( sps->i_mb_height + 1 ) & ~1;
    sps->b_mb_adaptive_frame_field = 0;
    sps->b_direct8x8_inference = 1;

    sps->crop.i_left   = param->crop_rect.i_left;
    sps->crop.i_top    = param->crop_rect.i_top;
    sps->crop.i_right  = param->crop_rect.i_right + sps->i_mb_width*16 - param->i_width;
    sps->crop.i_bottom = (param->crop_rect.i_bottom + sps->i_mb_height*16 - param->i_height) >> !sps->b_frame_mbs_only;
    sps->b_crop = sps->crop.i_left  || sps->crop.i_top ||
                  sps->crop.i_right || sps->crop.i_bottom;

    sps->vui.b_aspect_ratio_info_present = 0;
    if( param->vui.i_sar_width > 0 && param->vui.i_sar_height > 0 )
    {
        sps->vui.b_aspect_ratio_info_present = 1;
        sps->vui.i_sar_width = param->vui.i_sar_width;
        sps->vui.i_sar_height= param->vui.i_sar_height;
    }

    sps->vui.b_overscan_info_present = param->vui.i_overscan > 0 && param->vui.i_overscan <= 2;
    if( sps->vui.b_overscan_info_present )
        sps->vui.b_overscan_info = ( param->vui.i_overscan == 2 ? 1 : 0 );

    sps->vui.b_signal_type_present = 0;
    sps->vui.i_vidformat = ( param->vui.i_vidformat >= 0 && param->vui.i_vidformat <= 5 ? param->vui.i_vidformat : 5 );
    sps->vui.b_fullrange = ( param->vui.b_fullrange >= 0 && param->vui.b_fullrange <= 1 ? param->vui.b_fullrange :
                           ( 0 ) );
    sps->vui.b_color_description_present = 0;

    sps->vui.i_colorprim = ( param->vui.i_colorprim >= 0 && param->vui.i_colorprim <=  9 ? param->vui.i_colorprim : 2 );
    sps->vui.i_transfer  = ( param->vui.i_transfer  >= 0 && param->vui.i_transfer  <= 15 ? param->vui.i_transfer  : 2 );
    sps->vui.i_colmatrix = ( param->vui.i_colmatrix >= 0 && param->vui.i_colmatrix <= 10 ? param->vui.i_colmatrix :
                           ( 2 ) );
    if( sps->vui.i_colorprim != 2 ||
        sps->vui.i_transfer  != 2 ||
        sps->vui.i_colmatrix != 2 )
    {
        sps->vui.b_color_description_present = 1;
    }

    if( sps->vui.i_vidformat != 5 ||
        sps->vui.b_fullrange ||
        sps->vui.b_color_description_present )
    {
        sps->vui.b_signal_type_present = 1;
    }

    /* FIXME: not sufficient for interlaced video */
    sps->vui.b_chroma_loc_info_present = param->vui.i_chroma_loc > 0 && param->vui.i_chroma_loc <= 5;
    if( sps->vui.b_chroma_loc_info_present )
    {
        sps->vui.i_chroma_loc_top = param->vui.i_chroma_loc;
        sps->vui.i_chroma_loc_bottom = param->vui.i_chroma_loc;
    }

    sps->vui.b_timing_info_present = 0;
    sps->vui.b_vcl_hrd_parameters_present = 0; // we don't support VCL HRD
    sps->vui.b_nal_hrd_parameters_present = !!param->i_nal_hrd;
    sps->vui.b_pic_struct_present = 0;

    // NOTE: HRD related parts of the SPS are initialised in x264_ratecontrol_init_reconfigurable

    sps->vui.b_bitstream_restriction = param->i_keyint_max > 1;
    if( sps->vui.b_bitstream_restriction )
    {
        sps->vui.b_motion_vectors_over_pic_boundaries = 1;
        sps->vui.i_max_bytes_per_pic_denom = 0;
        sps->vui.i_max_bits_per_mb_denom = 0;
        sps->vui.i_log2_max_mv_length_horizontal =
        sps->vui.i_log2_max_mv_length_vertical = (int)log2f( X264_MAX( 1, param->analyse.i_mv_range*4-1 ) ) + 1;
    }
}

void x264_sps_write( bs_t *s, x264_sps_t *sps )
{
    bs_realign( s );
    bs_write( s, 8, sps->i_profile_idc );
    bs_write1( s, sps->b_constraint_set0 );
    bs_write1( s, sps->b_constraint_set1 );
    bs_write1( s, sps->b_constraint_set2 );
    bs_write1( s, sps->b_constraint_set3 );

    bs_write( s, 4, 0 );    /* reserved */

    bs_write( s, 8, sps->i_level_idc );

    bs_write_ue( s, sps->i_id );

    if( sps->i_profile_idc >= PROFILE_HIGH )
    {
        bs_write_ue( s, sps->i_chroma_format_idc );
        bs_write_ue( s, BIT_DEPTH-8 ); // bit_depth_luma_minus8
        bs_write_ue( s, BIT_DEPTH-8 ); // bit_depth_chroma_minus8
        bs_write1( s, sps->b_qpprime_y_zero_transform_bypass );
        bs_write1( s, 0 ); // seq_scaling_matrix_present_flag
    }

    bs_write_ue( s, sps->i_log2_max_frame_num - 4 );
    bs_write_ue( s, sps->i_poc_type );
    if( sps->i_poc_type == 0 )
        bs_write_ue( s, sps->i_log2_max_poc_lsb - 4 );
    bs_write_ue( s, sps->i_num_ref_frames );
    bs_write1( s, sps->b_gaps_in_frame_num_value_allowed );
    bs_write_ue( s, sps->i_mb_width - 1 );
    bs_write_ue( s, (sps->i_mb_height >> !sps->b_frame_mbs_only) - 1);
    bs_write1( s, sps->b_frame_mbs_only );
    if( !sps->b_frame_mbs_only )
        bs_write1( s, sps->b_mb_adaptive_frame_field );
    bs_write1( s, sps->b_direct8x8_inference );

    bs_write1( s, sps->b_crop );
    if( sps->b_crop )
    {
        int h_shift = sps->i_chroma_format_idc == CHROMA_420;
        int v_shift = sps->i_chroma_format_idc == CHROMA_420;
        bs_write_ue( s, sps->crop.i_left   >> h_shift );
        bs_write_ue( s, sps->crop.i_right  >> h_shift );
        bs_write_ue( s, sps->crop.i_top    >> v_shift );
        bs_write_ue( s, sps->crop.i_bottom >> v_shift );
    }

    bs_write1( s, sps->b_vui );
    if( sps->b_vui )
    {
        bs_write1( s, sps->vui.b_aspect_ratio_info_present );
        if( sps->vui.b_aspect_ratio_info_present )
        {
            int i;
            static const struct { uint8_t w, h, sar; } sar[] =
            {
                // aspect_ratio_idc = 0 -> unspecified
                {  1,  1, 1 }, { 12, 11, 2 }, { 10, 11, 3 }, { 16, 11, 4 },
                { 40, 33, 5 }, { 24, 11, 6 }, { 20, 11, 7 }, { 32, 11, 8 },
                { 80, 33, 9 }, { 18, 11, 10}, { 15, 11, 11}, { 64, 33, 12},
                {160, 99, 13}, {  4,  3, 14}, {  3,  2, 15}, {  2,  1, 16},
                // aspect_ratio_idc = [17..254] -> reserved
                { 0, 0, 255 }
            };
            for( i = 0; sar[i].sar != 255; i++ )
            {
                if( sar[i].w == sps->vui.i_sar_width &&
                    sar[i].h == sps->vui.i_sar_height )
                    break;
            }
            bs_write( s, 8, sar[i].sar );
            if( sar[i].sar == 255 ) /* aspect_ratio_idc (extended) */
            {
                bs_write( s, 16, sps->vui.i_sar_width );
                bs_write( s, 16, sps->vui.i_sar_height );
            }
        }

        bs_write1( s, sps->vui.b_overscan_info_present );
        if( sps->vui.b_overscan_info_present )
            bs_write1( s, sps->vui.b_overscan_info );

        bs_write1( s, sps->vui.b_signal_type_present );
        if( sps->vui.b_signal_type_present )
        {
            bs_write( s, 3, sps->vui.i_vidformat );
            bs_write1( s, sps->vui.b_fullrange );
            bs_write1( s, sps->vui.b_color_description_present );
            if( sps->vui.b_color_description_present )
            {
                bs_write( s, 8, sps->vui.i_colorprim );
                bs_write( s, 8, sps->vui.i_transfer );
                bs_write( s, 8, sps->vui.i_colmatrix );
            }
        }

        bs_write1( s, sps->vui.b_chroma_loc_info_present );
        if( sps->vui.b_chroma_loc_info_present )
        {
            bs_write_ue( s, sps->vui.i_chroma_loc_top );
            bs_write_ue( s, sps->vui.i_chroma_loc_bottom );
        }

        bs_write1( s, sps->vui.b_timing_info_present );//b_timing_info_present is 0

        bs_write1( s, sps->vui.b_nal_hrd_parameters_present );
        if( sps->vui.b_nal_hrd_parameters_present )
        {
            bs_write_ue( s, sps->vui.hrd.i_cpb_cnt - 1 );
            bs_write( s, 4, sps->vui.hrd.i_bit_rate_scale );
            bs_write( s, 4, sps->vui.hrd.i_cpb_size_scale );

            bs_write_ue( s, sps->vui.hrd.i_bit_rate_value - 1 );
            bs_write_ue( s, sps->vui.hrd.i_cpb_size_value - 1 );

            bs_write1( s, sps->vui.hrd.b_cbr_hrd );

            bs_write( s, 5, sps->vui.hrd.i_initial_cpb_removal_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_cpb_removal_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_dpb_output_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_time_offset_length );
        }

        bs_write1( s, sps->vui.b_vcl_hrd_parameters_present );

        if( sps->vui.b_nal_hrd_parameters_present || sps->vui.b_vcl_hrd_parameters_present )
            bs_write1( s, 0 );   /* low_delay_hrd_flag */

        bs_write1( s, sps->vui.b_pic_struct_present );
        bs_write1( s, sps->vui.b_bitstream_restriction );
        if( sps->vui.b_bitstream_restriction )
        {
            bs_write1( s, sps->vui.b_motion_vectors_over_pic_boundaries );
            bs_write_ue( s, sps->vui.i_max_bytes_per_pic_denom );
            bs_write_ue( s, sps->vui.i_max_bits_per_mb_denom );
            bs_write_ue( s, sps->vui.i_log2_max_mv_length_horizontal );
            bs_write_ue( s, sps->vui.i_log2_max_mv_length_vertical );
            bs_write_ue( s, sps->vui.i_num_reorder_frames );
            bs_write_ue( s, sps->vui.i_max_dec_frame_buffering );
        }
    }

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_pps_init( x264_pps_t *pps, int i_id, x264_param_t *param, x264_sps_t *sps )
{
    pps->i_id = i_id;
    pps->i_sps_id = sps->i_id;
    pps->b_cabac = 0;

    pps->b_pic_order = 0;
    pps->i_num_slice_groups = 1;

    pps->i_num_ref_idx_l0_default_active = 1;
    pps->i_num_ref_idx_l1_default_active = 1;

    pps->b_weighted_pred = 0;
    pps->b_weighted_bipred = 0;

    pps->i_pic_init_qp = param->rc.i_rc_method == X264_RC_ABR || 26  ;
    pps->i_pic_init_qs = 26 ;

    pps->i_chroma_qp_index_offset = param->analyse.i_chroma_qp_offset;
    pps->b_deblocking_filter_control = 1;
    pps->b_constrained_intra_pred = param->b_constrained_intra;
    pps->b_redundant_pic_cnt = 0;

    pps->b_transform_8x8_mode =  0;

    for( int i = 0; i < 8; i++ ){
        pps->scaling_list[i] = x264_cqm_flat16;
    }
}

void x264_pps_write( bs_t *s, x264_sps_t *sps, x264_pps_t *pps )
{
    bs_realign( s );
    bs_write_ue( s, pps->i_id );
    bs_write_ue( s, pps->i_sps_id );

    bs_write1( s, pps->b_cabac );
    bs_write1( s, pps->b_pic_order );
    bs_write_ue( s, pps->i_num_slice_groups - 1 );

    bs_write_ue( s, pps->i_num_ref_idx_l0_default_active - 1 );
    bs_write_ue( s, pps->i_num_ref_idx_l1_default_active - 1 );
    bs_write1( s, pps->b_weighted_pred );
    bs_write( s, 2, pps->b_weighted_bipred );

    bs_write_se( s, pps->i_pic_init_qp - 26  );
    bs_write_se( s, pps->i_pic_init_qs - 26  );
    bs_write_se( s, pps->i_chroma_qp_index_offset );

    bs_write1( s, pps->b_deblocking_filter_control );
    bs_write1( s, pps->b_constrained_intra_pred );
    bs_write1( s, pps->b_redundant_pic_cnt );

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sei_recovery_point_write( x264_t *h, bs_t *s, int recovery_frame_cnt )
{
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    bs_write_ue( &q, recovery_frame_cnt ); // recovery_frame_cnt
    bs_write1( &q, 1 );   //exact_match_flag 1
    bs_write1( &q, 0 );   //broken_link_flag 0
    bs_write( &q, 2, 0 ); //changing_slice_group 0

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_RECOVERY_POINT );
}

int x264_sei_version_write( x264_t *h, bs_t *s )
{
    // random ID number generated according to ISO-11578
    static const uint8_t uuid[16] =
    {
        0xdc, 0x45, 0xe9, 0xbd, 0xe6, 0xd9, 0x48, 0xb7,
        0x96, 0x2c, 0xd8, 0x20, 0xd9, 0x23, 0xee, 0xef
    };
    char *opts = x264_param2string( &h->param, 0 );
    char *payload;
    int length;

    if( !opts )
        return -1;
    CHECKED_MALLOC( payload, 200 + strlen( opts ) );

    memcpy( payload, uuid, 16 );
    sprintf( payload+16, "x264 - core %d%s - H.264/MPEG-4 AVC codec - "
             "Copy%s 2003-2013 - http://www.videolan.org/x264.html - options: %s",
             X264_BUILD, X264_VERSION, HAVE_GPL?"left":"right", opts );
    length = strlen(payload)+1;

    x264_sei_write( s, (uint8_t *)payload, length, SEI_USER_DATA_UNREGISTERED );

    x264_free( opts );
    x264_free( payload );
    return 0;
fail:
    x264_free( opts );
    return -1;
}

void x264_sei_buffering_period_write( x264_t *h, bs_t *s )
{
    x264_sps_t *sps = h->sps;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );
    bs_write_ue( &q, sps->i_id );

    if( sps->vui.b_nal_hrd_parameters_present )
    {
        bs_write( &q, sps->vui.hrd.i_initial_cpb_removal_delay_length, h->initial_cpb_removal_delay );
        bs_write( &q, sps->vui.hrd.i_initial_cpb_removal_delay_length, h->initial_cpb_removal_delay_offset );
    }

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_BUFFERING_PERIOD );
}

void x264_sei_pic_timing_write( x264_t *h, bs_t *s )
{
    x264_sps_t *sps = h->sps;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    if( sps->vui.b_nal_hrd_parameters_present || sps->vui.b_vcl_hrd_parameters_present )
    {
        bs_write( &q, sps->vui.hrd.i_cpb_removal_delay_length, h->fenc->i_cpb_delay - h->i_cpb_delay_pir_offset );
        bs_write( &q, sps->vui.hrd.i_dpb_output_delay_length, h->fenc->i_dpb_output_delay );
    }

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_PIC_TIMING );
}


void x264_filler_write( x264_t *h, bs_t *s, int filler )
{
    bs_realign( s );

    for( int i = 0; i < filler; i++ )
        bs_write( s, 8, 0xff );

    bs_rbsp_trailing( s );
    bs_flush( s );
}





const x264_level_t x264_levels[] =
{
    //level_idc, mbps, frame_size, dpb, bitrate, cpb, mv_range, mvs_per_2mb, slice_rate,mincr, bipred8x8, direct8x8, frame_only
    { 10,    1485,    99,    396,     64,    175,  64, 64,  0, 2, 0, 0, 1 },
    {  9,    1485,    99,    396,    128,    350,  64, 64,  0, 2, 0, 0, 1 }, /* "1b" */
    { 11,    3000,   396,    900,    192,    500, 128, 64,  0, 2, 0, 0, 1 },
    { 12,    6000,   396,   2376,    384,   1000, 128, 64,  0, 2, 0, 0, 1 },
    { 13,   11880,   396,   2376,    768,   2000, 128, 64,  0, 2, 0, 0, 1 },
    { 20,   11880,   396,   2376,   2000,   2000, 128, 64,  0, 2, 0, 0, 1 },
    { 21,   19800,   792,   4752,   4000,   4000, 256, 64,  0, 2, 0, 0, 0 },
    { 22,   20250,  1620,   8100,   4000,   4000, 256, 64,  0, 2, 0, 0, 0 },
    { 30,   40500,  1620,   8100,  10000,  10000, 256, 32, 22, 2, 0, 1, 0 },
    { 31,  108000,  3600,  18000,  14000,  14000, 512, 16, 60, 4, 1, 1, 0 },
    { 32,  216000,  5120,  20480,  20000,  20000, 512, 16, 60, 4, 1, 1, 0 },
    { 40,  245760,  8192,  32768,  20000,  25000, 512, 16, 60, 4, 1, 1, 0 },
    { 41,  245760,  8192,  32768,  50000,  62500, 512, 16, 24, 2, 1, 1, 0 },
    { 42,  522240,  8704,  34816,  50000,  62500, 512, 16, 24, 2, 1, 1, 1 },
    { 50,  589824, 22080, 110400, 135000, 135000, 512, 16, 24, 2, 1, 1, 1 },
    { 51,  983040, 36864, 184320, 240000, 240000, 512, 16, 24, 2, 1, 1, 1 },
    { 52, 2073600, 36864, 184320, 240000, 240000, 512, 16, 24, 2, 1, 1, 1 },
    { 0 }
};

#define ERROR(...)\
{\
    if( verbose )\
        x264_log( h, X264_LOG_WARNING, __VA_ARGS__ );\
    ret = 1;\
}

int x264_validate_levels( x264_t *h , int verbose)
{
    int ret = 0;
    int mbs = h->sps->i_mb_width * h->sps->i_mb_height;
    int dpb = mbs * h->sps->vui.i_max_dec_frame_buffering;
    int cbp_factor = 4;

    const x264_level_t *l = x264_levels;
    while( l->level_idc != 0 && l->level_idc != h->param.i_level_idc )
        l++;

    if( l->frame_size < mbs
        || l->frame_size*8 < h->sps->i_mb_width * h->sps->i_mb_width
        || l->frame_size*8 < h->sps->i_mb_height * h->sps->i_mb_height )
        ERROR( "frame MB size (%dx%d) > level limit (%d)\n",
               h->sps->i_mb_width, h->sps->i_mb_height, l->frame_size );
    if( dpb > l->dpb )
        ERROR( "DPB size (%d frames, %d mbs) > level limit (%d frames, %d mbs)\n",
                h->sps->vui.i_max_dec_frame_buffering, dpb, l->dpb / mbs, l->dpb );

#define CHECK( name, limit, val ) \
    if( (val) > (limit) ) \
        ERROR( name " (%"PRId64") > level limit (%d)\n", (int64_t)(val), (limit) );

    CHECK( "VBV bitrate", (l->bitrate * cbp_factor) / 4, h->param.rc.i_vbv_max_bitrate );
    CHECK( "VBV buffer", (l->cpb * cbp_factor) / 4, h->param.rc.i_vbv_buffer_size );
    CHECK( "MV range", l->mv_range, h->param.analyse.i_mv_range );

    if( h->param.i_fps_den > 0 )
        CHECK( "MB rate", l->mbps, (int64_t)mbs * h->param.i_fps_num / h->param.i_fps_den );

    /* TODO check the rest of the limits */
    return ret;
}
