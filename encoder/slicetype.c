/*****************************************************************************
 * slicetype.c: lookahead analysis
 *****************************************************************************
 * Copyright (C) 2005-2013 x264 project
 *
 * Authors: Jason Garrett-Glaser <darkshikari@gmail.com>
 *          Loren Merritt <lorenm@u.washington.edu>
 *          Dylan Yudaken <dyudaken@gmail.com>
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
#include "macroblock.h"
#include "me.h"

static int x264_slicetype_frame_cost( x264_t *h, x264_mb_analysis_t *a,
                                      x264_frame_t **frames, int p0 );

static void x264_lowres_context_init( x264_t *h, x264_mb_analysis_t *a )
{
    a->i_qp = X264_LOOKAHEAD_QP;
    a->i_lambda = x264_lambda_tab[ a->i_qp ];
    x264_mb_analyse_load_costs( h, a );
    if( h->param.analyse.i_subpel_refine > 1 )
    {
        h->mb.i_me_method = X264_MIN( X264_ME_HEX, h->param.analyse.i_me_method );
        h->mb.i_subpel_refine = 4;
    }
    else
    {
        h->mb.i_me_method = X264_ME_DIA;
        h->mb.i_subpel_refine = 2;
    }
    h->mb.b_chroma_me = 0;
}

/* Output buffers are separated by 128 bytes to avoid false sharing of cachelines
 * in multithreaded lookahead. */
#define PAD_SIZE 32
/* cost_est, cost_est_aq, intra_mbs, num rows */
#define NUM_INTS 4
#define COST_EST 0
#define COST_EST_AQ 1
#define INTRA_MBS 2
#define NUM_ROWS 3
#define ROW_SATD (NUM_INTS + (h->mb.i_mb_y - h->i_threadslice_start))

static void x264_slicetype_mb_cost( x264_t *h, x264_mb_analysis_t *a,
                                    x264_frame_t **frames, int p0, int b,
                                    int dist_scale_factor, int do_search,
                                    int *output_inter, int *output_intra )
{
    x264_frame_t *fref0 = frames[p0];
    x264_frame_t *fenc  = frames[b];
    const int i_mb_x = h->mb.i_mb_x;
    const int i_mb_y = h->mb.i_mb_y;
    const int i_mb_stride = h->mb.i_mb_width;
    const int i_mb_xy = i_mb_x + i_mb_y * i_mb_stride;
    const int i_stride = fenc->i_stride_lowres;
    const int i_pel_offset = 8 * (i_mb_x + i_mb_y * i_stride);
    int16_t (*fenc_mvs[2])[2] = { &fenc->lowres_mvs[0][b-p0-1][i_mb_xy], &fenc->lowres_mvs[1][-b][i_mb_xy] };
    int (*fenc_costs[2]) = { &fenc->lowres_mv_costs[0][b-p0-1][i_mb_xy], &fenc->lowres_mv_costs[1][-b][i_mb_xy] };
    int b_frame_score_mb = (i_mb_x > 0 && i_mb_x < h->mb.i_mb_width - 1 &&
                            i_mb_y > 0 && i_mb_y < h->mb.i_mb_height - 1) ||
                            h->mb.i_mb_width <= 2 || h->mb.i_mb_height <= 2;

    ALIGNED_ARRAY_16( pixel, pix1,[9*FDEC_STRIDE] );
    x264_me_t m[2];
    int i_bcost = COST_MAX;
    int list_used = 0;
    /* A small, arbitrary bias to avoid VBV problems caused by zero-residual lookahead blocks. */
    int lowres_penalty = 4;

    h->mb.pic.p_fenc[0] = h->mb.pic.fenc_buf;
    h->mc.copy[PIXEL_8x8]( h->mb.pic.p_fenc[0], FENC_STRIDE, &fenc->lowres[0][i_pel_offset], i_stride, 8 );

    if( p0 == 1 )
        goto lowres_intra_mb;

    // no need for h->mb.mv_min[]
    h->mb.mv_limit_fpel[0][0] = -8*h->mb.i_mb_x - 4;
    h->mb.mv_limit_fpel[1][0] = 8*( h->mb.i_mb_width - h->mb.i_mb_x - 1 ) + 4;
    h->mb.mv_min_spel[0] = 4*( h->mb.mv_limit_fpel[0][0] - 8 );
    h->mb.mv_max_spel[0] = 4*( h->mb.mv_limit_fpel[1][0] + 8 );
    if( h->mb.i_mb_x >= h->mb.i_mb_width - 2 )
    {
        h->mb.mv_limit_fpel[0][1] = -8*h->mb.i_mb_y - 4;
        h->mb.mv_limit_fpel[1][1] = 8*( h->mb.i_mb_height - h->mb.i_mb_y - 1 ) + 4;
        h->mb.mv_min_spel[1] = 4*( h->mb.mv_limit_fpel[0][1] - 8 );
        h->mb.mv_max_spel[1] = 4*( h->mb.mv_limit_fpel[1][1] + 8 );
    }

#define LOAD_HPELS_LUMA(dst, src) \
    { \
        (dst)[0] = &(src)[0][i_pel_offset]; \
        (dst)[1] = &(src)[1][i_pel_offset]; \
        (dst)[2] = &(src)[2][i_pel_offset]; \
        (dst)[3] = &(src)[3][i_pel_offset]; \
    }
#define LOAD_WPELS_LUMA(dst,src) \
    (dst) = &(src)[i_pel_offset];

#define CLIP_MV( mv ) \
    do{ \
        mv[0] = x264_clip3( mv[0], h->mb.mv_min_spel[0], h->mb.mv_max_spel[0] ); \
        mv[1] = x264_clip3( mv[1], h->mb.mv_min_spel[1], h->mb.mv_max_spel[1] ); \
    }while(0)

    m[0].i_pixel = PIXEL_8x8;
    m[0].p_cost_mv = a->p_cost_mv;
    m[0].i_stride[0] = i_stride;
    m[0].p_fenc[0] = h->mb.pic.p_fenc[0];
    m[0].i_ref = 0;
    LOAD_HPELS_LUMA( m[0].p_fref, fref0->lowres );


    {
        if( do_search )
        {
            int i_mvc = 0;
            int16_t (*fenc_mv)[2] = fenc_mvs[0];
            ALIGNED_4( int16_t mvc[4][2] );

            /* Reverse-order MV prediction. */
            M32( mvc[0] ) = 0;
            M32( mvc[2] ) = 0;
#define MVC(mv) { CP32( mvc[i_mvc], mv ); i_mvc++; }
            if( i_mb_x < h->mb.i_mb_width - 1 )
                MVC( fenc_mv[1] );
            if( i_mb_y < h->i_threadslice_end - 1 )
            {
                MVC( fenc_mv[i_mb_stride] );
                if( i_mb_x > 0 )
                    MVC( fenc_mv[i_mb_stride-1] );
                if( i_mb_x < h->mb.i_mb_width - 1 )
                    MVC( fenc_mv[i_mb_stride+1] );
            }
#undef MVC
            if( i_mvc <= 1 )
                CP32( m[0].mvp, mvc[0] );
            else
                x264_median_mv( m[0].mvp, mvc[0], mvc[1], mvc[2] );

            /* Fast skip for cases of near-zero residual.  Shortcut: don't bother except in the mv0 case,
             * since anything else is likely to have enough residual to not trigger the skip. */
            if( !M32( m[0].mvp ) )
            {
                m[0].cost = h->pixf.mbcmp[PIXEL_8x8]( m[0].p_fenc[0], FENC_STRIDE, m[0].p_fref[0], m[0].i_stride[0] );
                if( m[0].cost < 64 )
                {
                    M32( m[0].mv ) = 0;
                    goto skip_motionest;
                }
            }

            x264_me_search( h, &m[0], mvc, i_mvc );
            m[0].cost -= a->p_cost_mv[0]; // remove mvcost from skip mbs
            if( M32( m[0].mv ) )
                m[0].cost += 5 * a->i_lambda;

skip_motionest:
            CP32( fenc_mvs[0], m[0].mv );
            *fenc_costs[0] = m[0].cost;
        }
        else
        {
            CP32( m[0].mv, fenc_mvs[0] );
            m[0].cost = *fenc_costs[0];
        }
        COPY2_IF_LT( i_bcost, m[0].cost, list_used, 1 );
    }

lowres_intra_mb:
    if( !fenc->b_intra_calculated )
    {
        ALIGNED_ARRAY_16( pixel, edge,[36] );
        pixel *pix = &pix1[8+FDEC_STRIDE];
        pixel *src = &fenc->lowres[0][i_pel_offset];
        const int intra_penalty = 5 * a->i_lambda;
        int satds[3];
        int pixoff = 4 / sizeof(pixel);

        /* Avoid store forwarding stalls by writing larger chunks */
        memcpy( pix-FDEC_STRIDE, src-i_stride, 16 * sizeof(pixel) );
        for( int i = -1; i < 8; i++ )
            M32( &pix[i*FDEC_STRIDE-pixoff] ) = M32( &src[i*i_stride-pixoff] );

        h->pixf.intra_mbcmp_x3_8x8c( h->mb.pic.p_fenc[0], pix, satds );
        int i_icost = X264_MIN3( satds[0], satds[1], satds[2] );

        i_icost += intra_penalty + lowres_penalty;
        fenc->i_intra_cost[i_mb_xy] = i_icost;
        int i_icost_aq = i_icost;
        if( h->param.rc.i_aq_mode )
            i_icost_aq = (i_icost_aq * fenc->i_inv_qscale_factor[i_mb_xy] + 128) >> 8;
        output_intra[ROW_SATD] += i_icost_aq;
        if( b_frame_score_mb )
        {
            output_intra[COST_EST] += i_icost;
            output_intra[COST_EST_AQ] += i_icost_aq;
        }
    }
    i_bcost += lowres_penalty;

    /* forbid intra-mbs in B-frames, because it's rare and not worth checking */
    /* FIXME: Should we still forbid them now that we cache intra scores? */
    {
        int i_icost = fenc->i_intra_cost[i_mb_xy];
        int b_intra = i_icost < i_bcost;
        if( b_intra )
        {
            i_bcost = i_icost;
            list_used = 0;
        }
        if( b_frame_score_mb )
            output_inter[INTRA_MBS] += b_intra;
    }

    /* In an I-frame, we've already added the results above in the intra section. */
    if( p0 != 1 )
    {
        int i_bcost_aq = i_bcost;
        if( h->param.rc.i_aq_mode )
            i_bcost_aq = (i_bcost_aq * fenc->i_inv_qscale_factor[i_mb_xy] + 128) >> 8;
        output_inter[ROW_SATD] += i_bcost_aq;
        if( b_frame_score_mb )
        {
            /* Don't use AQ-weighted costs for slicetype decision, only for ratecontrol. */
            output_inter[COST_EST] += i_bcost;
            output_inter[COST_EST_AQ] += i_bcost_aq;
        }
    }

    fenc->lowres_costs[b-p0][1-b][i_mb_xy] = X264_MIN( i_bcost, LOWRES_COST_MASK ) + (list_used << LOWRES_COST_SHIFT);
}

#define NUM_MBS\
   (h->mb.i_mb_width > 2 && h->mb.i_mb_height > 2 ?\
   (h->mb.i_mb_width - 2) * (h->mb.i_mb_height - 2) :\
    h->mb.i_mb_width * h->mb.i_mb_height)

typedef struct
{
    x264_t *h;
    x264_mb_analysis_t *a;
    x264_frame_t **frames;
    int p0;
    int b;
    int dist_scale_factor;
    int do_search;
    int *output_inter;
    int *output_intra;
} x264_slicetype_slice_t;

static void x264_slicetype_slice_cost( x264_slicetype_slice_t *s )
{
    x264_t *h = s->h;

    /* Lowres lookahead goes backwards because the MVs are used as predictors in the main encode.
     * This considerably improves MV prediction overall. */

    /* The edge mbs seem to reduce the predictive quality of the
     * whole frame's score, but are needed for a spatial distribution. */
    int do_edges = h->param.rc.i_vbv_buffer_size || h->mb.i_mb_width <= 2 || h->mb.i_mb_height <= 2;

    int start_y = X264_MIN( h->i_threadslice_end - 1, h->mb.i_mb_height - 2 + do_edges );
    int end_y = X264_MAX( h->i_threadslice_start, 1 - do_edges );
    int start_x = h->mb.i_mb_width - 2 + do_edges;
    int end_x = 1 - do_edges;

    for( h->mb.i_mb_y = start_y; h->mb.i_mb_y >= end_y; h->mb.i_mb_y-- )
        for( h->mb.i_mb_x = start_x; h->mb.i_mb_x >= end_x; h->mb.i_mb_x-- )
            x264_slicetype_mb_cost( h, s->a, s->frames, s->p0, s->b, s->dist_scale_factor,
                                    s->do_search, s->output_inter, s->output_intra );
}

static int x264_slicetype_frame_cost( x264_t *h, x264_mb_analysis_t *a,
                                      x264_frame_t *frames[2], int p0  )
{
    int i_score = 0;
    int do_search;
    x264_frame_t *fenc = frames[1];

    /* Check whether we already evaluated this frame
     * If we have tried this frame as P, then we have also tried
     * the preceding frames as B. (is this still true?) */
    /* Also check that we already calculated the row SATDs for the current frame. */
    if( fenc->i_cost_est[1-p0] >= 0 && (!h->param.rc.i_vbv_buffer_size || fenc->i_row_satds[1-p0][0] != -1) )
        i_score = fenc->i_cost_est[1-p0];
    else
    {
        int dist_scale_factor = 128;

        /* For each list, check to see whether we have lowres motion-searched this reference frame before. */
        do_search = 1 != p0 && fenc->lowres_mvs[0][-p0][0][0] == 0x7FFF;
        if( do_search )
        {
            fenc->lowres_mvs[0][-p0][0][0] = 0;
        }

        if( 1 != p0 )
            dist_scale_factor = ( ((1-p0) << 8) + ((1-p0) >> 1) ) / (1-p0);

        int output_buf_size = h->mb.i_mb_height + (NUM_INTS + PAD_SIZE);
        int *output_inter[1];
        int *output_intra[1];
        output_inter[0] = h->scratch_buffer2;
        output_intra[0] = output_inter[0] + output_buf_size;

        {
            {
                h->i_threadslice_start = 0;
                h->i_threadslice_end = h->mb.i_mb_height;
                memset( output_inter[0], 0, (output_buf_size - PAD_SIZE) * sizeof(int) );
                memset( output_intra[0], 0, (output_buf_size - PAD_SIZE) * sizeof(int) );
                output_inter[0][NUM_ROWS] = output_intra[0][NUM_ROWS] = h->mb.i_mb_height;
                x264_slicetype_slice_t s = (x264_slicetype_slice_t){ h, a, frames, p0, 1, dist_scale_factor, do_search,
                    output_inter[0], output_intra[0] };
                x264_slicetype_slice_cost( &s );
            }

            /* Sum up accumulators */
            fenc->i_intra_mbs[1-p0] = 0;
            if( !fenc->b_intra_calculated )
            {
                fenc->i_cost_est[0] = 0;
                fenc->i_cost_est_aq[0] = 0;
            }
            fenc->i_cost_est[1-p0] = 0;
            fenc->i_cost_est_aq[1-p0] = 0;

            int *row_satd_inter = fenc->i_row_satds[1-p0];
            int *row_satd_intra = fenc->i_row_satds[0];
            {
                fenc->i_intra_mbs[1-p0] += output_inter[0][INTRA_MBS];
                if( !fenc->b_intra_calculated )
                {
                    fenc->i_cost_est[0] += output_intra[0][COST_EST];
                    fenc->i_cost_est_aq[0] += output_intra[0][COST_EST_AQ];
                }

                fenc->i_cost_est[1-p0] += output_inter[0][COST_EST];
                fenc->i_cost_est_aq[1-p0] += output_inter[0][COST_EST_AQ];

                if( h->param.rc.i_vbv_buffer_size )
                {
                    int row_count = output_inter[0][NUM_ROWS];
                    memcpy( row_satd_inter, output_inter[0] + NUM_INTS, row_count * sizeof(int) );
                    if( !fenc->b_intra_calculated )
                        memcpy( row_satd_intra, output_intra[0] + NUM_INTS, row_count * sizeof(int) );
                    row_satd_inter += row_count;
                    row_satd_intra += row_count;
                }
            }

            i_score = fenc->i_cost_est[1-p0];

            fenc->b_intra_calculated = 1;

            fenc->i_cost_est[1-p0] = i_score;

        }
    }

    return i_score;
}

static void x264_calculate_durations( x264_t *h, x264_frame_t *cur_frame, x264_frame_t *prev_frame, int64_t *i_cpb_delay, int64_t *i_coded_fields )
{
    cur_frame->i_cpb_delay = *i_cpb_delay;
    cur_frame->i_dpb_output_delay = cur_frame->i_field_cnt - *i_coded_fields;

    // add a correction term for frame reordering
    cur_frame->i_dpb_output_delay += h->sps->vui.i_num_reorder_frames*2;

    // fix possible negative dpb_output_delay because of pulldown changes and reordering
    if( cur_frame->i_dpb_output_delay < 0 )
    {
        cur_frame->i_cpb_delay += cur_frame->i_dpb_output_delay;
        cur_frame->i_dpb_output_delay = 0;
        if( prev_frame )
            prev_frame->i_cpb_duration += cur_frame->i_dpb_output_delay;
    }

    // don't reset cpb delay for IDR frames when using intra-refresh
    if( cur_frame->b_keyframe && !h->param.b_intra_refresh )
        *i_cpb_delay = 0;

    *i_cpb_delay += cur_frame->i_duration;
    *i_coded_fields += cur_frame->i_duration;
    cur_frame->i_cpb_duration = cur_frame->i_duration;
}

static int scenecut_internal( x264_t *h, x264_mb_analysis_t *a, x264_frame_t **frames, int real_scenecut )
{
    x264_frame_t *frame = frames[1];

    x264_slicetype_frame_cost( h, a, frames, 0 );

    int icost = frame->i_cost_est[0];
    int pcost = frame->i_cost_est[1];
    float f_bias;
    int i_gop_size = frame->i_frame - h->lookahead->i_last_keyframe;
    float f_thresh_max = h->param.i_scenecut_threshold / 100.0;
    /* magic numbers pulled out of thin air */
    float f_thresh_min = f_thresh_max * 0.25;
    int res;

    if( h->param.i_keyint_min == h->param.i_keyint_max )
        f_thresh_min = f_thresh_max;
    if( i_gop_size <= h->param.i_keyint_min / 4 || h->param.b_intra_refresh )
        f_bias = f_thresh_min / 4;
    else if( i_gop_size <= h->param.i_keyint_min )
        f_bias = f_thresh_min * i_gop_size / h->param.i_keyint_min;
    else
    {
        f_bias = f_thresh_min
                 + ( f_thresh_max - f_thresh_min )
                 * ( i_gop_size - h->param.i_keyint_min )
                 / ( h->param.i_keyint_max - h->param.i_keyint_min );
    }

    res = pcost >= (1.0 - f_bias) * icost;
    if( res && real_scenecut )
    {
        int imb = frame->i_intra_mbs[1];
        int pmb = NUM_MBS - imb;
        x264_log( h, X264_LOG_DEBUG, "scene cut at %d Icost:%d Pcost:%d ratio:%.4f bias:%.4f gop:%d (imb:%d pmb:%d)\n",
                  frame->i_frame,
                  icost, pcost, 1. - (double)pcost / icost,
                  f_bias, i_gop_size, imb, pmb );
    }
    return res;
}

static int scenecut( x264_t *h, x264_mb_analysis_t *a, x264_frame_t **frames, int real_scenecut )
{
    /* Ignore frames that are part of a flash, i.e. cannot be real scenecuts. */
    if( !frames[1]->b_scenecut )
        return 0;
    return scenecut_internal( h, a, frames, real_scenecut );
}

void x264_slicetype_analyse( x264_t *h, int intra_minigop )
{
    x264_mb_analysis_t a;
    x264_frame_t *frames[2] = { NULL, };
    int num_frames, keyint_limit;
    int framecnt = 0;

    int keyframe = !!intra_minigop;

    assert( h->frames.b_have_lowres );

    if( !h->lookahead->last_nonb )
        return;
    frames[0] = h->lookahead->last_nonb;
    if( h->lookahead->current_frame->i_type == X264_TYPE_AUTO ){
        framecnt++;
        frames[1] = h->lookahead->current_frame;
    }
    x264_lowres_context_init( h, &a );

    if( !framecnt )
    {
        return;
    }
    keyint_limit = h->param.i_keyint_max - frames[0]->i_frame + h->lookahead->i_last_keyframe - 1;
    num_frames = h->param.b_intra_refresh ? framecnt : X264_MIN( framecnt, keyint_limit );

    /* This is important psy-wise: if we have a non-scenecut keyframe,
     * there will be significant visual artifacts if the frames just before
     * go down in quality due to being referenced less, despite it being
     * more RD-optimal. */
    if( num_frames == 0 )
    {
        frames[1]->i_type = X264_TYPE_I;
        return;
    }

    int reset_start;
    if( h->param.i_scenecut_threshold && scenecut( h, &a, frames, 1  ) )
    {
        frames[1]->i_type = X264_TYPE_I;
        return;
    }


    {
        for( int j = 1; j <= num_frames; j++ )
            frames[j]->i_type = X264_TYPE_P;
        reset_start = !keyframe + 1;
    }

    /* Enforce keyframe limit. */
    if( !h->param.b_intra_refresh )
        for( int i = keyint_limit+1; i <= num_frames; i += h->param.i_keyint_max )
        {
            frames[i]->i_type = X264_TYPE_I;
            reset_start = X264_MIN( reset_start, i+1 );
        }

    /* Restore frametypes for all frames that haven't actually been decided yet. */
    for( int j = reset_start; j <= num_frames; j++ )
        frames[j]->i_type = X264_TYPE_AUTO;
}

void x264_slicetype_decide( x264_t *h )
{
    x264_frame_t *frames[2];
    x264_frame_t *frm;

    if( !h->lookahead->current_frame )
        return;

    if( h->param.i_scenecut_threshold )
        x264_slicetype_analyse( h, 0 );

    {
        frm = h->lookahead->current_frame;

        if( frm->i_type == X264_TYPE_KEYFRAME )
            frm->i_type = X264_TYPE_IDR;

        /* Limit GOP size */
        if( (!h->param.b_intra_refresh || frm->i_frame == 0) && frm->i_frame - h->lookahead->i_last_keyframe >= h->param.i_keyint_max )
        {
            if( frm->i_type == X264_TYPE_AUTO || frm->i_type == X264_TYPE_I )
                frm->i_type = X264_TYPE_IDR;
            int warn = frm->i_type != X264_TYPE_IDR;
            if( warn )
            {
                x264_log( h, X264_LOG_WARNING, "specified frame type (%d) at %d is not compatible with keyframe interval\n", frm->i_type, frm->i_frame );
                frm->i_type = X264_TYPE_IDR;
            }
        }
        if( frm->i_type == X264_TYPE_I && frm->i_frame - h->lookahead->i_last_keyframe >= h->param.i_keyint_min )
        {
                frm->i_type = X264_TYPE_IDR;
        }
        if( frm->i_type == X264_TYPE_IDR )
        {
            /* Close GOP */
            h->lookahead->i_last_keyframe = frm->i_frame;
            frm->b_keyframe = 1;
        }

        if( frm->i_type == X264_TYPE_AUTO )
            frm->i_type = X264_TYPE_P;
    }

    /* calculate the frame costs ahead of time for x264_rc_analyse_slice while we still have lowres */
    if( h->param.rc.i_rc_method != X264_RC_CQP )
    {
        x264_mb_analysis_t a;
        int p0;

        x264_lowres_context_init( h, &a );

        frames[0] = h->lookahead->last_nonb;
        memcpy( &frames[1], h->lookahead->current_frame, sizeof(x264_frame_t*) );
        if( IS_X264_TYPE_I( h->lookahead->current_frame->i_type ) )
            p0 = 1;
        else // P
            p0 = 0;

        x264_slicetype_frame_cost( h, &a, frames, p0 );

        if( (p0 != 1 ) && h->param.rc.i_vbv_buffer_size )
        {
            /* We need the intra costs for row SATDs. */
            x264_slicetype_frame_cost( h, &a, frames, 1 );
        }
    }

    /* shift sequence to coded order.
       use a small temporary list to avoid shifting the entire next buffer around */
    int i_coded = h->lookahead->current_frame->i_frame;

    {
        h->lookahead->current_frame->i_coded = i_coded++;

        x264_calculate_durations( h, h->lookahead->current_frame, NULL, &h->i_cpb_delay, &h->i_coded_fields );

        h->lookahead->current_frame->f_planned_cpb_duration[0] = (double)h->lookahead->current_frame->i_cpb_duration *
                                                                h->sps->vui.i_num_units_in_tick / h->sps->vui.i_time_scale;
    }
}

int x264_rc_analyse_slice( x264_t *h )
{
    int p0 = 0;
    int b;
    int cost;


    if( IS_X264_TYPE_I(h->fenc->i_type) )
        b = 0;
    else if( h->fenc->i_type == X264_TYPE_P )
        b = 1;
    else //B
    {
        assert(0);
    }
    /* We don't need to assign p0 since we are not performing any real analysis here. */
    x264_frame_t **frames = &h->fenc - b;

    /* cost should have been already calculated by x264_slicetype_decide */
    cost = frames[b]->i_cost_est[b-p0];
    assert( cost >= 0 );

    /* In AQ, use the weighted score instead. */
    if( h->param.rc.i_aq_mode )
        cost = frames[b]->i_cost_est_aq[b-p0];

    h->fenc->i_row_satd = h->fenc->i_row_satds[b-p0];
    h->fdec->i_row_satd = h->fdec->i_row_satds[b-p0];
    h->fdec->i_satd = cost;
    memcpy( h->fdec->i_row_satd, h->fenc->i_row_satd, h->mb.i_mb_height * sizeof(int) );
    if( !IS_X264_TYPE_I(h->fenc->i_type) )
        memcpy( h->fdec->i_row_satds[0], h->fenc->i_row_satds[0], h->mb.i_mb_height * sizeof(int) );

    if( h->param.b_intra_refresh && h->param.rc.i_vbv_buffer_size && h->fenc->i_type == X264_TYPE_P )
    {
        int ip_factor = 256 * h->param.rc.f_ip_factor; /* fix8 */
        for( int y = 0; y < h->mb.i_mb_height; y++ )
        {
            int mb_xy = y * h->mb.i_mb_stride + h->fdec->i_pir_start_col;
            for( int x = h->fdec->i_pir_start_col; x <= h->fdec->i_pir_end_col; x++, mb_xy++ )
            {
                int intra_cost = (h->fenc->i_intra_cost[mb_xy] * ip_factor + 128) >> 8;
                int inter_cost = h->fenc->lowres_costs[b-p0][0][mb_xy] & LOWRES_COST_MASK;
                int diff = intra_cost - inter_cost;
                if( h->param.rc.i_aq_mode )
                    h->fdec->i_row_satd[y] += (diff * frames[b]->i_inv_qscale_factor[mb_xy] + 128) >> 8;
                else
                    h->fdec->i_row_satd[y] += diff;
                cost += diff;
            }
        }
    }

    return cost;
}
