/*****************************************************************************
 * analyse.c: macroblock analysis
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

#define _ISOC99_SOURCE

#include "common/common.h"
#include "macroblock.h"
#include "me.h"
#include "ratecontrol.h"
#include "analyse.h"
#include "rdo.c"

typedef struct
{
    /* 16x16 */
    int       i_rd16x16;
    x264_me_t me16x16;

    /* 8x8 */
    int       i_cost8x8;
    /* [ref][0] is 16x16 mv, [ref][1..4] are 8x8 mv from partition [0..3] */
    int16_t mvc[32][5][2];
    x264_me_t me8x8[4];

    /* Sub 4x4 */
    int       i_cost4x4[4]; /* cost per 8x8 partition */
    x264_me_t me4x4[4][4];

    /* Sub 8x4 */
    int       i_cost8x4[4]; /* cost per 8x8 partition */
    x264_me_t me8x4[4][2];

    /* Sub 4x8 */
    int       i_cost4x8[4]; /* cost per 8x8 partition */
    x264_me_t me4x8[4][2];

    /* 16x8 */
    int       i_cost16x8;
    x264_me_t me16x8[2];

    /* 8x16 */
    int       i_cost8x16;
    x264_me_t me8x16[2];

} x264_mb_analysis_list_t;

typedef struct
{
    /* conduct the analysis using this lamda and QP */
    int i_lambda;
    int i_lambda2;
    int i_qp;
    uint16_t *p_cost_mv;
    uint16_t *p_cost_ref;
    int i_mbrd;


    /* I: Intra part */
    /* Take some shortcuts in intra search if intra is deemed unlikely */
    int b_fast_intra;
    int b_force_intra; /* For Periodic Intra Refresh.  Only supported in P-frames. */
    int b_avoid_topright; /* For Periodic Intra Refresh: don't predict from top-right pixels. */
    int b_try_skip;

    /* Luma part */
    int i_satd_i16x16;  //this mb Luma satd
    int i_satd_i16x16_dir[7];//7-intra16x16_pred_e
    int i_predict16x16; //the mode of i_satd_i16x16

    int i_satd_i4x4;
    int i_predict4x4[16];//16 ?

    int i_satd_pcm;

    /* Chroma part */
    int i_satd_chroma; //this mb chroma satd
    int i_satd_chroma_dir[7]; //7-intra_chroma_pred_e
    int i_predict8x8chroma; //which predict mode, one from intra_chroma_pred_e

    /* II: Inter part P frame */
    x264_mb_analysis_list_t l0;

    int i_satd8x8[4]; /* L0 [8x8 0..3] SATD only */
    int i_cost_est16x8[2]; /* Per-partition estimated cost */
    int i_cost_est8x16[2];

    int b_early_terminate;

} x264_mb_analysis_t;

/* lambda = pow(2,qp/6-2) */
const uint16_t x264_lambda_tab[QP_MAX_MAX+1] =
{
   1,   1,   1,   1,   1,   1,   1,   1, /*  0- 7 */
   1,   1,   1,   1,   1,   1,   1,   1, /*  8-15 */
   2,   2,   2,   2,   3,   3,   3,   4, /* 16-23 */
   4,   4,   5,   6,   6,   7,   8,   9, /* 24-31 */
  10,  11,  13,  14,  16,  18,  20,  23, /* 32-39 */
  25,  29,  32,  36,  40,  45,  51,  57, /* 40-47 */
  64,  72,  81,  91, 102, 114, 128, 144, /* 48-55 */
 161, 181, 203, 228, 256, 287, 323, 362, /* 56-63 */
 406, 456, 512, 575, 645, 724,           /* 64-69 */
};

/* lambda2 = pow(lambda,2) * .9 * 256 */
/* Capped to avoid overflow */
const int x264_lambda2_tab[QP_MAX_MAX+1] =
{
       14,       18,       22,       28,       36,       45,      57,      72, /*  0- 7 */
       91,      115,      145,      182,      230,      290,     365,     460, /*  8-15 */
      580,      731,      921,     1161,     1462,     1843,    2322,    2925, /* 16-23 */
     3686,     4644,     5851,     7372,     9289,    11703,   14745,   18578, /* 24-31 */
    23407,    29491,    37156,    46814,    58982,    74313,   93628,  117964, /* 32-39 */
   148626,   187257,   235929,   297252,   374514,   471859,  594505,  749029, /* 40-47 */
   943718,  1189010,  1498059,  1887436,  2378021,  2996119, 3774873, 4756042, /* 48-55 */
  5992238,  7549747,  9512085, 11984476, 15099494, 19024170,23968953,30198988, /* 56-63 */
 38048341, 47937906, 60397977, 76096683, 95875813,120795955,                   /* 64-69 */
};
//for x264_exp2fix8
const uint8_t x264_exp2_lut[64] =
{
      0,   3,   6,   8,  11,  14,  17,  20,  23,  26,  29,  32,  36,  39,  42,  45,
     48,  52,  55,  58,  62,  65,  69,  72,  76,  80,  83,  87,  91,  94,  98, 102,
    106, 110, 114, 118, 122, 126, 130, 135, 139, 143, 147, 152, 156, 161, 165, 170,
    175, 179, 184, 189, 194, 198, 203, 208, 214, 219, 224, 229, 234, 240, 245, 250
};
//for x264_log2
const float x264_log2_lut[128] =
{
    0.00000, 0.01123, 0.02237, 0.03342, 0.04439, 0.05528, 0.06609, 0.07682,
    0.08746, 0.09803, 0.10852, 0.11894, 0.12928, 0.13955, 0.14975, 0.15987,
    0.16993, 0.17991, 0.18982, 0.19967, 0.20945, 0.21917, 0.22882, 0.23840,
    0.24793, 0.25739, 0.26679, 0.27612, 0.28540, 0.29462, 0.30378, 0.31288,
    0.32193, 0.33092, 0.33985, 0.34873, 0.35755, 0.36632, 0.37504, 0.38370,
    0.39232, 0.40088, 0.40939, 0.41785, 0.42626, 0.43463, 0.44294, 0.45121,
    0.45943, 0.46761, 0.47573, 0.48382, 0.49185, 0.49985, 0.50779, 0.51570,
    0.52356, 0.53138, 0.53916, 0.54689, 0.55459, 0.56224, 0.56986, 0.57743,
    0.58496, 0.59246, 0.59991, 0.60733, 0.61471, 0.62205, 0.62936, 0.63662,
    0.64386, 0.65105, 0.65821, 0.66534, 0.67243, 0.67948, 0.68650, 0.69349,
    0.70044, 0.70736, 0.71425, 0.72110, 0.72792, 0.73471, 0.74147, 0.74819,
    0.75489, 0.76155, 0.76818, 0.77479, 0.78136, 0.78790, 0.79442, 0.80090,
    0.80735, 0.81378, 0.82018, 0.82655, 0.83289, 0.83920, 0.84549, 0.85175,
    0.85798, 0.86419, 0.87036, 0.87652, 0.88264, 0.88874, 0.89482, 0.90087,
    0.90689, 0.91289, 0.91886, 0.92481, 0.93074, 0.93664, 0.94251, 0.94837,
    0.95420, 0.96000, 0.96578, 0.97154, 0.97728, 0.98299, 0.98868, 0.99435,
};

/* Avoid an int/float conversion. */
const float x264_log2_lz_lut[32] =
{
    31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
};

#define MAX_CHROMA_LAMBDA_OFFSET 36
static const uint16_t x264_chroma_lambda2_offset_tab[MAX_CHROMA_LAMBDA_OFFSET+1] =
{
       16,    20,    25,    32,    40,    50,
       64,    80,   101,   128,   161,   203,
      256,   322,   406,   512,   645,   812,
     1024,  1290,  1625,  2048,  2580,  3250,
     4096,  5160,  6501,  8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015,
    65535
};

//4-D_from L0_4x4 to D_L0_8x8,
static const uint8_t i_sub_mb_p_cost_table[4] =
{
    5, 3, 3, 1
};

static void x264_analyse_update_cache( x264_t *h, x264_mb_analysis_t *a );

static uint16_t x264_cost_ref[QP_MAX+1][3][33]; //3 - ref index?
static uint16_t x264_cost_i4x4_mode[(QP_MAX+2)*32];

float *x264_analyse_prepare_costs( x264_t *h )
{
    float *logs = x264_malloc( (2*4*2048+1)*sizeof(float) );
    if( !logs )
        return NULL;
    logs[0] = 0.718f;
    for( int i = 1; i <= 2*4*2048; i++ )
        logs[i] = log2f(i+1)*2 + 1.718f;
    return logs;
}

int x264_analyse_init_costs( x264_t *h, float *logs, int qp )
{
    int lambda = x264_lambda_tab[qp];
    if( h->cost_mv[qp] )
        return 0;
    /* factor of 4 from qpel, 2 from sign, and 2 because mv can be opposite from mvp */
    CHECKED_MALLOC( h->cost_mv[qp], (4*4*2048 + 1) * sizeof(uint16_t) );
    h->cost_mv[qp] += 2*4*2048;
    for( int i = 0; i <= 2*4*2048; i++ )
    {
        h->cost_mv[qp][-i] =
        h->cost_mv[qp][i]  = X264_MIN( lambda * logs[i] + .5f, (1<<16)-1 );
    }
    for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 33; j++ )
            x264_cost_ref[qp][i][j] = X264_MIN( i ? lambda * bs_size_te( i, j ) : 0, (1<<16)-1 );
    //me_method is X264_ME_ESA or X264_ME_TESA need four pixel cost
    if( h->param.analyse.i_me_method >= X264_ME_ESA && !h->cost_mv_fpel[qp][0] )
    {
        for( int j = 0; j < 4; j++ )
        {
            CHECKED_MALLOC( h->cost_mv_fpel[qp][j], (4*2048 + 1) * sizeof(uint16_t) );
            h->cost_mv_fpel[qp][j] += 2*2048;
            for( int i = -2*2048; i < 2*2048; i++ )
                h->cost_mv_fpel[qp][j][i] = h->cost_mv[qp][i*4+j];
        }
    }
    uint16_t *cost_i4x4_mode = (uint16_t*)ALIGN((intptr_t)x264_cost_i4x4_mode,64) + qp*32;
    for( int i = 0; i < 17; i++ )
        cost_i4x4_mode[i] = 3*lambda*(i!=8);
    return 0;
fail:
    return -1;
}

void x264_analyse_free_costs( x264_t *h )
{
    for( int i = 0; i < QP_MAX+1; i++ )
    {
        if( h->cost_mv[i] )
            x264_free( h->cost_mv[i] - 2*4*2048 );
        if( h->cost_mv_fpel[i][0] )
            for( int j = 0; j < 4; j++ )
                x264_free( h->cost_mv_fpel[i][j] - 2*2048 );
    }
}

/* initialize an array of lambda*nbits for all possible mvs */
static void x264_mb_analyse_load_costs( x264_t *h, x264_mb_analysis_t *a )
{
    a->p_cost_mv = h->cost_mv[a->i_qp];
    a->p_cost_ref = x264_cost_ref[a->i_qp][x264_clip3(h->sh.i_num_ref_idx_l0_active-1,0,2)];
}

static void x264_mb_analyse_init_qp( x264_t *h, x264_mb_analysis_t *a, int qp )
{
    int effective_chroma_qp = h->chroma_qp_table[SPEC_QP(qp)] + X264_MAX( qp - QP_MAX_SPEC, 0 );
    a->i_lambda = x264_lambda_tab[qp];
    a->i_lambda2 = x264_lambda2_tab[qp];

    h->mb.i_psy_rd_lambda = a->i_lambda;
    /* Adjusting chroma lambda based on QP offset hurts PSNR but improves visual quality. */
    int chroma_offset_idx = X264_MIN( qp-effective_chroma_qp+12, MAX_CHROMA_LAMBDA_OFFSET );
    h->mb.i_chroma_lambda2_offset = h->param.analyse.b_psy ? x264_chroma_lambda2_offset_tab[chroma_offset_idx] : 256;

    if( qp > QP_MAX_SPEC )
    {
        h->nr_offset = h->nr_offset_emergency[qp-QP_MAX_SPEC-1];
        h->nr_residual_sum = h->nr_residual_sum_buf[1];
        h->nr_count = h->nr_count_buf[1];
        h->mb.b_noise_reduction = 1;
        qp = QP_MAX_SPEC; /* Out-of-spec QPs are just used for calculating lambda values. */
    }
    else
    {
        h->nr_offset = h->nr_offset_denoise;
        h->nr_residual_sum = h->nr_residual_sum_buf[0];
        h->nr_count = h->nr_count_buf[0];
        h->mb.b_noise_reduction = 0;
    }

    a->i_qp = h->mb.i_qp = qp;
    h->mb.i_chroma_qp = h->chroma_qp_table[qp];
}

static void x264_mb_analyse_init( x264_t *h, x264_mb_analysis_t *a, int qp )
{
    int subme = h->param.analyse.i_subpel_refine;

    /* mbrd == 1 -> RD mode decision */
    /* mbrd == 2 -> RD refinement */
    /* mbrd == 3 -> QPRD */
    a->i_mbrd = (subme>=6) + (subme>=8) + (h->param.analyse.i_subpel_refine>=10);
    h->mb.b_deblock_rdo = h->param.analyse.i_subpel_refine >= 9 && h->sh.i_disable_deblocking_filter_idc != 1;
    a->b_early_terminate = h->param.analyse.i_subpel_refine < 11;

    x264_mb_analyse_init_qp( h, a, qp );

    /* I: Intra part */
    a->i_satd_i16x16 =
    a->i_satd_i4x4   =
    a->i_satd_chroma = COST_MAX;

    /* non-RD PCM decision is inaccurate (as is psy-rd), so don't do it.
     * PCM cost can overflow with high lambda2, so cap it at COST_MAX. */
    uint64_t pcm_cost = ((uint64_t)X264_PCM_COST*a->i_lambda2 + 128) >> 8;
    a->i_satd_pcm = !h->mb.i_psy_rd && a->i_mbrd && pcm_cost < COST_MAX ? pcm_cost : COST_MAX;

    a->b_fast_intra = 0;
    a->b_avoid_topright = 0;
    h->mb.i_skip_intra =
        h->mb.b_lossless ? 0 :
        a->i_mbrd ? 2 : 1;

    /* II: Inter part P/B frame */
    if( h->sh.i_type != SLICE_TYPE_I )
    {
        int i_fmv_range = 4 * h->param.analyse.i_mv_range;
        // limit motion search to a slightly smaller range than the theoretical limit,
        // since the search may go a few iterations past its given range
        int i_fpel_border = 6; // umh: 1 for diamond, 2 for octagon, 2 for hpel

        /* Calculate max allowed MV range */
#define CLIP_FMV(mv) x264_clip3( mv, -i_fmv_range, i_fmv_range-1 )
        h->mb.mv_min[0] = 4*( -16*h->mb.i_mb_x - 24 );
        h->mb.mv_max[0] = 4*( 16*( h->mb.i_mb_width - h->mb.i_mb_x - 1 ) + 24 );
        h->mb.mv_min_spel[0] = CLIP_FMV( h->mb.mv_min[0] );
        h->mb.mv_max_spel[0] = CLIP_FMV( h->mb.mv_max[0] );
        if( h->param.b_intra_refresh && h->sh.i_type == SLICE_TYPE_P )
        {
            int max_x = (h->fref[0]->i_pir_end_col * 16 - 3)*4; /* 3 pixels of hpel border */
            int max_mv = max_x - 4*16*h->mb.i_mb_x;
            /* If we're left of the refresh bar, don't reference right of it. */
            if( max_mv > 0 && h->mb.i_mb_x < h->fdec->i_pir_start_col )
                h->mb.mv_max_spel[0] = X264_MIN( h->mb.mv_max_spel[0], max_mv );
        }
        h->mb.mv_limit_fpel[0][0] = (h->mb.mv_min_spel[0]>>2) + i_fpel_border;
        h->mb.mv_limit_fpel[1][0] = (h->mb.mv_max_spel[0]>>2) - i_fpel_border;
        if( h->mb.i_mb_x == 0 )
        {
            int mb_y = h->mb.i_mb_y ;
            int thread_mvy_range = i_fmv_range;

            h->mb.mv_min[1] = 4*( -16*mb_y - 24 );
            h->mb.mv_max[1] = 4*( 16*( h->mb.i_mb_height - mb_y - 1 ) + 24 );
            h->mb.mv_min_spel[1] = x264_clip3( h->mb.mv_min[1], -i_fmv_range, i_fmv_range );
            h->mb.mv_max_spel[1] = CLIP_FMV( h->mb.mv_max[1] );
            h->mb.mv_max_spel[1] = X264_MIN( h->mb.mv_max_spel[1], thread_mvy_range*4 );
            h->mb.mv_limit_fpel[0][1] = (h->mb.mv_min_spel[1]>>2) + i_fpel_border;
            h->mb.mv_limit_fpel[1][1] = (h->mb.mv_max_spel[1]>>2) - i_fpel_border;
        }
#undef CLIP_FMV

        a->l0.me16x16.cost =
        a->l0.i_rd16x16    =
        a->l0.i_cost8x8    =
        a->l0.i_cost16x8   =
        a->l0.i_cost8x16   = COST_MAX;
        if( h->param.analyse.inter & X264_ANALYSE_PSUB8x8 )
            for( int i = 0; i < 4; i++ )
            {
                a->l0.i_cost4x4[i] =
                a->l0.i_cost8x4[i] =
                a->l0.i_cost4x8[i] = COST_MAX;
            }

        /* Fast intra decision */
        if( a->b_early_terminate && h->mb.i_mb_xy - h->sh.i_first_mb > 4 )
        {
            /* Always run in fast-intra mode for subme < 3 */
            if( h->mb.i_subpel_refine > 2 &&
              ( IS_INTRA( h->mb.i_mb_type_left[0] ) ||
                IS_INTRA( h->mb.i_mb_type_top ) ||
                IS_INTRA( h->mb.i_mb_type_topleft ) ||
                IS_INTRA( h->mb.i_mb_type_topright ) ||
                (h->sh.i_type == SLICE_TYPE_P && IS_INTRA( h->fref[0]->mb_type[h->mb.i_mb_xy] )) ||
                (h->mb.i_mb_xy - h->sh.i_first_mb < 3*(h->stat.frame.i_mb_count[I_4x4] + h->stat.frame.i_mb_count[I_8x8] + h->stat.frame.i_mb_count[I_16x16])) ) )
            { /* intra is likely */ }
            else
            {
                a->b_fast_intra = 1;
            }
        }
        h->mb.b_skip_mc = 0;
        //if intra refresh and mb_x in pir range, force intra
        if( h->param.b_intra_refresh && h->sh.i_type == SLICE_TYPE_P &&
            h->mb.i_mb_x >= h->fdec->i_pir_start_col && h->mb.i_mb_x <= h->fdec->i_pir_end_col )
        {
            a->b_force_intra = 1;
            a->b_fast_intra = 0;
            a->b_avoid_topright = h->mb.i_mb_x == h->fdec->i_pir_end_col;
        }
        else
            a->b_force_intra = 0;
    }
}

/* Prediction modes allowed for various combinations of neighbors. */
/* Terminated by a -1. */
/* In order, no neighbors, left, top, top/left, top/left/topleft */
static const int8_t i16x16_mode_available[5][5] =
{
    {I_PRED_16x16_DC_128, -1, -1, -1, -1},
    {I_PRED_16x16_DC_LEFT, I_PRED_16x16_H, -1, -1, -1},
    {I_PRED_16x16_DC_TOP, I_PRED_16x16_V, -1, -1, -1},
    {I_PRED_16x16_V, I_PRED_16x16_H, I_PRED_16x16_DC, -1, -1},
    {I_PRED_16x16_V, I_PRED_16x16_H, I_PRED_16x16_DC, I_PRED_16x16_P, -1},
};
//same as i16x16_mode_available
static const int8_t chroma_mode_available[5][5] =
{
    {I_PRED_CHROMA_DC_128, -1, -1, -1, -1},
    {I_PRED_CHROMA_DC_LEFT, I_PRED_CHROMA_H, -1, -1, -1},
    {I_PRED_CHROMA_DC_TOP, I_PRED_CHROMA_V, -1, -1, -1},
    {I_PRED_CHROMA_V, I_PRED_CHROMA_H, I_PRED_CHROMA_DC, -1, -1},
    {I_PRED_CHROMA_V, I_PRED_CHROMA_H, I_PRED_CHROMA_DC, I_PRED_CHROMA_P, -1},
};

//2:avoid_topright 0=top include topright 1 = top do not include topright
//5: In order, no neighbors, left, top, top/left, top/left/topleft.   top include topright?
static const int8_t i4x4_mode_available[2][5][10] =
{
    {
        {I_PRED_4x4_DC_128, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC_LEFT, I_PRED_4x4_H, I_PRED_4x4_HU, -1, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC_TOP, I_PRED_4x4_V, I_PRED_4x4_DDL, I_PRED_4x4_VL, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC, I_PRED_4x4_H, I_PRED_4x4_V, I_PRED_4x4_DDL, I_PRED_4x4_VL, I_PRED_4x4_HU, -1, -1, -1, -1},
        {I_PRED_4x4_DC, I_PRED_4x4_H, I_PRED_4x4_V, I_PRED_4x4_DDL, I_PRED_4x4_DDR, I_PRED_4x4_VR, I_PRED_4x4_HD, I_PRED_4x4_VL, I_PRED_4x4_HU, -1},
    },
    {
        {I_PRED_4x4_DC_128, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC_LEFT, I_PRED_4x4_H, I_PRED_4x4_HU, -1, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC_TOP, I_PRED_4x4_V, -1, -1, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC, I_PRED_4x4_H, I_PRED_4x4_V, I_PRED_4x4_HU, -1, -1, -1, -1, -1, -1},
        {I_PRED_4x4_DC, I_PRED_4x4_H, I_PRED_4x4_V, I_PRED_4x4_DDR, I_PRED_4x4_VR, I_PRED_4x4_HD, I_PRED_4x4_HU, -1, -1, -1},
    }
};

static  const int8_t *predict_16x16_mode_available( int i_neighbour )
{
    int idx = i_neighbour & (MB_TOP|MB_LEFT|MB_TOPLEFT);
    idx = (idx == (MB_TOP|MB_LEFT|MB_TOPLEFT)) ? 4 : idx & (MB_TOP|MB_LEFT);
    return i16x16_mode_available[idx];
}

static  const int8_t *predict_chroma_mode_available( int i_neighbour )
{
    int idx = i_neighbour & (MB_TOP|MB_LEFT|MB_TOPLEFT);
    idx = (idx == (MB_TOP|MB_LEFT|MB_TOPLEFT)) ? 4 : idx & (MB_TOP|MB_LEFT);
    return chroma_mode_available[idx];
}

static  const int8_t *predict_4x4_mode_available( int force_intra, int i_neighbour, int i )
{
    int avoid_topright = force_intra && ((i&5) == 5);
    int idx = i_neighbour & (MB_TOP|MB_LEFT|MB_TOPLEFT);
    idx = (idx == (MB_TOP|MB_LEFT|MB_TOPLEFT)) ? 4 : idx & (MB_TOP|MB_LEFT);
    return i4x4_mode_available[avoid_topright][idx];
}

/* Reset fenc satd scores cache for psy RD */
static  void x264_mb_init_fenc_cache( x264_t *h, int b_satd )
{
    if( !h->mb.i_psy_rd )
        return;
    /* Writes beyond the end of the array, but not a problem since fenc_satd_cache is right after. */
    h->mc.memzero_aligned( h->mb.pic.fenc_hadamard_cache, sizeof(h->mb.pic.fenc_hadamard_cache) );
    if( b_satd )
        h->mc.memzero_aligned( h->mb.pic.fenc_satd_cache, sizeof(h->mb.pic.fenc_satd_cache) );
}

static void x264_mb_analyse_intra_chroma( x264_t *h, x264_mb_analysis_t *a )
{
    if( a->i_satd_chroma < COST_MAX )
        return;
    //get the avalibe mode array by the neighbour
    const int8_t *predict_mode = predict_chroma_mode_available( h->mb.i_neighbour_intra );
    int chromapix = h->luma2chroma_pixel[PIXEL_16x16];

    /* Prediction selection for chroma */
    //if predict_mode[3] = I_PRED_CHROMA_P
    if( predict_mode[3] >= 0 && !h->mb.b_lossless )
    {
        int satdu[4], satdv[4];
        //caculate the first three mode cost
        h->pixf.intra_mbcmp_x3_chroma( h->mb.pic.p_fenc[1], h->mb.pic.p_fdec[1], satdu );
        h->pixf.intra_mbcmp_x3_chroma( h->mb.pic.p_fenc[2], h->mb.pic.p_fdec[2], satdv );
        //then caculate the I_PRED_CHROMA_P cost
        h->predict_chroma[I_PRED_CHROMA_P]( h->mb.pic.p_fdec[1] );
        h->predict_chroma[I_PRED_CHROMA_P]( h->mb.pic.p_fdec[2] );
        satdu[I_PRED_CHROMA_P] = h->pixf.mbcmp[chromapix]( h->mb.pic.p_fdec[1], FDEC_STRIDE, h->mb.pic.p_fenc[1], FENC_STRIDE );
        satdv[I_PRED_CHROMA_P] = h->pixf.mbcmp[chromapix]( h->mb.pic.p_fdec[2], FDEC_STRIDE, h->mb.pic.p_fenc[2], FENC_STRIDE );
        //pick the least mode
        for( ; *predict_mode >= 0; predict_mode++ )
        {
            int i_mode = *predict_mode;
            int i_satd = satdu[i_mode] + satdv[i_mode] + a->i_lambda * bs_size_ue( i_mode );

            a->i_satd_chroma_dir[i_mode] = i_satd;
            COPY2_IF_LT( a->i_satd_chroma, i_satd, a->i_predict8x8chroma, i_mode );
        }
    }
    else
    {
        for( ; *predict_mode >= 0; predict_mode++ )
        {
            int i_satd;
            int i_mode = *predict_mode;

            /* we do the prediction */
            if( h->mb.b_lossless )
                x264_predict_lossless_chroma( h, i_mode );
            else
            {
                h->predict_chroma[i_mode]( h->mb.pic.p_fdec[1] );
                h->predict_chroma[i_mode]( h->mb.pic.p_fdec[2] );
            }

            /* we calculate the cost */
            i_satd = h->pixf.mbcmp[chromapix]( h->mb.pic.p_fdec[1], FDEC_STRIDE, h->mb.pic.p_fenc[1], FENC_STRIDE ) +
                     h->pixf.mbcmp[chromapix]( h->mb.pic.p_fdec[2], FDEC_STRIDE, h->mb.pic.p_fenc[2], FENC_STRIDE ) +
                     a->i_lambda * bs_size_ue( x264_mb_chroma_pred_mode_fix[i_mode] );

            a->i_satd_chroma_dir[i_mode] = i_satd;
            COPY2_IF_LT( a->i_satd_chroma, i_satd, a->i_predict8x8chroma, i_mode );
        }
    }

    h->mb.i_chroma_pred_mode = a->i_predict8x8chroma;
}

/* FIXME: should we do any sort of merged chroma analysis with 4:4:4? */
static void x264_mb_analyse_intra( x264_t *h, x264_mb_analysis_t *a, int i_satd_inter )
{
    const unsigned int flags = h->sh.i_type == SLICE_TYPE_I ? h->param.analyse.intra : h->param.analyse.inter;
    pixel *p_src = h->mb.pic.p_fenc[0];
    pixel *p_dst = h->mb.pic.p_fdec[0];
    static const int8_t intra_analysis_shortcut[2][2][2][5] =
    {
        {{{I_PRED_4x4_HU, -1, -1, -1, -1},
          {I_PRED_4x4_DDL, I_PRED_4x4_VL, -1, -1, -1}},
         {{I_PRED_4x4_DDR, I_PRED_4x4_HD, I_PRED_4x4_HU, -1, -1},
          {I_PRED_4x4_DDL, I_PRED_4x4_DDR, I_PRED_4x4_VR, I_PRED_4x4_VL, -1}}},
        {{{I_PRED_4x4_HU, -1, -1, -1, -1},
          {-1, -1, -1, -1, -1}},
         {{I_PRED_4x4_DDR, I_PRED_4x4_HD, I_PRED_4x4_HU, -1, -1},
          {I_PRED_4x4_DDR, I_PRED_4x4_VR, -1, -1, -1}}},
    };

    int idx;
    int lambda = a->i_lambda;

    /*---------------- Try all mode and calculate their score ---------------*/
    {
        const int8_t *predict_mode = predict_16x16_mode_available( h->mb.i_neighbour_intra );

        /* Not heavily tuned */
        static const uint8_t i16x16_thresh_lut[11] = { 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4 };
        int i16x16_thresh = a->b_fast_intra ? (i16x16_thresh_lut[h->mb.i_subpel_refine]*i_satd_inter)>>1 : COST_MAX;

        //if predict_mode[3] = I_PRED_CHROMA_P
        if( !h->mb.b_lossless && predict_mode[3] >= 0 )
        {
            h->pixf.intra_mbcmp_x3_16x16( p_src, p_dst, a->i_satd_i16x16_dir );
            a->i_satd_i16x16_dir[I_PRED_16x16_V] += lambda * bs_size_ue(0);
            a->i_satd_i16x16_dir[I_PRED_16x16_H] += lambda * bs_size_ue(1);
            a->i_satd_i16x16_dir[I_PRED_16x16_DC] += lambda * bs_size_ue(2);
            COPY2_IF_LT( a->i_satd_i16x16, a->i_satd_i16x16_dir[I_PRED_16x16_V], a->i_predict16x16, I_PRED_16x16_V );
            COPY2_IF_LT( a->i_satd_i16x16, a->i_satd_i16x16_dir[I_PRED_16x16_H], a->i_predict16x16, I_PRED_16x16_H );
            COPY2_IF_LT( a->i_satd_i16x16, a->i_satd_i16x16_dir[I_PRED_16x16_DC], a->i_predict16x16, I_PRED_16x16_DC );

            /* Plane is expensive, so don't check it unless one of the previous modes was useful. */
            if( a->i_satd_i16x16 <= i16x16_thresh )
            {
                h->predict_16x16[I_PRED_16x16_P]( p_dst );
                a->i_satd_i16x16_dir[I_PRED_16x16_P] = h->pixf.mbcmp[PIXEL_16x16]( p_dst, FDEC_STRIDE, p_src, FENC_STRIDE );
                a->i_satd_i16x16_dir[I_PRED_16x16_P] += lambda * bs_size_ue(3);
                COPY2_IF_LT( a->i_satd_i16x16, a->i_satd_i16x16_dir[I_PRED_16x16_P], a->i_predict16x16, 3 );
            }
        }
        else
        {
            for( ; *predict_mode >= 0; predict_mode++ )
            {
                int i_satd;
                int i_mode = *predict_mode;

                if( h->mb.b_lossless )
                    x264_predict_lossless_16x16( h, 0, i_mode );
                else
                    h->predict_16x16[i_mode]( p_dst );

                i_satd = h->pixf.mbcmp[PIXEL_16x16]( p_dst, FDEC_STRIDE, p_src, FENC_STRIDE ) +
                         lambda * bs_size_ue( x264_mb_pred_mode16x16_fix[i_mode] );
                COPY2_IF_LT( a->i_satd_i16x16, i_satd, a->i_predict16x16, i_mode );
                a->i_satd_i16x16_dir[i_mode] = i_satd;
            }
        }

        if( a->i_satd_i16x16 > i16x16_thresh )
            return;
    }

    /* 4x4 prediction selection */
    if( flags & X264_ANALYSE_I4x4 )
    {
        int i_cost = lambda * (24+16); /* 24from JVT (SATD0), 16 from base predmode costs */
        int i_satd_thresh = a->b_early_terminate ? X264_MIN( i_satd_inter, a->i_satd_i16x16 ) : COST_MAX;
        h->mb.i_cbp_luma = 0;

        if( a->b_early_terminate && a->i_mbrd )
            i_satd_thresh = i_satd_thresh * (10-a->b_fast_intra)/8;

        for( idx = 0;; idx++ )
        {
            pixel *p_src_by = p_src + block_idx_xy_fenc[idx];
            pixel *p_dst_by = p_dst + block_idx_xy_fdec[idx];
            int i_best = COST_MAX;
            int i_pred_mode = x264_mb_predict_intra4x4_mode( h, idx );

            const int8_t *predict_mode = predict_4x4_mode_available( a->b_avoid_topright, h->mb.i_neighbour4[idx], idx );

            if( (h->mb.i_neighbour4[idx] & (MB_TOPRIGHT|MB_TOP)) == MB_TOP )
                /* emulate missing topright samples */
                MPIXEL_X4( &p_dst_by[4 - FDEC_STRIDE] ) = PIXEL_SPLAT_X4( p_dst_by[3 - FDEC_STRIDE] );

            {
                if( !h->mb.b_lossless && predict_mode[5] >= 0 )
                {
                    ALIGNED_ARRAY_16( int32_t, satd,[9] );
                    h->pixf.intra_mbcmp_x3_4x4( p_src_by, p_dst_by, satd );
                    int favor_vertical = satd[I_PRED_4x4_H] > satd[I_PRED_4x4_V];
                    satd[i_pred_mode] -= 3 * lambda;
                    i_best = satd[I_PRED_4x4_DC]; a->i_predict4x4[idx] = I_PRED_4x4_DC;
                    COPY2_IF_LT( i_best, satd[I_PRED_4x4_H], a->i_predict4x4[idx], I_PRED_4x4_H );
                    COPY2_IF_LT( i_best, satd[I_PRED_4x4_V], a->i_predict4x4[idx], I_PRED_4x4_V );

                    /* Take analysis shortcuts: don't analyse modes that are too
                     * far away direction-wise from the favored mode. */
                    if( a->i_mbrd < 1 + a->b_fast_intra )
                        predict_mode = intra_analysis_shortcut[a->b_avoid_topright][predict_mode[8] >= 0][favor_vertical];
                    else
                        predict_mode += 3;
                }

                if( i_best > 0 )
                {
                    for( ; *predict_mode >= 0; predict_mode++ )
                    {
                        int i_satd;
                        int i_mode = *predict_mode;

                        if( h->mb.b_lossless )
                            x264_predict_lossless_4x4( h, p_dst_by, 0, idx, i_mode );
                        else
                            h->predict_4x4[i_mode]( p_dst_by );

                        i_satd = h->pixf.mbcmp[PIXEL_4x4]( p_dst_by, FDEC_STRIDE, p_src_by, FENC_STRIDE );
                        if( i_pred_mode == x264_mb_pred_mode4x4_fix(i_mode) )
                        {
                            i_satd -= lambda * 3;
                            if( i_satd <= 0 )
                            {
                                i_best = i_satd;
                                a->i_predict4x4[idx] = i_mode;
                                break;
                            }
                        }

                        COPY2_IF_LT( i_best, i_satd, a->i_predict4x4[idx], i_mode );
                    }
                }

                i_cost += i_best + 3 * lambda;
                if( i_cost > i_satd_thresh || idx == 15 )
                    break;
                if( h->mb.b_lossless )
                    x264_predict_lossless_4x4( h, p_dst_by, 0, idx, a->i_predict4x4[idx] );
                else
                    h->predict_4x4[a->i_predict4x4[idx]]( p_dst_by );
                h->mb.cache.intra4x4_pred_mode[x264_scan8[idx]] = a->i_predict4x4[idx];
            }
            /* we need to encode this block now (for next ones) */
            x264_mb_encode_i4x4( h, 0, idx, a->i_qp, a->i_predict4x4[idx], 0 );
        }
        if( idx == 15 )
        {
            a->i_satd_i4x4 = i_cost;
            if( h->mb.i_skip_intra )
            {
                h->mc.copy[PIXEL_16x16]( h->mb.pic.i4x4_fdec_buf, 16, p_dst, FDEC_STRIDE, 16 );
                h->mb.pic.i4x4_nnz_buf[0] = M32( &h->mb.cache.non_zero_count[x264_scan8[ 0]] );
                h->mb.pic.i4x4_nnz_buf[1] = M32( &h->mb.cache.non_zero_count[x264_scan8[ 2]] );
                h->mb.pic.i4x4_nnz_buf[2] = M32( &h->mb.cache.non_zero_count[x264_scan8[ 8]] );
                h->mb.pic.i4x4_nnz_buf[3] = M32( &h->mb.cache.non_zero_count[x264_scan8[10]] );
                h->mb.pic.i4x4_cbp = h->mb.i_cbp_luma;
                if( h->mb.i_skip_intra == 2 )
                    h->mc.memcpy_aligned( h->mb.pic.i4x4_dct_buf, h->dct.luma4x4, sizeof(h->mb.pic.i4x4_dct_buf) );
            }
        }
        else
            a->i_satd_i4x4 = COST_MAX;
    }
}

static void x264_intra_rd( x264_t *h, x264_mb_analysis_t *a, int i_satd_thresh )
{
    if( !a->b_early_terminate )
        i_satd_thresh = COST_MAX;

    if( a->i_satd_i16x16 < i_satd_thresh )
    {
        h->mb.i_type = I_16x16;
        x264_analyse_update_cache( h, a );
        a->i_satd_i16x16 = x264_rd_cost_mb( h, a->i_lambda2 );
    }
    else
        a->i_satd_i16x16 = COST_MAX;

    if( a->i_satd_i4x4 < i_satd_thresh )
    {
        h->mb.i_type = I_4x4;
        x264_analyse_update_cache( h, a );
        a->i_satd_i4x4 = x264_rd_cost_mb( h, a->i_lambda2 );
    }
    else
        a->i_satd_i4x4 = COST_MAX;
}

static void x264_intra_rd_refine( x264_t *h, x264_mb_analysis_t *a )
{
    uint64_t i_satd, i_best;
    int plane_count = 1;
    h->mb.i_skip_intra = 0;

    if( h->mb.i_type == I_16x16 )
    {
        int old_pred_mode = a->i_predict16x16;
        const int8_t *predict_mode = predict_16x16_mode_available( h->mb.i_neighbour_intra );
        int i_thresh = a->b_early_terminate ? a->i_satd_i16x16_dir[old_pred_mode] * 9/8 : COST_MAX;
        i_best = a->i_satd_i16x16;
        for( ; *predict_mode >= 0; predict_mode++ )
        {
            int i_mode = *predict_mode;
            if( i_mode == old_pred_mode || a->i_satd_i16x16_dir[i_mode] > i_thresh )
                continue;
            h->mb.i_intra16x16_pred_mode = i_mode;
            i_satd = x264_rd_cost_mb( h, a->i_lambda2 );
            COPY2_IF_LT( i_best, i_satd, a->i_predict16x16, i_mode );
        }
    }

    /* RD selection for chroma prediction */
    if( 1 )
    {
        const int8_t *predict_mode = predict_chroma_mode_available( h->mb.i_neighbour_intra );
        if( predict_mode[1] >= 0 )
        {
            int8_t predict_mode_sorted[4];
            int i_max;
            int i_thresh = a->b_early_terminate ? a->i_satd_chroma * 5/4 : COST_MAX;

            for( i_max = 0; *predict_mode >= 0; predict_mode++ )
            {
                int i_mode = *predict_mode;
                if( a->i_satd_chroma_dir[i_mode] < i_thresh && i_mode != a->i_predict8x8chroma )
                    predict_mode_sorted[i_max++] = i_mode;
            }

            if( i_max > 0 )
            {
                int i_cbp_chroma_best = h->mb.i_cbp_chroma;
                int i_chroma_lambda = x264_lambda2_tab[h->mb.i_chroma_qp];
                /* the previous thing encoded was x264_intra_rd(), so the pixels and
                 * coefs for the current chroma mode are still around, so we only
                 * have to recount the bits. */
                i_best = x264_rd_cost_chroma( h, i_chroma_lambda, a->i_predict8x8chroma, 0 );
                for( int i = 0; i < i_max; i++ )
                {
                    int i_mode = predict_mode_sorted[i];
                    if( h->mb.b_lossless )
                        x264_predict_lossless_chroma( h, i_mode );
                    else
                    {
                        h->predict_chroma[i_mode]( h->mb.pic.p_fdec[1] );
                        h->predict_chroma[i_mode]( h->mb.pic.p_fdec[2] );
                    }
                    /* if we've already found a mode that needs no residual, then
                     * probably any mode with a residual will be worse.
                     * so avoid dct on the remaining modes to improve speed. */
                    i_satd = x264_rd_cost_chroma( h, i_chroma_lambda, i_mode, h->mb.i_cbp_chroma != 0x00 );
                    COPY3_IF_LT( i_best, i_satd, a->i_predict8x8chroma, i_mode, i_cbp_chroma_best, h->mb.i_cbp_chroma );
                }
                h->mb.i_chroma_pred_mode = a->i_predict8x8chroma;
                h->mb.i_cbp_chroma = i_cbp_chroma_best;
            }
        }
    }

    if( h->mb.i_type == I_4x4 )
    {
        pixel4 pels[3][4] = {{0}}; // doesn't need initting, just shuts up a gcc warning
        int nnz[3] = {0};
        for( int idx = 0; idx < 16; idx++ )
        {
            pixel *dst[3] = {h->mb.pic.p_fdec[0] + block_idx_xy_fdec[idx],
                             h->mb.pic.p_fdec[1] + block_idx_xy_fdec[idx],
                             h->mb.pic.p_fdec[2] + block_idx_xy_fdec[idx]};
            i_best = COST_MAX64;

            const int8_t *predict_mode = predict_4x4_mode_available( a->b_avoid_topright, h->mb.i_neighbour4[idx], idx );

            if( (h->mb.i_neighbour4[idx] & (MB_TOPRIGHT|MB_TOP)) == MB_TOP )
                for( int p = 0; p < plane_count; p++ )
                    /* emulate missing topright samples */
                    MPIXEL_X4( dst[p]+4-FDEC_STRIDE ) = PIXEL_SPLAT_X4( dst[p][3-FDEC_STRIDE] );

            for( ; *predict_mode >= 0; predict_mode++ )
            {
                int i_mode = *predict_mode;
                i_satd = x264_rd_cost_i4x4( h, a->i_lambda2, idx, i_mode );

                if( i_best > i_satd )
                {
                    a->i_predict4x4[idx] = i_mode;
                    i_best = i_satd;
                    for( int p = 0; p < plane_count; p++ )
                    {
                        pels[p][0] = MPIXEL_X4( dst[p]+0*FDEC_STRIDE );
                        pels[p][1] = MPIXEL_X4( dst[p]+1*FDEC_STRIDE );
                        pels[p][2] = MPIXEL_X4( dst[p]+2*FDEC_STRIDE );
                        pels[p][3] = MPIXEL_X4( dst[p]+3*FDEC_STRIDE );
                        nnz[p] = h->mb.cache.non_zero_count[x264_scan8[idx+p*16]];
                    }
                }
            }

            for( int p = 0; p < plane_count; p++ )
            {
                MPIXEL_X4( dst[p]+0*FDEC_STRIDE ) = pels[p][0];
                MPIXEL_X4( dst[p]+1*FDEC_STRIDE ) = pels[p][1];
                MPIXEL_X4( dst[p]+2*FDEC_STRIDE ) = pels[p][2];
                MPIXEL_X4( dst[p]+3*FDEC_STRIDE ) = pels[p][3];
                h->mb.cache.non_zero_count[x264_scan8[idx+p*16]] = nnz[p];
            }

            h->mb.cache.intra4x4_pred_mode[x264_scan8[idx]] = a->i_predict4x4[idx];
        }
    }
}

#define LOAD_FENC(m, src, xoff, yoff) \
{ \
    (m)->p_cost_mv = a->p_cost_mv; \
    (m)->i_stride[0] = h->mb.pic.i_stride[0]; \
    (m)->i_stride[1] = h->mb.pic.i_stride[1]; \
    (m)->i_stride[2] = h->mb.pic.i_stride[2]; \
    (m)->p_fenc[0] = &(src)[0][(xoff)+(yoff)*FENC_STRIDE]; \
    (m)->p_fenc[1] = &(src)[1][((xoff)>>CHROMA_H_SHIFT)+((yoff)>>CHROMA_V_SHIFT)*FENC_STRIDE]; \
    (m)->p_fenc[2] = &(src)[2][((xoff)>>CHROMA_H_SHIFT)+((yoff)>>CHROMA_V_SHIFT)*FENC_STRIDE]; \
}

#define LOAD_HPELS(m, src,  ref, xoff, yoff) \
{ \
    (m)->p_fref[1] = &(src)[1][(xoff)+(yoff)*(m)->i_stride[0]]; \
    (m)->p_fref[2] = &(src)[2][(xoff)+(yoff)*(m)->i_stride[0]]; \
    (m)->p_fref[3] = &(src)[3][(xoff)+(yoff)*(m)->i_stride[0]]; \
    (m)->p_fref[4] = &(src)[4][(xoff)+((yoff)>>CHROMA_V_SHIFT)*(m)->i_stride[1]]; \
    (m)->integral = &h->mb.pic.p_integral[ref][(xoff)+(yoff)*(m)->i_stride[0]]; \
    (m)->i_ref = ref; \
}

#define REF_COST(ref) (a->p_cost_ref[ref])

static void x264_mb_analyse_inter_p16x16( x264_t *h, x264_mb_analysis_t *a )
{
    x264_me_t m;
    int i_mvc;
    ALIGNED_4( int16_t mvc[8][2] );
    int i_halfpel_thresh = INT_MAX;
    int *p_halfpel_thresh = (a->b_early_terminate && h->mb.pic.i_fref>1) ? &i_halfpel_thresh : NULL;

    /* 16x16 Search on all ref frame */
    m.i_pixel = PIXEL_16x16;
    LOAD_FENC( &m, h->mb.pic.p_fenc, 0, 0 );

    a->l0.me16x16.cost = INT_MAX;
    for( int i_ref = 0; i_ref < h->mb.pic.i_fref; i_ref++ )
    {
        m.i_ref_cost = REF_COST(i_ref );
        i_halfpel_thresh -= m.i_ref_cost;

        /* search with ref */
        LOAD_HPELS( &m, h->mb.pic.p_fref[i_ref], i_ref, 0, 0 );

        x264_mb_predict_mv_16x16( h, i_ref, m.mvp );

        if( h->mb.ref_blind_dupe == i_ref )
        {
            CP32( m.mv, a->l0.mvc[0][0] );
            x264_me_refine_qpel_refdupe( h, &m, p_halfpel_thresh );
        }
        else
        {
            x264_mb_predict_mv_ref16x16( h, i_ref, mvc, &i_mvc );
            x264_me_search_ref( h, &m, mvc, i_mvc, p_halfpel_thresh );
        }

        /* save mv for predicting neighbors */
        CP32( h->mb.mvr[i_ref][h->mb.i_mb_xy], m.mv );
        CP32( a->l0.mvc[i_ref][0], m.mv );

        /* early termination
         * SSD threshold would probably be better than SATD */
        if( i_ref == 0
            && a->b_try_skip
            && m.cost-m.cost_mv < 300*a->i_lambda
            &&  abs(m.mv[0]-h->mb.cache.pskip_mv[0])
              + abs(m.mv[1]-h->mb.cache.pskip_mv[1]) <= 1
            && x264_macroblock_probe_pskip( h ) )
        {
            h->mb.i_type = P_SKIP;
            x264_analyse_update_cache( h, a );
            return;
        }

        m.cost += m.i_ref_cost;
        i_halfpel_thresh += m.i_ref_cost;

        if( m.cost < a->l0.me16x16.cost )
            h->mc.memcpy_aligned( &a->l0.me16x16, &m, sizeof(x264_me_t) );
    }

    x264_macroblock_cache_ref( h, 0, 0, 4, 4, a->l0.me16x16.i_ref );

    h->mb.i_type = P_L0;
    if( a->i_mbrd )
    {
        x264_mb_init_fenc_cache( h, a->i_mbrd >= 2 || h->param.analyse.inter & X264_ANALYSE_PSUB8x8 );
        if( a->l0.me16x16.i_ref == 0 && M32( a->l0.me16x16.mv ) == M32( h->mb.cache.pskip_mv ) && !a->b_force_intra )
        {
            h->mb.i_partition = D_16x16;
            x264_macroblock_cache_mv_ptr( h, 0, 0, 4, 4, a->l0.me16x16.mv );
            a->l0.i_rd16x16 = x264_rd_cost_mb( h, a->i_lambda2 );
            if( !(h->mb.i_cbp_luma|h->mb.i_cbp_chroma) )
                h->mb.i_type = P_SKIP;
        }
    }
}

static void x264_mb_analyse_inter_p8x8( x264_t *h, x264_mb_analysis_t *a )
{
    /* Duplicate refs are rarely useful in p8x8 due to the high cost of the
     * reference frame flags.  Thus, if we're not doing mixedrefs, just
     * don't bother analysing the dupes. */
    const int i_ref = h->mb.ref_blind_dupe == a->l0.me16x16.i_ref ? 0 : a->l0.me16x16.i_ref;
    const int i_ref_cost = i_ref ? REF_COST( i_ref ) : 0;
    pixel **p_fenc = h->mb.pic.p_fenc;
    int i_mvc;
    int16_t (*mvc)[2] = a->l0.mvc[i_ref];

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_8x8;

    i_mvc = 1;
    CP32( mvc[0], a->l0.me16x16.mv );

    for( int i = 0; i < 4; i++ )
    {
        x264_me_t *m = &a->l0.me8x8[i];
        int x8 = i&1;
        int y8 = i>>1;

        m->i_pixel = PIXEL_8x8;
        m->i_ref_cost = i_ref_cost;

        LOAD_FENC( m, p_fenc, 8*x8, 8*y8 );
        LOAD_HPELS( m, h->mb.pic.p_fref[i_ref], i_ref, 8*x8, 8*y8 );

        x264_mb_predict_mv( h, 4*i, 2, m->mvp );
        x264_me_search( h, m, mvc, i_mvc );

        x264_macroblock_cache_mv_ptr( h, 2*x8, 2*y8, 2, 2, m->mv );

        CP32( mvc[i_mvc], m->mv );
        i_mvc++;

        a->i_satd8x8[i] = m->cost - m->cost_mv;

        /* mb type cost */
        m->cost += i_ref_cost;
        if( 1 || (h->param.analyse.inter & X264_ANALYSE_PSUB8x8) )
            m->cost += a->i_lambda * i_sub_mb_p_cost_table[D_L0_8x8];
    }

    a->l0.i_cost8x8 = a->l0.me8x8[0].cost + a->l0.me8x8[1].cost +
                      a->l0.me8x8[2].cost + a->l0.me8x8[3].cost;

    h->mb.i_sub_partition[0] = h->mb.i_sub_partition[1] =
    h->mb.i_sub_partition[2] = h->mb.i_sub_partition[3] = D_L0_8x8;
}

static void x264_mb_analyse_inter_p16x8( x264_t *h, x264_mb_analysis_t *a, int i_best_satd )
{
    x264_me_t m;
    pixel **p_fenc = h->mb.pic.p_fenc;
    ALIGNED_4( int16_t mvc[3][2] );

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_16x8;

    for( int i = 0; i < 2; i++ )
    {
        x264_me_t *l0m = &a->l0.me16x8[i];
        const int minref = X264_MIN( a->l0.me8x8[2*i].i_ref, a->l0.me8x8[2*i+1].i_ref );
        const int maxref = X264_MAX( a->l0.me8x8[2*i].i_ref, a->l0.me8x8[2*i+1].i_ref );
        const int ref8[2] = { minref, maxref };
        const int i_ref8s = ( ref8[0] == ref8[1] ) ? 1 : 2;

        m.i_pixel = PIXEL_16x8;

        LOAD_FENC( &m, p_fenc, 0, 8*i );
        l0m->cost = INT_MAX;
        for( int j = 0; j < i_ref8s; j++ )
        {
            const int i_ref = ref8[j];
            m.i_ref_cost = REF_COST( i_ref );

            /* if we skipped the 16x16 predictor, we wouldn't have to copy anything... */
            CP32( mvc[0], a->l0.mvc[i_ref][0] );
            CP32( mvc[1], a->l0.mvc[i_ref][2*i+1] );
            CP32( mvc[2], a->l0.mvc[i_ref][2*i+2] );

            LOAD_HPELS( &m, h->mb.pic.p_fref[i_ref], i_ref, 0, 8*i );

            x264_macroblock_cache_ref( h, 0, 2*i, 4, 2, i_ref );
            x264_mb_predict_mv( h, 8*i, 4, m.mvp );
            /* We can only take this shortcut if the first search was performed on ref0. */
            if( h->mb.ref_blind_dupe == i_ref && !ref8[0] )
            {
                /* We can just leave the MV from the previous ref search. */
                x264_me_refine_qpel_refdupe( h, &m, NULL );
            }
            else
                x264_me_search( h, &m, mvc, 3 );

            m.cost += m.i_ref_cost;

            if( m.cost < l0m->cost )
                h->mc.memcpy_aligned( l0m, &m, sizeof(x264_me_t) );
        }

        /* Early termination based on the current SATD score of partition[0]
           plus the estimated SATD score of partition[1] */
        if( a->b_early_terminate && (!i && l0m->cost + a->i_cost_est16x8[1] > i_best_satd * (4 + !!a->i_mbrd) / 4) )
        {
            a->l0.i_cost16x8 = COST_MAX;
            return;
        }

        x264_macroblock_cache_mv_ptr( h, 0, 2*i, 4, 2, l0m->mv );
        x264_macroblock_cache_ref( h, 0, 2*i, 4, 2, l0m->i_ref );
    }

    a->l0.i_cost16x8 = a->l0.me16x8[0].cost + a->l0.me16x8[1].cost;
}

static void x264_mb_analyse_inter_p8x16( x264_t *h, x264_mb_analysis_t *a, int i_best_satd )
{
    x264_me_t m;
    pixel **p_fenc = h->mb.pic.p_fenc;
    ALIGNED_4( int16_t mvc[3][2] );

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_8x16;

    for( int i = 0; i < 2; i++ )
    {
        x264_me_t *l0m = &a->l0.me8x16[i];
        const int minref = X264_MIN( a->l0.me8x8[i].i_ref, a->l0.me8x8[i+2].i_ref );
        const int maxref = X264_MAX( a->l0.me8x8[i].i_ref, a->l0.me8x8[i+2].i_ref );
        const int ref8[2] = { minref, maxref };
        const int i_ref8s = ( ref8[0] == ref8[1] ) ? 1 : 2;

        m.i_pixel = PIXEL_8x16;

        LOAD_FENC( &m, p_fenc, 8*i, 0 );
        l0m->cost = INT_MAX;
        for( int j = 0; j < i_ref8s; j++ )
        {
            const int i_ref = ref8[j];
            m.i_ref_cost = REF_COST( i_ref );

            CP32( mvc[0], a->l0.mvc[i_ref][0] );
            CP32( mvc[1], a->l0.mvc[i_ref][i+1] );
            CP32( mvc[2], a->l0.mvc[i_ref][i+3] );

            LOAD_HPELS( &m, h->mb.pic.p_fref[i_ref], i_ref, 8*i, 0 );

            x264_macroblock_cache_ref( h, 2*i, 0, 2, 4, i_ref );
            x264_mb_predict_mv( h, 4*i, 2, m.mvp );
            /* We can only take this shortcut if the first search was performed on ref0. */
            if( h->mb.ref_blind_dupe == i_ref && !ref8[0] )
            {
                /* We can just leave the MV from the previous ref search. */
                x264_me_refine_qpel_refdupe( h, &m, NULL );
            }
            else
                x264_me_search( h, &m, mvc, 3 );

            m.cost += m.i_ref_cost;

            if( m.cost < l0m->cost )
                h->mc.memcpy_aligned( l0m, &m, sizeof(x264_me_t) );
        }

        /* Early termination based on the current SATD score of partition[0]
           plus the estimated SATD score of partition[1] */
        if( a->b_early_terminate && (!i && l0m->cost + a->i_cost_est8x16[1] > i_best_satd * (4 + !!a->i_mbrd) / 4) )
        {
            a->l0.i_cost8x16 = COST_MAX;
            return;
        }

        x264_macroblock_cache_mv_ptr( h, 2*i, 0, 2, 4, l0m->mv );
        x264_macroblock_cache_ref( h, 2*i, 0, 2, 4, l0m->i_ref );
    }

    a->l0.i_cost8x16 = a->l0.me8x16[0].cost + a->l0.me8x16[1].cost;
}


static  int x264_mb_analyse_inter_p4x4_chroma_internal( x264_t *h, x264_mb_analysis_t *a,
                                                                     pixel **p_fref, int i8x8, int size)
{
    ALIGNED_ARRAY_N( pixel, pix1,[16*16] );
    pixel *pix2 = pix1+8;
    int i_stride = h->mb.pic.i_stride[1];
    int chroma_h_shift = 1;
    int chroma_v_shift = 1;
    int or = 8*(i8x8&1) + (4>>chroma_v_shift)*(i8x8&2)*i_stride;

#define CHROMA4x4MC( width, height, me, x, y ) \
    do{ \
        int offset = x + (2>>chroma_v_shift)*16*y; \
        int chroma_height = (2>>chroma_v_shift)*height; \
        h->mc.mc_chroma( &pix1[offset], &pix2[offset], 16, &p_fref[4][or+2*x+(2>>chroma_v_shift)*y*i_stride], i_stride, \
                         (me).mv[0], (2>>chroma_v_shift)*((me).mv[1]+0), width, chroma_height ); \
    }while(0)

    if( size == PIXEL_4x4 )
    {
        x264_me_t *m = a->l0.me4x4[i8x8];
        CHROMA4x4MC( 2,2, m[0], 0,0 );
        CHROMA4x4MC( 2,2, m[1], 2,0 );
        CHROMA4x4MC( 2,2, m[2], 0,2 );
        CHROMA4x4MC( 2,2, m[3], 2,2 );
    }
    else if( size == PIXEL_8x4 )
    {
        x264_me_t *m = a->l0.me8x4[i8x8];
        CHROMA4x4MC( 4,2, m[0], 0,0 );
        CHROMA4x4MC( 4,2, m[1], 0,2 );
    }
    else
    {
        x264_me_t *m = a->l0.me4x8[i8x8];
        CHROMA4x4MC( 2,4, m[0], 0,0 );
        CHROMA4x4MC( 2,4, m[1], 2,0 );
    }
#undef CHROMA4x4MC

    int oe = (8>>chroma_h_shift)*(i8x8&1) + (4>>chroma_v_shift)*(i8x8&2)*FENC_STRIDE;
    int chromapix = PIXEL_4x4;
    return h->pixf.mbcmp[chromapix]( &h->mb.pic.p_fenc[1][oe], FENC_STRIDE, pix1, 16 )
         + h->pixf.mbcmp[chromapix]( &h->mb.pic.p_fenc[2][oe], FENC_STRIDE, pix2, 16 );
}

static int x264_mb_analyse_inter_p4x4_chroma( x264_t *h, x264_mb_analysis_t *a, pixel **p_fref, int i8x8, int size )
{
    return x264_mb_analyse_inter_p4x4_chroma_internal( h, a, p_fref, i8x8, size );
}

static void x264_mb_analyse_inter_p4x4( x264_t *h, x264_mb_analysis_t *a, int i8x8 )
{
    pixel **p_fref = h->mb.pic.p_fref[a->l0.me8x8[i8x8].i_ref];
    pixel **p_fenc = h->mb.pic.p_fenc;
    const int i_ref = a->l0.me8x8[i8x8].i_ref;

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_8x8;

    for( int i4x4 = 0; i4x4 < 4; i4x4++ )
    {
        const int idx = 4*i8x8 + i4x4;
        const int x4 = block_idx_x[idx];
        const int y4 = block_idx_y[idx];
        const int i_mvc = (i4x4 == 0);

        x264_me_t *m = &a->l0.me4x4[i8x8][i4x4];

        m->i_pixel = PIXEL_4x4;

        LOAD_FENC( m, p_fenc, 4*x4, 4*y4 );
        LOAD_HPELS( m, p_fref, i_ref, 4*x4, 4*y4 );

        x264_mb_predict_mv( h, idx, 1, m->mvp );
        x264_me_search( h, m, &a->l0.me8x8[i8x8].mv, i_mvc );

        x264_macroblock_cache_mv_ptr( h, x4, y4, 1, 1, m->mv );
    }
    a->l0.i_cost4x4[i8x8] = a->l0.me4x4[i8x8][0].cost +
                            a->l0.me4x4[i8x8][1].cost +
                            a->l0.me4x4[i8x8][2].cost +
                            a->l0.me4x4[i8x8][3].cost +
                            REF_COST( i_ref ) +
                            a->i_lambda * i_sub_mb_p_cost_table[D_L0_4x4];
    if( h->mb.b_chroma_me )
        a->l0.i_cost4x4[i8x8] += x264_mb_analyse_inter_p4x4_chroma( h, a, p_fref, i8x8, PIXEL_4x4 );
}

static void x264_mb_analyse_inter_p8x4( x264_t *h, x264_mb_analysis_t *a, int i8x8 )
{
    pixel **p_fref = h->mb.pic.p_fref[a->l0.me8x8[i8x8].i_ref];
    pixel **p_fenc = h->mb.pic.p_fenc;
    const int i_ref = a->l0.me8x8[i8x8].i_ref;

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_8x8;

    for( int i8x4 = 0; i8x4 < 2; i8x4++ )
    {
        const int idx = 4*i8x8 + 2*i8x4;
        const int x4 = block_idx_x[idx];
        const int y4 = block_idx_y[idx];
        const int i_mvc = (i8x4 == 0);

        x264_me_t *m = &a->l0.me8x4[i8x8][i8x4];

        m->i_pixel = PIXEL_8x4;

        LOAD_FENC( m, p_fenc, 4*x4, 4*y4 );
        LOAD_HPELS( m, p_fref, i_ref, 4*x4, 4*y4 );

        x264_mb_predict_mv( h, idx, 2, m->mvp );
        x264_me_search( h, m, &a->l0.me4x4[i8x8][0].mv, i_mvc );

        x264_macroblock_cache_mv_ptr( h, x4, y4, 2, 1, m->mv );
    }
    a->l0.i_cost8x4[i8x8] = a->l0.me8x4[i8x8][0].cost + a->l0.me8x4[i8x8][1].cost +
                            REF_COST( i_ref ) +
                            a->i_lambda * i_sub_mb_p_cost_table[D_L0_8x4];
    if( h->mb.b_chroma_me )
        a->l0.i_cost8x4[i8x8] += x264_mb_analyse_inter_p4x4_chroma( h, a, p_fref, i8x8, PIXEL_8x4 );
}

static void x264_mb_analyse_inter_p4x8( x264_t *h, x264_mb_analysis_t *a, int i8x8 )
{
    pixel **p_fref = h->mb.pic.p_fref[a->l0.me8x8[i8x8].i_ref];
    pixel **p_fenc = h->mb.pic.p_fenc;
    const int i_ref = a->l0.me8x8[i8x8].i_ref;

    /* XXX Needed for x264_mb_predict_mv */
    h->mb.i_partition = D_8x8;

    for( int i4x8 = 0; i4x8 < 2; i4x8++ )
    {
        const int idx = 4*i8x8 + i4x8;
        const int x4 = block_idx_x[idx];
        const int y4 = block_idx_y[idx];
        const int i_mvc = (i4x8 == 0);

        x264_me_t *m = &a->l0.me4x8[i8x8][i4x8];

        m->i_pixel = PIXEL_4x8;

        LOAD_FENC( m, p_fenc, 4*x4, 4*y4 );
        LOAD_HPELS( m, p_fref, i_ref, 4*x4, 4*y4 );

        x264_mb_predict_mv( h, idx, 1, m->mvp );
        x264_me_search( h, m, &a->l0.me4x4[i8x8][0].mv, i_mvc );

        x264_macroblock_cache_mv_ptr( h, x4, y4, 1, 2, m->mv );
    }
    a->l0.i_cost4x8[i8x8] = a->l0.me4x8[i8x8][0].cost + a->l0.me4x8[i8x8][1].cost +
                            REF_COST( i_ref ) +
                            a->i_lambda * i_sub_mb_p_cost_table[D_L0_4x8];
    if( h->mb.b_chroma_me )
        a->l0.i_cost4x8[i8x8] += x264_mb_analyse_inter_p4x4_chroma( h, a, p_fref, i8x8, PIXEL_4x8 );
}

static  void x264_mb_cache_mv_p8x8( x264_t *h, x264_mb_analysis_t *a, int i )
{
    int x = 2*(i&1);
    int y = i&2;

    switch( h->mb.i_sub_partition[i] )
    {
        case D_L0_8x8:
            x264_macroblock_cache_mv_ptr( h, x, y, 2, 2, a->l0.me8x8[i].mv );
            break;
        case D_L0_8x4:
            x264_macroblock_cache_mv_ptr( h, x, y+0, 2, 1, a->l0.me8x4[i][0].mv );
            x264_macroblock_cache_mv_ptr( h, x, y+1, 2, 1, a->l0.me8x4[i][1].mv );
            break;
        case D_L0_4x8:
            x264_macroblock_cache_mv_ptr( h, x+0, y, 1, 2, a->l0.me4x8[i][0].mv );
            x264_macroblock_cache_mv_ptr( h, x+1, y, 1, 2, a->l0.me4x8[i][1].mv );
            break;
        case D_L0_4x4:
            x264_macroblock_cache_mv_ptr( h, x+0, y+0, 1, 1, a->l0.me4x4[i][0].mv );
            x264_macroblock_cache_mv_ptr( h, x+1, y+0, 1, 1, a->l0.me4x4[i][1].mv );
            x264_macroblock_cache_mv_ptr( h, x+0, y+1, 1, 1, a->l0.me4x4[i][2].mv );
            x264_macroblock_cache_mv_ptr( h, x+1, y+1, 1, 1, a->l0.me4x4[i][3].mv );
            break;
        default:
            x264_log( h, X264_LOG_ERROR, "internal error\n" );
            break;
    }
}

static void x264_mb_analyse_p_rd( x264_t *h, x264_mb_analysis_t *a, int i_satd )
{
    int thresh = a->b_early_terminate ? i_satd * 5/4 + 1 : COST_MAX;

    h->mb.i_type = P_L0;
    if( a->l0.i_rd16x16 == COST_MAX && (!a->b_early_terminate || a->l0.me16x16.cost <= i_satd * 3/2) )
    {
        h->mb.i_partition = D_16x16;
        x264_analyse_update_cache( h, a );
        a->l0.i_rd16x16 = x264_rd_cost_mb( h, a->i_lambda2 );
    }

    if( a->l0.i_cost16x8 < thresh )
    {
        h->mb.i_partition = D_16x8;
        x264_analyse_update_cache( h, a );
        a->l0.i_cost16x8 = x264_rd_cost_mb( h, a->i_lambda2 );
    }
    else
        a->l0.i_cost16x8 = COST_MAX;

    if( a->l0.i_cost8x16 < thresh )
    {
        h->mb.i_partition = D_8x16;
        x264_analyse_update_cache( h, a );
        a->l0.i_cost8x16 = x264_rd_cost_mb( h, a->i_lambda2 );
    }
    else
        a->l0.i_cost8x16 = COST_MAX;

    if( a->l0.i_cost8x8 < thresh )
    {
        h->mb.i_type = P_8x8;
        h->mb.i_partition = D_8x8;
        if( h->param.analyse.inter & X264_ANALYSE_PSUB8x8 )
        {
            x264_macroblock_cache_ref( h, 0, 0, 2, 2, a->l0.me8x8[0].i_ref );
            x264_macroblock_cache_ref( h, 2, 0, 2, 2, a->l0.me8x8[1].i_ref );
            x264_macroblock_cache_ref( h, 0, 2, 2, 2, a->l0.me8x8[2].i_ref );
            x264_macroblock_cache_ref( h, 2, 2, 2, 2, a->l0.me8x8[3].i_ref );
            /* FIXME: In the 8x8 blocks where RDO isn't run, the NNZ values used for context selection
             * for future blocks are those left over from previous RDO calls. */
            for( int i = 0; i < 4; i++ )
            {
                int costs[4] = {a->l0.i_cost4x4[i], a->l0.i_cost8x4[i], a->l0.i_cost4x8[i], a->l0.me8x8[i].cost};
                int sub8x8_thresh = a->b_early_terminate ? X264_MIN4( costs[0], costs[1], costs[2], costs[3] ) * 5 / 4 : COST_MAX;
                int subtype, btype = D_L0_8x8;
                uint64_t bcost = COST_MAX64;
                for( subtype = D_L0_4x4; subtype <= D_L0_8x8; subtype++ )
                {
                    uint64_t cost;
                    if( costs[subtype] > sub8x8_thresh )
                        continue;
                    h->mb.i_sub_partition[i] = subtype;
                    x264_mb_cache_mv_p8x8( h, a, i );
                    if( subtype == btype )
                        continue;
                    cost = x264_rd_cost_part( h, a->i_lambda2, i<<2, PIXEL_8x8 );
                    COPY2_IF_LT( bcost, cost, btype, subtype );
                }
                if( h->mb.i_sub_partition[i] != btype )
                {
                    h->mb.i_sub_partition[i] = btype;
                    x264_mb_cache_mv_p8x8( h, a, i );
                }
            }
        }
        else
            x264_analyse_update_cache( h, a );
        a->l0.i_cost8x8 = x264_rd_cost_mb( h, a->i_lambda2 );
    }
    else
        a->l0.i_cost8x8 = COST_MAX;
}

/* Rate-distortion optimal QP selection.
 * FIXME: More than half of the benefit of this function seems to be
 * in the way it improves the coding of chroma DC (by decimating or
 * finding a better way to code a single DC coefficient.)
 * There must be a more efficient way to get that portion of the benefit
 * without doing full QP-RD, but RD-decimation doesn't seem to do the
 * trick. */
static  void x264_mb_analyse_qp_rd( x264_t *h, x264_mb_analysis_t *a )
{
    int bcost, cost, failures, prevcost, origcost;
    int orig_qp = h->mb.i_qp, bqp = h->mb.i_qp;
    int last_qp_tried = 0;
    origcost = bcost = x264_rd_cost_mb( h, a->i_lambda2 );
    int origcbp = h->mb.cbp[h->mb.i_mb_xy];

    /* If CBP is already zero, don't raise the quantizer any higher. */
    for( int direction = origcbp ? 1 : -1; direction >= -1; direction-=2 )
    {
        /* Without psy-RD, require monotonicity when moving quant away from previous
         * macroblock's quant; allow 1 failure when moving quant towards previous quant.
         * With psy-RD, allow 1 failure when moving quant away from previous quant,
         * allow 2 failures when moving quant towards previous quant.
         * Psy-RD generally seems to result in more chaotic RD score-vs-quantizer curves. */
        int threshold = (!!h->mb.i_psy_rd);
        /* Raise the threshold for failures if we're moving towards the last QP. */
        if( ( h->mb.i_last_qp < orig_qp && direction == -1 ) ||
            ( h->mb.i_last_qp > orig_qp && direction ==  1 ) )
            threshold++;
        h->mb.i_qp = orig_qp;
        failures = 0;
        prevcost = origcost;

        /* If the current QP results in an empty CBP, it's highly likely that lower QPs
         * (up to a point) will too.  So, jump down to where the threshold will kick in
         * and check the QP there.  If the CBP is still empty, skip the main loop.
         * If it isn't empty, we would have ended up having to check this QP anyways,
         * so as long as we store it for later lookup, we lose nothing. */
        int already_checked_qp = -1;
        int already_checked_cost = COST_MAX;
        if( direction == -1 )
        {
            if( !origcbp )
            {
                h->mb.i_qp = X264_MAX( h->mb.i_qp - threshold - 1, SPEC_QP( h->param.rc.i_qp_min ) );
                h->mb.i_chroma_qp = h->chroma_qp_table[h->mb.i_qp];
                already_checked_cost = x264_rd_cost_mb( h, a->i_lambda2 );
                if( !h->mb.cbp[h->mb.i_mb_xy] )
                {
                    /* If our empty-CBP block is lower QP than the last QP,
                     * the last QP almost surely doesn't have a CBP either. */
                    if( h->mb.i_last_qp > h->mb.i_qp )
                        last_qp_tried = 1;
                    break;
                }
                already_checked_qp = h->mb.i_qp;
                h->mb.i_qp = orig_qp;
            }
        }

        h->mb.i_qp += direction;
        while( h->mb.i_qp >= h->param.rc.i_qp_min && h->mb.i_qp <= SPEC_QP( h->param.rc.i_qp_max ) )
        {
            if( h->mb.i_last_qp == h->mb.i_qp )
                last_qp_tried = 1;
            if( h->mb.i_qp == already_checked_qp )
                cost = already_checked_cost;
            else
            {
                h->mb.i_chroma_qp = h->chroma_qp_table[h->mb.i_qp];
                cost = x264_rd_cost_mb( h, a->i_lambda2 );
                COPY2_IF_LT( bcost, cost, bqp, h->mb.i_qp );
            }

            /* We can't assume that the costs are monotonic over QPs.
             * Tie case-as-failure seems to give better results. */
            if( cost < prevcost )
                failures = 0;
            else
                failures++;
            prevcost = cost;

            if( failures > threshold )
                break;
            if( direction == 1 && !h->mb.cbp[h->mb.i_mb_xy] )
                break;
            h->mb.i_qp += direction;
        }
    }

    /* Always try the last block's QP. */
    if( !last_qp_tried )
    {
        h->mb.i_qp = h->mb.i_last_qp;
        h->mb.i_chroma_qp = h->chroma_qp_table[h->mb.i_qp];
        cost = x264_rd_cost_mb( h, a->i_lambda2 );
        COPY2_IF_LT( bcost, cost, bqp, h->mb.i_qp );
    }

    h->mb.i_qp = bqp;
    h->mb.i_chroma_qp = h->chroma_qp_table[h->mb.i_qp];

}

/*****************************************************************************
 * x264_macroblock_analyse:
 *****************************************************************************/
void x264_macroblock_analyse( x264_t *h )
{
    x264_mb_analysis_t analysis;
    int i_cost = COST_MAX;

    h->mb.i_qp = x264_ratecontrol_mb_qp( h );
    /* If the QP of this MB is within 1 of the previous MB, code the same QP as the previous MB,
     * to lower the bit cost of the qp_delta.  Don't do this if QPRD is enabled. */
    if( h->param.rc.i_aq_mode && h->param.analyse.i_subpel_refine < 10 )
        h->mb.i_qp = abs(h->mb.i_qp - h->mb.i_last_qp) == 1 ? h->mb.i_last_qp : h->mb.i_qp;

    x264_mb_analyse_init( h, &analysis, h->mb.i_qp );

    /*--------------------------- Do the analysis ---------------------------*/
    if( h->sh.i_type == SLICE_TYPE_I )
    {
intra_analysis:
        if( analysis.i_mbrd )
            x264_mb_init_fenc_cache( h, analysis.i_mbrd >= 2 );
        x264_mb_analyse_intra( h, &analysis, COST_MAX );
        if( analysis.i_mbrd )
            x264_intra_rd( h, &analysis, COST_MAX );

        i_cost = analysis.i_satd_i16x16;
        h->mb.i_type = I_16x16;
        COPY2_IF_LT( i_cost, analysis.i_satd_i4x4, h->mb.i_type, I_4x4 );
        if( analysis.i_satd_pcm < i_cost )
            h->mb.i_type = I_PCM;

        else if( analysis.i_mbrd >= 2 )
            x264_intra_rd_refine( h, &analysis );
    }
    else if( h->sh.i_type == SLICE_TYPE_P )
    {
        int b_skip = 0;

        analysis.b_try_skip = 0;
        if( analysis.b_force_intra )
        {
            if( !h->param.analyse.b_psy )
            {
                x264_mb_analyse_init_qp( h, &analysis, X264_MAX( h->mb.i_qp - h->mb.ip_offset, h->param.rc.i_qp_min ) );
                goto intra_analysis;
            }
        }
        else
        {
            int skip_invalid = 0;
            /* Fast P_SKIP detection */
            if( h->param.analyse.b_fast_pskip )
            {
                if( skip_invalid )
                    // FIXME don't need to check this if the reference frame is done
                    {}
                else if( h->param.analyse.i_subpel_refine >= 3 )
                    analysis.b_try_skip = 1;
                else if( h->mb.i_mb_type_left[0] == P_SKIP ||
                         h->mb.i_mb_type_top == P_SKIP ||
                         h->mb.i_mb_type_topleft == P_SKIP ||
                         h->mb.i_mb_type_topright == P_SKIP )
                    b_skip = x264_macroblock_probe_pskip( h );
            }
        }

        if( b_skip )
        {
            h->mb.i_type = P_SKIP;
            h->mb.i_partition = D_16x16;
            /* Set up MVs for future predictors */
            for( int i = 0; i < h->mb.pic.i_fref; i++ )
                M32( h->mb.mvr[i][h->mb.i_mb_xy] ) = 0;
        }
        else
        {
            const unsigned int flags = h->param.analyse.inter;
            int i_type;
            int i_partition;
            int i_satd_inter, i_satd_intra;

            x264_mb_analyse_load_costs( h, &analysis );

            x264_mb_analyse_inter_p16x16( h, &analysis );

            if( h->mb.i_type == P_SKIP )
            {
                for( int i = 1; i < h->mb.pic.i_fref; i++ )
                    M32( h->mb.mvr[i][h->mb.i_mb_xy] ) = 0;
                return;
            }

            if( flags & X264_ANALYSE_PSUB16x16 )
            {
                x264_mb_analyse_inter_p8x8( h, &analysis );
            }

            /* Select best inter mode */
            i_type = P_L0;
            i_partition = D_16x16;
            i_cost = analysis.l0.me16x16.cost;

            if( ( flags & X264_ANALYSE_PSUB16x16 ) && (!analysis.b_early_terminate ||
                analysis.l0.i_cost8x8 < analysis.l0.me16x16.cost) )
            {
                i_type = P_8x8;
                i_partition = D_8x8;
                i_cost = analysis.l0.i_cost8x8;

                /* Do sub 8x8 */
                if( flags & X264_ANALYSE_PSUB8x8 )
                {
                    for( int i = 0; i < 4; i++ )
                    {
                        x264_mb_analyse_inter_p4x4( h, &analysis, i );
                        int i_thresh8x4 = analysis.l0.me4x4[i][1].cost_mv + analysis.l0.me4x4[i][2].cost_mv;
                        if( !analysis.b_early_terminate || analysis.l0.i_cost4x4[i] < analysis.l0.me8x8[i].cost + i_thresh8x4 )
                        {
                            int i_cost8x8 = analysis.l0.i_cost4x4[i];
                            h->mb.i_sub_partition[i] = D_L0_4x4;

                            x264_mb_analyse_inter_p8x4( h, &analysis, i );
                            COPY2_IF_LT( i_cost8x8, analysis.l0.i_cost8x4[i],
                                         h->mb.i_sub_partition[i], D_L0_8x4 );

                            x264_mb_analyse_inter_p4x8( h, &analysis, i );
                            COPY2_IF_LT( i_cost8x8, analysis.l0.i_cost4x8[i],
                                         h->mb.i_sub_partition[i], D_L0_4x8 );

                            i_cost += i_cost8x8 - analysis.l0.me8x8[i].cost;
                        }
                        x264_mb_cache_mv_p8x8( h, &analysis, i );
                    }
                    analysis.l0.i_cost8x8 = i_cost;
                }
            }

            /* Now do 16x8/8x16 */
            int i_thresh16x8 = analysis.l0.me8x8[1].cost_mv + analysis.l0.me8x8[2].cost_mv;
            if( ( flags & X264_ANALYSE_PSUB16x16 ) && (!analysis.b_early_terminate ||
                analysis.l0.i_cost8x8 < analysis.l0.me16x16.cost + i_thresh16x8) )
            {
                int i_avg_mv_ref_cost = (analysis.l0.me8x8[2].cost_mv + analysis.l0.me8x8[2].i_ref_cost
                                      + analysis.l0.me8x8[3].cost_mv + analysis.l0.me8x8[3].i_ref_cost + 1) >> 1;
                analysis.i_cost_est16x8[1] = analysis.i_satd8x8[2] + analysis.i_satd8x8[3] + i_avg_mv_ref_cost;

                x264_mb_analyse_inter_p16x8( h, &analysis, i_cost );
                COPY3_IF_LT( i_cost, analysis.l0.i_cost16x8, i_type, P_L0, i_partition, D_16x8 );

                i_avg_mv_ref_cost = (analysis.l0.me8x8[1].cost_mv + analysis.l0.me8x8[1].i_ref_cost
                                  + analysis.l0.me8x8[3].cost_mv + analysis.l0.me8x8[3].i_ref_cost + 1) >> 1;
                analysis.i_cost_est8x16[1] = analysis.i_satd8x8[1] + analysis.i_satd8x8[3] + i_avg_mv_ref_cost;

                x264_mb_analyse_inter_p8x16( h, &analysis, i_cost );
                COPY3_IF_LT( i_cost, analysis.l0.i_cost8x16, i_type, P_L0, i_partition, D_8x16 );
            }

            h->mb.i_partition = i_partition;

            /* refine qpel */
            //FIXME mb_type costs?
            if( analysis.i_mbrd || !h->mb.i_subpel_refine )
            {
                /* refine later */
            }
            else if( i_partition == D_16x16 )
            {
                x264_me_refine_qpel( h, &analysis.l0.me16x16 );
                i_cost = analysis.l0.me16x16.cost;
            }
            else if( i_partition == D_16x8 )
            {
                x264_me_refine_qpel( h, &analysis.l0.me16x8[0] );
                x264_me_refine_qpel( h, &analysis.l0.me16x8[1] );
                i_cost = analysis.l0.me16x8[0].cost + analysis.l0.me16x8[1].cost;
            }
            else if( i_partition == D_8x16 )
            {
                x264_me_refine_qpel( h, &analysis.l0.me8x16[0] );
                x264_me_refine_qpel( h, &analysis.l0.me8x16[1] );
                i_cost = analysis.l0.me8x16[0].cost + analysis.l0.me8x16[1].cost;
            }
            else if( i_partition == D_8x8 )
            {
                i_cost = 0;
                for( int i8x8 = 0; i8x8 < 4; i8x8++ )
                {
                    switch( h->mb.i_sub_partition[i8x8] )
                    {
                        case D_L0_8x8:
                            x264_me_refine_qpel( h, &analysis.l0.me8x8[i8x8] );
                            i_cost += analysis.l0.me8x8[i8x8].cost;
                            break;
                        case D_L0_8x4:
                            x264_me_refine_qpel( h, &analysis.l0.me8x4[i8x8][0] );
                            x264_me_refine_qpel( h, &analysis.l0.me8x4[i8x8][1] );
                            i_cost += analysis.l0.me8x4[i8x8][0].cost +
                                      analysis.l0.me8x4[i8x8][1].cost;
                            break;
                        case D_L0_4x8:
                            x264_me_refine_qpel( h, &analysis.l0.me4x8[i8x8][0] );
                            x264_me_refine_qpel( h, &analysis.l0.me4x8[i8x8][1] );
                            i_cost += analysis.l0.me4x8[i8x8][0].cost +
                                      analysis.l0.me4x8[i8x8][1].cost;
                            break;

                        case D_L0_4x4:
                            x264_me_refine_qpel( h, &analysis.l0.me4x4[i8x8][0] );
                            x264_me_refine_qpel( h, &analysis.l0.me4x4[i8x8][1] );
                            x264_me_refine_qpel( h, &analysis.l0.me4x4[i8x8][2] );
                            x264_me_refine_qpel( h, &analysis.l0.me4x4[i8x8][3] );
                            i_cost += analysis.l0.me4x4[i8x8][0].cost +
                                      analysis.l0.me4x4[i8x8][1].cost +
                                      analysis.l0.me4x4[i8x8][2].cost +
                                      analysis.l0.me4x4[i8x8][3].cost;
                            break;
                        default:
                            x264_log( h, X264_LOG_ERROR, "internal error (!8x8 && !4x4)\n" );
                            break;
                    }
                }
            }

            if( h->mb.b_chroma_me )
            {
                {
                    x264_mb_analyse_intra_chroma( h, &analysis );
                    x264_mb_analyse_intra( h, &analysis, i_cost - analysis.i_satd_chroma );
                }
                analysis.i_satd_i16x16 += analysis.i_satd_chroma;
                analysis.i_satd_i4x4   += analysis.i_satd_chroma;
            }
            else
                x264_mb_analyse_intra( h, &analysis, i_cost );

            i_satd_inter = i_cost;
            i_satd_intra = X264_MIN( analysis.i_satd_i16x16,
                                      analysis.i_satd_i4x4 );

            if( analysis.i_mbrd )
            {
                x264_mb_analyse_p_rd( h, &analysis, X264_MIN(i_satd_inter, i_satd_intra) );
                i_type = P_L0;
                i_partition = D_16x16;
                i_cost = analysis.l0.i_rd16x16;
                COPY2_IF_LT( i_cost, analysis.l0.i_cost16x8, i_partition, D_16x8 );
                COPY2_IF_LT( i_cost, analysis.l0.i_cost8x16, i_partition, D_8x16 );
                COPY3_IF_LT( i_cost, analysis.l0.i_cost8x8, i_partition, D_8x8, i_type, P_8x8 );
                h->mb.i_type = i_type;
                h->mb.i_partition = i_partition;
                x264_intra_rd( h, &analysis, i_satd_inter * 5/4 + 1 );
            }

            COPY2_IF_LT( i_cost, analysis.i_satd_i16x16, i_type, I_16x16 );
            COPY2_IF_LT( i_cost, analysis.i_satd_i4x4, i_type, I_4x4 );
            COPY2_IF_LT( i_cost, analysis.i_satd_pcm, i_type, I_PCM );

            h->mb.i_type = i_type;

            if( analysis.b_force_intra && !IS_INTRA(i_type) )
            {
                /* Intra masking: copy fdec to fenc and re-encode the block as intra in order to make it appear as if
                 * it was an inter block. */
                x264_analyse_update_cache( h, &analysis );
                x264_macroblock_encode( h );
                for( int p = 0; p < (1); p++ )
                    h->mc.copy[PIXEL_16x16]( h->mb.pic.p_fenc[p], FENC_STRIDE, h->mb.pic.p_fdec[p], FDEC_STRIDE, 16 );
                {
                    int height = 16 >> CHROMA_V_SHIFT;
                    h->mc.copy[PIXEL_8x8]  ( h->mb.pic.p_fenc[1], FENC_STRIDE, h->mb.pic.p_fdec[1], FDEC_STRIDE, height );
                    h->mc.copy[PIXEL_8x8]  ( h->mb.pic.p_fenc[2], FENC_STRIDE, h->mb.pic.p_fdec[2], FDEC_STRIDE, height );
                }
                x264_mb_analyse_init_qp( h, &analysis, X264_MAX( h->mb.i_qp - h->mb.ip_offset, h->param.rc.i_qp_min ) );
                goto intra_analysis;
            }

            if( analysis.i_mbrd >= 2 && h->mb.i_type != I_PCM )
            {
                if( IS_INTRA( h->mb.i_type ) )
                {
                    x264_intra_rd_refine( h, &analysis );
                }
                else if( i_partition == D_16x16 )
                {
                    x264_macroblock_cache_ref( h, 0, 0, 4, 4, analysis.l0.me16x16.i_ref );
                    analysis.l0.me16x16.cost = i_cost;
                    x264_me_refine_qpel_rd( h, &analysis.l0.me16x16, analysis.i_lambda2, 0 );
                }
                else if( i_partition == D_16x8 )
                {
                    h->mb.i_sub_partition[0] = h->mb.i_sub_partition[1] =
                    h->mb.i_sub_partition[2] = h->mb.i_sub_partition[3] = D_L0_8x8;
                    x264_macroblock_cache_ref( h, 0, 0, 4, 2, analysis.l0.me16x8[0].i_ref );
                    x264_macroblock_cache_ref( h, 0, 2, 4, 2, analysis.l0.me16x8[1].i_ref );
                    x264_me_refine_qpel_rd( h, &analysis.l0.me16x8[0], analysis.i_lambda2, 0 );
                    x264_me_refine_qpel_rd( h, &analysis.l0.me16x8[1], analysis.i_lambda2, 8 );
                }
                else if( i_partition == D_8x16 )
                {
                    h->mb.i_sub_partition[0] = h->mb.i_sub_partition[1] =
                    h->mb.i_sub_partition[2] = h->mb.i_sub_partition[3] = D_L0_8x8;
                    x264_macroblock_cache_ref( h, 0, 0, 2, 4, analysis.l0.me8x16[0].i_ref );
                    x264_macroblock_cache_ref( h, 2, 0, 2, 4, analysis.l0.me8x16[1].i_ref );
                    x264_me_refine_qpel_rd( h, &analysis.l0.me8x16[0], analysis.i_lambda2, 0 );
                    x264_me_refine_qpel_rd( h, &analysis.l0.me8x16[1], analysis.i_lambda2, 4 );
                }
                else if( i_partition == D_8x8 )
                {
                    x264_analyse_update_cache( h, &analysis );
                    for( int i8x8 = 0; i8x8 < 4; i8x8++ )
                    {
                        if( h->mb.i_sub_partition[i8x8] == D_L0_8x8 )
                        {
                            x264_me_refine_qpel_rd( h, &analysis.l0.me8x8[i8x8], analysis.i_lambda2, i8x8*4 );
                        }
                        else if( h->mb.i_sub_partition[i8x8] == D_L0_8x4 )
                        {
                            x264_me_refine_qpel_rd( h, &analysis.l0.me8x4[i8x8][0], analysis.i_lambda2, i8x8*4+0 );
                            x264_me_refine_qpel_rd( h, &analysis.l0.me8x4[i8x8][1], analysis.i_lambda2, i8x8*4+2 );
                        }
                        else if( h->mb.i_sub_partition[i8x8] == D_L0_4x8 )
                        {
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x8[i8x8][0], analysis.i_lambda2, i8x8*4+0 );
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x8[i8x8][1], analysis.i_lambda2, i8x8*4+1 );
                        }
                        else if( h->mb.i_sub_partition[i8x8] == D_L0_4x4 )
                        {
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x4[i8x8][0], analysis.i_lambda2, i8x8*4+0 );
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x4[i8x8][1], analysis.i_lambda2, i8x8*4+1 );
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x4[i8x8][2], analysis.i_lambda2, i8x8*4+2 );
                            x264_me_refine_qpel_rd( h, &analysis.l0.me4x4[i8x8][3], analysis.i_lambda2, i8x8*4+3 );
                        }
                    }
                }
            }
        }
    }

    x264_analyse_update_cache( h, &analysis );

    /* In rare cases we can end up qpel-RDing our way back to a larger partition size
     * without realizing it.  Check for this and account for it if necessary. */
    if( analysis.i_mbrd >= 2 )
    {
        /* Don't bother with bipred or 8x8-and-below, the odds are incredibly low. */
        if( h->mb.i_type == P_L0 && h->mb.i_partition != D_16x16 &&
            M32( &h->mb.cache.mv[x264_scan8[0]] ) == M32( &h->mb.cache.mv[x264_scan8[12]] ) &&
            h->mb.cache.ref[x264_scan8[0]] == h->mb.cache.ref[x264_scan8[12]] )
                h->mb.i_partition = D_16x16;
    }


    if( analysis.i_mbrd == 3 && !IS_SKIP(h->mb.i_type) )
        x264_mb_analyse_qp_rd( h, &analysis );

    h->mb.b_noise_reduction = h->mb.b_noise_reduction || (!!h->param.analyse.i_noise_reduction && !IS_INTRA( h->mb.i_type ));

    if( h->mb.b_noise_reduction )
        h->mb.i_skip_intra = 0;
}

/*-------------------- Update MB from the analysis ----------------------*/
static void x264_analyse_update_cache( x264_t *h, x264_mb_analysis_t *a  )
{
    switch( h->mb.i_type )
    {
        case I_4x4:
            for( int i = 0; i < 16; i++ )
                h->mb.cache.intra4x4_pred_mode[x264_scan8[i]] = a->i_predict4x4[i];

            x264_mb_analyse_intra_chroma( h, a );
            break;
        case I_16x16:
            h->mb.i_intra16x16_pred_mode = a->i_predict16x16;
            x264_mb_analyse_intra_chroma( h, a );
            break;

        case I_PCM:
            break;

        case P_L0:
            switch( h->mb.i_partition )
            {
                case D_16x16:
                    x264_macroblock_cache_ref( h, 0, 0, 4, 4, a->l0.me16x16.i_ref );
                    x264_macroblock_cache_mv_ptr( h, 0, 0, 4, 4, a->l0.me16x16.mv );
                    break;

                case D_16x8:
                    x264_macroblock_cache_ref( h, 0, 0, 4, 2, a->l0.me16x8[0].i_ref );
                    x264_macroblock_cache_ref( h, 0, 2, 4, 2, a->l0.me16x8[1].i_ref );
                    x264_macroblock_cache_mv_ptr( h, 0, 0, 4, 2, a->l0.me16x8[0].mv );
                    x264_macroblock_cache_mv_ptr( h, 0, 2, 4, 2, a->l0.me16x8[1].mv );
                    break;

                case D_8x16:
                    x264_macroblock_cache_ref( h, 0, 0, 2, 4, a->l0.me8x16[0].i_ref );
                    x264_macroblock_cache_ref( h, 2, 0, 2, 4, a->l0.me8x16[1].i_ref );
                    x264_macroblock_cache_mv_ptr( h, 0, 0, 2, 4, a->l0.me8x16[0].mv );
                    x264_macroblock_cache_mv_ptr( h, 2, 0, 2, 4, a->l0.me8x16[1].mv );
                    break;

                default:
                    x264_log( h, X264_LOG_ERROR, "internal error P_L0 and partition=%d\n", h->mb.i_partition );
                    break;
            }
            break;

        case P_8x8:
            x264_macroblock_cache_ref( h, 0, 0, 2, 2, a->l0.me8x8[0].i_ref );
            x264_macroblock_cache_ref( h, 2, 0, 2, 2, a->l0.me8x8[1].i_ref );
            x264_macroblock_cache_ref( h, 0, 2, 2, 2, a->l0.me8x8[2].i_ref );
            x264_macroblock_cache_ref( h, 2, 2, 2, 2, a->l0.me8x8[3].i_ref );
            for( int i = 0; i < 4; i++ )
                x264_mb_cache_mv_p8x8( h, a, i );
            break;

        case P_SKIP:
        {
            h->mb.i_partition = D_16x16;
            x264_macroblock_cache_ref( h, 0, 0, 4, 4, 0 );
            x264_macroblock_cache_mv_ptr( h, 0, 0, 4, 4, h->mb.cache.pskip_mv );
            break;
        }

        default: /* the rest of the B types */
            x264_log( h, X264_LOG_ERROR, "internal error (invalid MB type)\n" );
            break;
    }
}

#include "slicetype.c"

