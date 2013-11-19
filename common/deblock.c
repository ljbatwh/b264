/*****************************************************************************
 * deblock.c: deblocking
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *          Henrik Gramner <hengar-6@student.ltu.se>
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

#include "common.h"

/* Deblocking filter */
static const uint8_t i_alpha_table[52+12*3] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  4,  4,  5,  6,
     7,  8,  9, 10, 12, 13, 15, 17, 20, 22,
    25, 28, 32, 36, 40, 45, 50, 56, 63, 71,
    80, 90,101,113,127,144,162,182,203,226,
   255,255,
   255,255,255,255,255,255,255,255,255,255,255,255,
};
static const uint8_t i_beta_table[52+12*3] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  2,  2,  3,
     3,  3,  3,  4,  4,  4,  6,  6,  7,  7,
     8,  8,  9,  9, 10, 10, 11, 11, 12, 12,
    13, 13, 14, 14, 15, 15, 16, 16, 17, 17,
    18, 18,
    18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
};
static const int8_t i_tc0_table[52+12*3][4] =
{
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 },
    {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 0 }, {-1, 0, 0, 1 },
    {-1, 0, 0, 1 }, {-1, 0, 0, 1 }, {-1, 0, 0, 1 }, {-1, 0, 1, 1 }, {-1, 0, 1, 1 }, {-1, 1, 1, 1 },
    {-1, 1, 1, 1 }, {-1, 1, 1, 1 }, {-1, 1, 1, 1 }, {-1, 1, 1, 2 }, {-1, 1, 1, 2 }, {-1, 1, 1, 2 },
    {-1, 1, 1, 2 }, {-1, 1, 2, 3 }, {-1, 1, 2, 3 }, {-1, 2, 2, 3 }, {-1, 2, 2, 4 }, {-1, 2, 3, 4 },
    {-1, 2, 3, 4 }, {-1, 3, 3, 5 }, {-1, 3, 4, 6 }, {-1, 3, 4, 6 }, {-1, 4, 5, 7 }, {-1, 4, 5, 8 },
    {-1, 4, 6, 9 }, {-1, 5, 7,10 }, {-1, 6, 8,11 }, {-1, 6, 8,13 }, {-1, 7,10,14 }, {-1, 8,11,16 },
    {-1, 9,12,18 }, {-1,10,13,20 }, {-1,11,15,23 }, {-1,13,17,25 },
    {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 },
    {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 }, {-1,13,17,25 },
};
#define alpha_table(x) i_alpha_table[(x)+24]
#define beta_table(x)  i_beta_table[(x)+24]
#define tc0_table(x)   i_tc0_table[(x)+24]

/* From ffmpeg */
static  void deblock_edge_luma_c( pixel *pix, intptr_t xstride, int alpha, int beta, int8_t tc0 )
{
    int p2 = pix[-3*xstride];
    int p1 = pix[-2*xstride];
    int p0 = pix[-1*xstride];
    int q0 = pix[ 0*xstride];
    int q1 = pix[ 1*xstride];
    int q2 = pix[ 2*xstride];

    if( abs( p0 - q0 ) < alpha && abs( p1 - p0 ) < beta && abs( q1 - q0 ) < beta )
    {
        int tc = tc0;
        int delta;
        if( abs( p2 - p0 ) < beta )
        {
            if( tc0 )
                pix[-2*xstride] = p1 + x264_clip3( (( p2 + ((p0 + q0 + 1) >> 1)) >> 1) - p1, -tc0, tc0 );
            tc++;
        }
        if( abs( q2 - q0 ) < beta )
        {
            if( tc0 )
                pix[ 1*xstride] = q1 + x264_clip3( (( q2 + ((p0 + q0 + 1) >> 1)) >> 1) - q1, -tc0, tc0 );
            tc++;
        }

        delta = x264_clip3( (((q0 - p0 ) << 2) + (p1 - q1) + 4) >> 3, -tc, tc );
        pix[-1*xstride] = x264_clip_pixel( p0 + delta );    /* p0' */
        pix[ 0*xstride] = x264_clip_pixel( q0 - delta );    /* q0' */
    }
}
static  void deblock_luma_c( pixel *pix, intptr_t xstride, intptr_t ystride, int alpha, int beta, int8_t *tc0 )
{
    for( int i = 0; i < 4; i++ )
    {
        if( tc0[i] < 0 )
        {
            pix += 4*ystride;
            continue;
        }
        for( int d = 0; d < 4; d++, pix += ystride )
            deblock_edge_luma_c( pix, xstride, alpha, beta, tc0[i] );
    }
}
static void deblock_h_luma_mbaff_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    for( int d = 0; d < 8; d++, pix += stride )
        deblock_edge_luma_c( pix, 1, alpha, beta, tc0[d>>1] );
}
static void deblock_v_luma_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    deblock_luma_c( pix, stride, 1, alpha, beta, tc0 );
}
static void deblock_h_luma_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    deblock_luma_c( pix, 1, stride, alpha, beta, tc0 );
}

static  void deblock_edge_chroma_c( pixel *pix, intptr_t xstride, int alpha, int beta, int8_t tc )
{
    int p1 = pix[-2*xstride];
    int p0 = pix[-1*xstride];
    int q0 = pix[ 0*xstride];
    int q1 = pix[ 1*xstride];

    if( abs( p0 - q0 ) < alpha && abs( p1 - p0 ) < beta && abs( q1 - q0 ) < beta )
    {
        int delta = x264_clip3( (((q0 - p0 ) << 2) + (p1 - q1) + 4) >> 3, -tc, tc );
        pix[-1*xstride] = x264_clip_pixel( p0 + delta );    /* p0' */
        pix[ 0*xstride] = x264_clip_pixel( q0 - delta );    /* q0' */
    }
}
static  void deblock_chroma_c( pixel *pix, int height, intptr_t xstride, intptr_t ystride, int alpha, int beta, int8_t *tc0 )
{
    for( int i = 0; i < 4; i++ )
    {
        int tc = tc0[i];
        if( tc <= 0 )
        {
            pix += height*ystride;
            continue;
        }
        for( int d = 0; d < height; d++, pix += ystride-2 )
            for( int e = 0; e < 2; e++, pix++ )
                deblock_edge_chroma_c( pix, xstride, alpha, beta, tc0[i] );
    }
}
static void deblock_h_chroma_mbaff_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    deblock_chroma_c( pix, 1, 2, stride, alpha, beta, tc0 );
}
static void deblock_v_chroma_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    deblock_chroma_c( pix, 2, stride, 2, alpha, beta, tc0 );
}
static void deblock_h_chroma_c( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 )
{
    deblock_chroma_c( pix, 2, 2, stride, alpha, beta, tc0 );
}

static  void deblock_edge_luma_intra_c( pixel *pix, intptr_t xstride, int alpha, int beta )
{
    int p2 = pix[-3*xstride];
    int p1 = pix[-2*xstride];
    int p0 = pix[-1*xstride];
    int q0 = pix[ 0*xstride];
    int q1 = pix[ 1*xstride];
    int q2 = pix[ 2*xstride];

    if( abs( p0 - q0 ) < alpha && abs( p1 - p0 ) < beta && abs( q1 - q0 ) < beta )
    {
        if( abs( p0 - q0 ) < ((alpha >> 2) + 2) )
        {
            if( abs( p2 - p0 ) < beta ) /* p0', p1', p2' */
            {
                const int p3 = pix[-4*xstride];
                pix[-1*xstride] = ( p2 + 2*p1 + 2*p0 + 2*q0 + q1 + 4 ) >> 3;
                pix[-2*xstride] = ( p2 + p1 + p0 + q0 + 2 ) >> 2;
                pix[-3*xstride] = ( 2*p3 + 3*p2 + p1 + p0 + q0 + 4 ) >> 3;
            }
            else /* p0' */
                pix[-1*xstride] = ( 2*p1 + p0 + q1 + 2 ) >> 2;
            if( abs( q2 - q0 ) < beta ) /* q0', q1', q2' */
            {
                const int q3 = pix[3*xstride];
                pix[0*xstride] = ( p1 + 2*p0 + 2*q0 + 2*q1 + q2 + 4 ) >> 3;
                pix[1*xstride] = ( p0 + q0 + q1 + q2 + 2 ) >> 2;
                pix[2*xstride] = ( 2*q3 + 3*q2 + q1 + q0 + p0 + 4 ) >> 3;
            }
            else /* q0' */
                pix[0*xstride] = ( 2*q1 + q0 + p1 + 2 ) >> 2;
        }
        else /* p0', q0' */
        {
            pix[-1*xstride] = ( 2*p1 + p0 + q1 + 2 ) >> 2;
            pix[ 0*xstride] = ( 2*q1 + q0 + p1 + 2 ) >> 2;
        }
    }
}
static  void deblock_luma_intra_c( pixel *pix, intptr_t xstride, intptr_t ystride, int alpha, int beta )
{
    for( int d = 0; d < 16; d++, pix += ystride )
        deblock_edge_luma_intra_c( pix, xstride, alpha, beta );
}
static void deblock_h_luma_intra_mbaff_c( pixel *pix, intptr_t ystride, int alpha, int beta )
{
    for( int d = 0; d < 8; d++, pix += ystride )
        deblock_edge_luma_intra_c( pix, 1, alpha, beta );
}
static void deblock_v_luma_intra_c( pixel *pix, intptr_t stride, int alpha, int beta )
{
    deblock_luma_intra_c( pix, stride, 1, alpha, beta );
}
static void deblock_h_luma_intra_c( pixel *pix, intptr_t stride, int alpha, int beta )
{
    deblock_luma_intra_c( pix, 1, stride, alpha, beta );
}

static  void deblock_edge_chroma_intra_c( pixel *pix, intptr_t xstride, int alpha, int beta )
{
    int p1 = pix[-2*xstride];
    int p0 = pix[-1*xstride];
    int q0 = pix[ 0*xstride];
    int q1 = pix[ 1*xstride];

    if( abs( p0 - q0 ) < alpha && abs( p1 - p0 ) < beta && abs( q1 - q0 ) < beta )
    {
        pix[-1*xstride] = (2*p1 + p0 + q1 + 2) >> 2;   /* p0' */
        pix[ 0*xstride] = (2*q1 + q0 + p1 + 2) >> 2;   /* q0' */
    }
}
static  void deblock_chroma_intra_c( pixel *pix, int width, int height, intptr_t xstride, intptr_t ystride, int alpha, int beta )
{
    for( int d = 0; d < height; d++, pix += ystride-2 )
        for( int e = 0; e < width; e++, pix++ )
            deblock_edge_chroma_intra_c( pix, xstride, alpha, beta );
}
static void deblock_h_chroma_intra_mbaff_c( pixel *pix, intptr_t stride, int alpha, int beta )
{
    deblock_chroma_intra_c( pix, 2, 4, 2, stride, alpha, beta );
}
static void deblock_v_chroma_intra_c( pixel *pix, intptr_t stride, int alpha, int beta )
{
    deblock_chroma_intra_c( pix, 1, 16, stride, 2, alpha, beta );
}
static void deblock_h_chroma_intra_c( pixel *pix, intptr_t stride, int alpha, int beta )
{
    deblock_chroma_intra_c( pix, 2, 8, 2, stride, alpha, beta );
}

static void deblock_strength_c( uint8_t nnz[X264_SCAN8_SIZE], int8_t ref[X264_SCAN8_LUMA_SIZE],
                                int16_t mv[X264_SCAN8_LUMA_SIZE][2], uint8_t bs[8][4], int mvy_limit)
{
    {
        int s1 = 8;
        int s2 = 1;
        for( int edge = 0; edge < 4; edge++ )
            for( int i = 0, loc = X264_SCAN8_0+edge*s2; i < 4; i++, loc += s1 )
            {
                int locn = loc - s2;
                if( nnz[loc] || nnz[locn] )
                    bs[edge][i] = 2;
                else if( ref[loc] != ref[locn] ||
                         abs( mv[loc][0] - mv[locn][0] ) >= 4 ||
                         abs( mv[loc][1] - mv[locn][1] ) >= mvy_limit )
                {
                    bs[edge][i] = 1;
                }
                else
                    bs[edge][i] = 0;
            }
    }
}

static  void deblock_edge( x264_t *h, pixel *pix, intptr_t i_stride, uint8_t bS[4], int i_qp,
                                        int a, int b, int b_chroma, x264_deblock_inter_t pf_inter )
{
    int index_a = i_qp + a;
    int index_b = i_qp + b;
    int alpha = alpha_table(index_a) << (BIT_DEPTH-8);
    int beta  = beta_table(index_b) << (BIT_DEPTH-8);
    int8_t tc[4];

    if( !M32(bS) || !alpha || !beta )
        return;

    tc[0] = (tc0_table(index_a)[bS[0]] << (BIT_DEPTH-8)) + b_chroma;
    tc[1] = (tc0_table(index_a)[bS[1]] << (BIT_DEPTH-8)) + b_chroma;
    tc[2] = (tc0_table(index_a)[bS[2]] << (BIT_DEPTH-8)) + b_chroma;
    tc[3] = (tc0_table(index_a)[bS[3]] << (BIT_DEPTH-8)) + b_chroma;

    pf_inter( pix, i_stride, alpha, beta, tc );
}

static  void deblock_edge_intra( x264_t *h, pixel *pix, intptr_t i_stride, uint8_t bS[4], int i_qp,
                                              int a, int b, int b_chroma, x264_deblock_intra_t pf_intra )
{
    int index_a = i_qp + a;
    int index_b = i_qp + b;
    int alpha = alpha_table(index_a) << (BIT_DEPTH-8);
    int beta  = beta_table(index_b) << (BIT_DEPTH-8);

    if( !alpha || !beta )
        return;

    pf_intra( pix, i_stride, alpha, beta );
}

static  void x264_macroblock_cache_load_neighbours_deblock( x264_t *h, int mb_x, int mb_y )
{
    int deblock_on_slice_edges = h->sh.i_disable_deblocking_filter_idc != 2;

    h->mb.i_neighbour = 0;
    h->mb.i_mb_xy = mb_y * h->mb.i_mb_stride + mb_x;
    h->mb.i_mb_top_y = mb_y - 1;
    h->mb.i_mb_top_xy = mb_x + h->mb.i_mb_stride*h->mb.i_mb_top_y;
    h->mb.i_mb_left_xy[1] =
    h->mb.i_mb_left_xy[0] = h->mb.i_mb_xy - 1;

    if( mb_x > 0 && (deblock_on_slice_edges ||
        h->mb.slice_table[h->mb.i_mb_left_xy[0]] == h->mb.slice_table[h->mb.i_mb_xy]) )
        h->mb.i_neighbour |= MB_LEFT;
    if( mb_y > 0 && (deblock_on_slice_edges
        || h->mb.slice_table[h->mb.i_mb_top_xy] == h->mb.slice_table[h->mb.i_mb_xy]) )
        h->mb.i_neighbour |= MB_TOP;
}

void x264_frame_deblock_row( x264_t *h, int mb_y )
{
    int a = h->sh.i_alpha_c0_offset ;
    int b = h->sh.i_beta_offset ;
    int qp_thresh = 15 - X264_MIN( a, b ) - X264_MAX( 0, h->pps->i_chroma_qp_index_offset );
    int stridey   = h->fdec->i_stride[0];
    int strideuv  = h->fdec->i_stride[1];
    int chroma_height = 16 >> CHROMA_V_SHIFT;

    for( int mb_x = 0; mb_x < h->mb.i_mb_width; mb_x += (~0 | mb_y)&1, mb_y ^= 0 )
    {
        x264_macroblock_cache_load_neighbours_deblock( h, mb_x, mb_y );

        int mb_xy = h->mb.i_mb_xy;
        int intra_cur = IS_INTRA( h->mb.type[mb_xy] );
        uint8_t (*bs)[8][4] = h->deblock_strength[mb_y&1][mb_x];

        pixel *pixy = h->fdec->plane[0] + 16*mb_y*stridey  + 16*mb_x;
        pixel *pixuv = h->fdec->plane[1] + chroma_height*mb_y*strideuv + 16*mb_x;

        int stride2y  = stridey;
        int stride2uv = strideuv;
        int qp = h->mb.qp[mb_xy];
        int qpc = h->chroma_qp_table[qp];
        int first_edge_only = (h->mb.partition[mb_xy] == D_16x16 && !h->mb.cbp[mb_xy] && !intra_cur) || qp <= qp_thresh;

        #define FILTER( intra, dir, edge, qp, chroma_qp )\
        do\
        {\
            {\
                deblock_edge##intra( h, pixy + 4*edge*(dir?stride2y:1),\
                                     stride2y, bs[dir][edge], qp, a, b, 0,\
                                     h->loopf.deblock_luma##intra[dir] );\
                if( CHROMA_FORMAT == CHROMA_420 && !(edge & 1) )\
                {\
                    deblock_edge##intra( h, pixuv + edge*(dir?2*stride2uv:4),\
                                         stride2uv, bs[dir][edge], chroma_qp, a, b, 1,\
                                         h->loopf.deblock_chroma##intra[dir] );\
                }\
            }\
        } while(0)

        if( h->mb.i_neighbour & MB_LEFT )
        {
            {
                int luma_qp[2];
                int chroma_qp[2];
                int left_qp[2];
                x264_deblock_inter_t luma_deblock = h->loopf.deblock_luma_mbaff;
                x264_deblock_inter_t chroma_deblock = h->loopf.deblock_chroma_mbaff;
                x264_deblock_intra_t luma_intra_deblock = h->loopf.deblock_luma_intra_mbaff;
                x264_deblock_intra_t chroma_intra_deblock = h->loopf.deblock_chroma_intra_mbaff;
                int c = 1;

                left_qp[0] = h->mb.qp[h->mb.i_mb_left_xy[0]];
                luma_qp[0] = (qp + left_qp[0] + 1) >> 1;
                chroma_qp[0] = (qpc + h->chroma_qp_table[left_qp[0]] + 1) >> 1;
                if( intra_cur || IS_INTRA( h->mb.type[h->mb.i_mb_left_xy[0]] ) )
                {
                    deblock_edge_intra( h, pixy,           2*stridey,  bs[0][0], luma_qp[0],   a, b, 0, luma_intra_deblock );
                    deblock_edge_intra( h, pixuv,          2*strideuv, bs[0][0], chroma_qp[0], a, b, c, chroma_intra_deblock );
                }
                else
                {
                    deblock_edge( h, pixy,           2*stridey,  bs[0][0], luma_qp[0],   a, b, 0, luma_deblock );
                    deblock_edge( h, pixuv,          2*strideuv, bs[0][0], chroma_qp[0], a, b, c, chroma_deblock );
                }

                int offy =  0;
                int offuv = 0;
                left_qp[1] = h->mb.qp[h->mb.i_mb_left_xy[1]];
                luma_qp[1] = (qp + left_qp[1] + 1) >> 1;
                chroma_qp[1] = (qpc + h->chroma_qp_table[left_qp[1]] + 1) >> 1;
                if( intra_cur || IS_INTRA( h->mb.type[h->mb.i_mb_left_xy[1]] ) )
                {
                    deblock_edge_intra( h, pixy           + (stridey<<offy),   2*stridey,  bs[0][4], luma_qp[1],   a, b, 0, luma_intra_deblock );
                    deblock_edge_intra( h, pixuv          + (strideuv<<offuv), 2*strideuv, bs[0][4], chroma_qp[1], a, b, c, chroma_intra_deblock );
                }
                else
                {
                    deblock_edge( h, pixy           + (stridey<<offy),   2*stridey,  bs[0][4], luma_qp[1],   a, b, 0, luma_deblock );
                    deblock_edge( h, pixuv          + (strideuv<<offuv), 2*strideuv, bs[0][4], chroma_qp[1], a, b, c, chroma_deblock );
                }
            }
        }
        if( !first_edge_only )
        {
            FILTER( , 0, 1, qp, qpc );
            FILTER( , 0, 2, qp, qpc );
            FILTER( , 0, 3, qp, qpc );
        }

        if( h->mb.i_neighbour & MB_TOP )
        {
            {
                int mbn_xy = mb_xy - 2 * h->mb.i_mb_stride;

                for( int j = 0; j < 2; j++, mbn_xy += h->mb.i_mb_stride )
                {
                    int qpt = h->mb.qp[mbn_xy];
                    int qp_top = (qp + qpt + 1) >> 1;
                    int qpc_top = (qpc + h->chroma_qp_table[qpt] + 1) >> 1;
                    int intra_top = IS_INTRA( h->mb.type[mbn_xy] );
                    if( intra_cur || intra_top )
                        M32( bs[1][4*j] ) = 0x03030303;

                    // deblock the first horizontal edge of the even rows, then the first horizontal edge of the odd rows
                    deblock_edge( h, pixy      + j*stridey,  2* stridey, bs[1][4*j], qp_top, a, b, 0, h->loopf.deblock_luma[1] );
                    deblock_edge( h, pixuv          + j*strideuv, 2*strideuv, bs[1][4*j], qpc_top, a, b, 1, h->loopf.deblock_chroma[1] );
                }
            }
        }

        if( !first_edge_only )
        {
            FILTER( , 1, 1, qp, qpc );
            FILTER( , 1, 2, qp, qpc );
            FILTER( , 1, 3, qp, qpc );
        }

        #undef FILTER
    }
}

/* For deblock-aware RD.
 * TODO:
 *  deblock macroblock edges
 *  support analysis partitions smaller than 16x16
 *  deblock chroma for 4:2:0/4:2:2
 *  handle duplicate refs correctly
 */
void x264_macroblock_deblock( x264_t *h )
{
    int a = h->sh.i_alpha_c0_offset ;
    int b = h->sh.i_beta_offset ;
    int qp_thresh = 15 - X264_MIN( a, b ) - X264_MAX( 0, h->pps->i_chroma_qp_index_offset );
    int intra_cur = IS_INTRA( h->mb.i_type );
    int qp = h->mb.i_qp;
    if( (h->mb.i_partition == D_16x16 && !h->mb.i_cbp_luma && !intra_cur) || qp <= qp_thresh )
        return;

    uint8_t (*bs)[8][4] = h->mb.cache.deblock_strength;
    if( intra_cur )
    {
        memset( &bs[0][1], 3, 3*4*sizeof(uint8_t) );
        memset( &bs[1][1], 3, 3*4*sizeof(uint8_t) );
    }
    else
        h->loopf.deblock_strength( h->mb.cache.non_zero_count, h->mb.cache.ref, h->mb.cache.mv,
                                   *bs, 4 );

    #define FILTER( edge )\
    do\
    {\
        deblock_edge( h, h->mb.pic.p_fdec[0] + 4*edge*(1),\
                      FDEC_STRIDE, bs[0][edge], qp, a, b, 0,\
                      h->loopf.deblock_luma[0] );\
    } while(0)

                        FILTER( 1 );
                        FILTER( 2 );
                        FILTER( 3 );

    #undef FILTER
}

void x264_deblock_init( x264_deblock_function_t *pf)
{
    pf->deblock_luma[1] = deblock_v_luma_c;
    pf->deblock_luma[0] = deblock_h_luma_c;
    pf->deblock_chroma[1] = deblock_v_chroma_c;
    pf->deblock_h_chroma_420 = deblock_h_chroma_c;
    pf->deblock_luma_intra[1] = deblock_v_luma_intra_c;
    pf->deblock_luma_intra[0] = deblock_h_luma_intra_c;
    pf->deblock_chroma_intra[1] = deblock_v_chroma_intra_c;
    pf->deblock_h_chroma_420_intra = deblock_h_chroma_intra_c;
    pf->deblock_luma_mbaff = deblock_h_luma_mbaff_c;
    pf->deblock_chroma_420_mbaff = deblock_h_chroma_mbaff_c;
    pf->deblock_luma_intra_mbaff = deblock_h_luma_intra_mbaff_c;
    pf->deblock_chroma_420_intra_mbaff = deblock_h_chroma_intra_mbaff_c;
    pf->deblock_strength = deblock_strength_c;
}
