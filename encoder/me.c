/*****************************************************************************
 * me.c: motion estimation
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
 *          Laurent Aimar <fenrir@via.ecp.fr>
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
#include "macroblock.h"
#include "me.h"

/* presets selected from good points on the speed-vs-quality curve of several test videos
 * subpel_iters[i_subpel_refine] = { refine_hpel, refine_qpel, me_hpel, me_qpel }
 * where me_* are the number of EPZS iterations run on all candidate block types,
 * and refine_* are run only on the winner.
 * the subme=8,9 values are much higher because any amount of satd search makes
 * up its time by reducing the number of qpel-rd iterations. */
static const uint8_t subpel_iterations[][4] =
   {{0,0,0,0},
    {1,1,0,0},
    {0,1,1,0},
    {0,2,1,0},
    {0,2,1,1},
    {0,2,1,2},
    {0,0,2,2},
    {0,0,2,2},
    {0,0,4,10},
    {0,0,4,10},
    {0,0,4,10},
    {0,0,4,10}};

/* (x-1)%6 */
static const uint8_t mod6m1[8] = {5,0,1,2,3,4,5,0};
/* radius 2 hexagon. repeated entries are to avoid having to compute mod6 every time. */
static const int8_t hex2[8][2] = {{-1,-2}, {-2,0}, {-1,2}, {1,2}, {2,0}, {1,-2}, {-1,-2}, {-2,0}};
static const int8_t square1[9][2] = {{0,0}, {0,-1}, {0,1}, {-1,0}, {1,0}, {-1,-1}, {-1,1}, {1,-1}, {1,1}};

static void refine_subpel( x264_t *h, x264_me_t *m, int hpel_iters, int qpel_iters, int *p_halfpel_thresh, int b_refine_qpel );

#define BITS_MVD( mx, my )\
    (p_cost_mvx[(mx)<<2] + p_cost_mvy[(my)<<2])

#define COST_MV( mx, my )\
do\
{\
    int cost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE,\
                   &p_fref_w[(my)*stride+(mx)], stride )\
             + BITS_MVD(mx,my);\
    COPY3_IF_LT( bcost, cost, bmx, mx, bmy, my );\
} while(0)


#define COST_MV_X4_DIR( m0x, m0y, m1x, m1y, m2x, m2y, m3x, m3y, costs )\
{\
    pixel *pix_base = p_fref_w + bmx + bmy*stride;\
    h->pixf.fpelcmp_x4[i_pixel]( p_fenc,\
        pix_base + (m0x) + (m0y)*stride,\
        pix_base + (m1x) + (m1y)*stride,\
        pix_base + (m2x) + (m2y)*stride,\
        pix_base + (m3x) + (m3y)*stride,\
        stride, costs );\
    (costs)[0] += BITS_MVD( bmx+(m0x), bmy+(m0y) );\
    (costs)[1] += BITS_MVD( bmx+(m1x), bmy+(m1y) );\
    (costs)[2] += BITS_MVD( bmx+(m2x), bmy+(m2y) );\
    (costs)[3] += BITS_MVD( bmx+(m3x), bmy+(m3y) );\
}



#define FPEL(mv) (((mv)+2)>>2) /* Convert subpel MV to fullpel with rounding... */
#define SPEL(mv) ((mv)<<2)     /* ... and the reverse. */
#define SPELx2(mv) (SPEL(mv)&0xFFFCFFFC) /* for two packed MVs */

void x264_me_search_ref( x264_t *h, x264_me_t *m, int16_t (*mvc)[2], int i_mvc, int *p_halfpel_thresh )
{
    const int i_pixel = m->i_pixel;
    const int stride = m->i_stride[0];
    int i_me_range = h->param.analyse.i_me_range;
    int bmx, bmy, bcost = COST_MAX;
    int pmx, pmy;
    pixel *p_fenc = m->p_fenc[0];
    pixel *p_fref_w = m->p_fref_w;
    ALIGNED_ARRAY_8( int16_t, mvc_temp,[16],[2] );

    ALIGNED_ARRAY_16( int, costs,[16] );

    int mv_x_min = h->mb.mv_limit_fpel[0][0];
    int mv_y_min = h->mb.mv_limit_fpel[0][1];
    int mv_x_max = h->mb.mv_limit_fpel[1][0];
    int mv_y_max = h->mb.mv_limit_fpel[1][1];
/* Special version of pack to allow shortcuts in CHECK_MVRANGE */
#define pack16to32_mask2(mx,my) ((mx<<16)|(my&0x7FFF))
    uint32_t mv_min = pack16to32_mask2( -mv_x_min, -mv_y_min );
    uint32_t mv_max = pack16to32_mask2( mv_x_max, mv_y_max )|0x8000;
    uint32_t pmv = 0;

#define CHECK_MVRANGE(mx,my) (!(((pack16to32_mask2(mx,my) + mv_min) | (mv_max - pack16to32_mask2(mx,my))) & 0x80004000))

    const uint16_t *p_cost_mvx = m->p_cost_mv - m->mvp[0];
    const uint16_t *p_cost_mvy = m->p_cost_mv - m->mvp[1];

    /* Try extra predictors if provided.  If subme >= 3, check subpel predictors,
     * otherwise round them to fullpel. */
    {
        /* Calculate and check the fullpel MVP first */
        bmx = pmx = x264_clip3( FPEL(m->mvp[0]), mv_x_min, mv_x_max );
        bmy = pmy = x264_clip3( FPEL(m->mvp[1]), mv_y_min, mv_y_max );
        pmv = pack16to32_mask( bmx, bmy );

        /* Because we are rounding the predicted motion vector to fullpel, there will be
         * an extra MV cost in 15 out of 16 cases.  However, when the predicted MV is
         * chosen as the best predictor, it is often the case that the subpel search will
         * result in a vector at or next to the predicted motion vector.  Therefore, we omit
         * the cost of the MV from the rounded MVP to avoid unfairly biasing against use of
         * the predicted motion vector.
         *
         * Disclaimer: this is a post-hoc rationalization for why this hack works. */
        bcost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE, &p_fref_w[bmy*stride+bmx], stride );

        if( i_mvc > 0 )
        {
            /* Like in subme>=3, except we also round the candidates to fullpel. */
            int valid_mvcs = x264_predictor_roundclip( mvc_temp+2, mvc, i_mvc, h->mb.mv_limit_fpel, pmv );
            if( valid_mvcs > 0 )
            {
                int i = 1, cost;
                M32( mvc_temp[1] ) = pmv;
                bcost <<= 4;
                do
                {
                    int mx = mvc_temp[i+1][0];
                    int my = mvc_temp[i+1][1];
                    cost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE, &p_fref_w[my*stride+mx], stride ) + BITS_MVD( mx, my );
                    COPY1_IF_LT( bcost, (cost << 4) + i );
                } while( ++i <= valid_mvcs );
                bmx = mvc_temp[(bcost&15)+1][0];
                bmy = mvc_temp[(bcost&15)+1][1];
                bcost >>= 4;
            }
        }

        /* Same as above, except the condition is simpler. */
        if( pmv )
            COST_MV( 0, 0 );
    }

    switch( h->mb.i_me_method )
    {
        case X264_ME_DIA:
        {
            /* diamond search, radius 1 */
            bcost <<= 4;
            int i = i_me_range;
            do
            {
                COST_MV_X4_DIR( 0,-1, 0,1, -1,0, 1,0, costs );
                COPY1_IF_LT( bcost, (costs[0]<<4)+1 );
                COPY1_IF_LT( bcost, (costs[1]<<4)+3 );
                COPY1_IF_LT( bcost, (costs[2]<<4)+4 );
                COPY1_IF_LT( bcost, (costs[3]<<4)+12 );
                if( !(bcost&15) )
                    break;
                bmx -= (bcost<<28)>>30;
                bmy -= (bcost<<30)>>30;
                bcost &= ~15;
            } while( --i && CHECK_MVRANGE(bmx, bmy) );
            bcost >>= 4;
            break;
        }

        break;
    }

    /* -> qpel mv */
    uint32_t bmv = pack16to32_mask(bmx,bmy);
    uint32_t bmv_spel = SPELx2(bmv);
    {
        m->cost_mv = p_cost_mvx[bmx<<2] + p_cost_mvy[bmy<<2];
        m->cost = bcost;
        /* compute the real cost */
        if( bmv == pmv ) m->cost += m->cost_mv;
        M32( m->mv ) = bmv_spel;
    }

}
#undef COST_MV

void x264_me_refine_qpel( x264_t *h, x264_me_t *m )
{
    int hpel = subpel_iterations[h->mb.i_subpel_refine][0];
    int qpel = subpel_iterations[h->mb.i_subpel_refine][1];

    if( m->i_pixel <= PIXEL_8x8 )
        m->cost -= m->i_ref_cost;

    refine_subpel( h, m, hpel, qpel, NULL, 1 );
}

void x264_me_refine_qpel_refdupe( x264_t *h, x264_me_t *m, int *p_halfpel_thresh )
{
    refine_subpel( h, m, 0, X264_MIN( 2, subpel_iterations[h->mb.i_subpel_refine][3] ), p_halfpel_thresh, 0 );
}

#define COST_MV_SATD( mx, my, dir ) \
if( b_refine_qpel || (dir^1) != odir ) \
{ \
    intptr_t stride = 16; \
    pixel *src = h->mc.get_ref( pix, &stride, &m->p_fref[0], m->i_stride[0], mx, my, bw, bh ); \
    int cost = h->pixf.mbcmp_unaligned[i_pixel]( m->p_fenc[0], FENC_STRIDE, src, stride ) \
             + p_cost_mvx[ mx ] + p_cost_mvy[ my ]; \
    if( b_chroma_me && cost < bcost ) \
    { \
        { \
            h->mc.mc_chroma( pix, pix+8, 16, m->p_fref[4], m->i_stride[1], \
                             mx, 2*(my+mvy_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift ); \
            cost += h->pixf.mbcmp[chromapix]( m->p_fenc[1], FENC_STRIDE, pix, 16 ); \
            if( cost < bcost ) \
            { \
                cost += h->pixf.mbcmp[chromapix]( m->p_fenc[2], FENC_STRIDE, pix+8, 16 ); \
            } \
        } \
    } \
    COPY4_IF_LT( bcost, cost, bmx, mx, bmy, my, bdir, dir ); \
}

static void refine_subpel( x264_t *h, x264_me_t *m, int hpel_iters, int qpel_iters, int *p_halfpel_thresh, int b_refine_qpel )
{
    const int bw = x264_pixel_size[m->i_pixel].w;
    const int bh = x264_pixel_size[m->i_pixel].h;
    const uint16_t *p_cost_mvx = m->p_cost_mv - m->mvp[0];
    const uint16_t *p_cost_mvy = m->p_cost_mv - m->mvp[1];
    const int i_pixel = m->i_pixel;
    const int b_chroma_me = h->mb.b_chroma_me && (i_pixel <= PIXEL_8x8);
    int chromapix = h->luma2chroma_pixel[i_pixel];
    int chroma_v_shift = CHROMA_V_SHIFT;
    int mvy_offset = chroma_v_shift & 0 & m->i_ref ? (h->mb.i_mb_y & 1)*4 - 2 : 0;

    ALIGNED_ARRAY_N( pixel, pix,[64*18] ); // really 17x17x2, but round up for alignment
    ALIGNED_ARRAY_16( int, costs,[4] );

    int bmx = m->mv[0];
    int bmy = m->mv[1];
    int bcost = m->cost;
    int odir = -1, bdir;

    /* halfpel diamond search */
    if( hpel_iters )
    {
        /* try the subpel component of the predicted mv */
        if( h->mb.i_subpel_refine < 3 )
        {
            int mx = x264_clip3( m->mvp[0], h->mb.mv_min_spel[0]+2, h->mb.mv_max_spel[0]-2 );
            int my = x264_clip3( m->mvp[1], h->mb.mv_min_spel[1]+2, h->mb.mv_max_spel[1]-2 );
            if( (mx-bmx)|(my-bmy) ){
                intptr_t stride = 16;
                pixel *src = h->mc.get_ref( pix, &stride, m->p_fref, m->i_stride[0], mx, my, bw, bh );
                int cost = h->pixf.fpelcmp[i_pixel]( m->p_fenc[0], FENC_STRIDE, src, stride )
                         + p_cost_mvx[ mx ] + p_cost_mvy[ my ];
                COPY3_IF_LT( bcost, cost, bmx, mx, bmy, my );
            }
        }

        bcost <<= 6;
        for( int i = hpel_iters; i > 0; i-- )
        {
            int omx = bmx, omy = bmy;
            intptr_t stride = 64; // candidates are either all hpel or all qpel, so one stride is enough
            pixel *src0, *src1, *src2, *src3;
            src0 = h->mc.get_ref( pix,    &stride, m->p_fref, m->i_stride[0], omx, omy-2, bw, bh+1 );
            src2 = h->mc.get_ref( pix+32, &stride, m->p_fref, m->i_stride[0], omx-2, omy, bw+4, bh );
            src1 = src0 + stride;
            src3 = src2 + 1;
            h->pixf.fpelcmp_x4[i_pixel]( m->p_fenc[0], src0, src1, src2, src3, stride, costs );
            costs[0] += p_cost_mvx[omx  ] + p_cost_mvy[omy-2];
            costs[1] += p_cost_mvx[omx  ] + p_cost_mvy[omy+2];
            costs[2] += p_cost_mvx[omx-2] + p_cost_mvy[omy  ];
            costs[3] += p_cost_mvx[omx+2] + p_cost_mvy[omy  ];
            COPY1_IF_LT( bcost, (costs[0]<<6)+2 );
            COPY1_IF_LT( bcost, (costs[1]<<6)+6 );
            COPY1_IF_LT( bcost, (costs[2]<<6)+16 );
            COPY1_IF_LT( bcost, (costs[3]<<6)+48 );
            if( !(bcost&63) )
                break;
            bmx -= (bcost<<26)>>29;
            bmy -= (bcost<<29)>>29;
            bcost &= ~63;
        }
        bcost >>= 6;
    }

    if( !b_refine_qpel && (h->pixf.mbcmp_unaligned[0] != h->pixf.fpelcmp[0] || b_chroma_me) )
    {
        bcost = COST_MAX;
        COST_MV_SATD( bmx, bmy, -1 );
    }

    /* early termination when examining multiple reference frames */
    if( p_halfpel_thresh )
    {
        if( (bcost*7)>>3 > *p_halfpel_thresh )
        {
            m->cost = bcost;
            m->mv[0] = bmx;
            m->mv[1] = bmy;
            // don't need cost_mv
            return;
        }
        else if( bcost < *p_halfpel_thresh )
            *p_halfpel_thresh = bcost;
    }

    /* quarterpel diamond search */
    if( h->mb.i_subpel_refine != 1 )
    {
        bdir = -1;
        for( int i = qpel_iters; i > 0; i-- )
        {
            if( bmy <= h->mb.mv_min_spel[1] || bmy >= h->mb.mv_max_spel[1] || bmx <= h->mb.mv_min_spel[0] || bmx >= h->mb.mv_max_spel[0] )
                break;
            odir = bdir;
            int omx = bmx, omy = bmy;
            COST_MV_SATD( omx, omy - 1, 0 );
            COST_MV_SATD( omx, omy + 1, 1 );
            COST_MV_SATD( omx - 1, omy, 2 );
            COST_MV_SATD( omx + 1, omy, 3 );
            if( (bmx == omx) & (bmy == omy) )
                break;
        }
    }
    /* Special simplified case for subme=1 */
    else if( bmy > h->mb.mv_min_spel[1] && bmy < h->mb.mv_max_spel[1] && bmx > h->mb.mv_min_spel[0] && bmx < h->mb.mv_max_spel[0] )
    {
        int omx = bmx, omy = bmy;
        /* We have to use mc_luma because all strides must be the same to use fpelcmp_x4 */
        h->mc.mc_luma( pix   , 64, m->p_fref, m->i_stride[0], omx, omy-1, bw, bh );
        h->mc.mc_luma( pix+16, 64, m->p_fref, m->i_stride[0], omx, omy+1, bw, bh );
        h->mc.mc_luma( pix+32, 64, m->p_fref, m->i_stride[0], omx-1, omy, bw, bh );
        h->mc.mc_luma( pix+48, 64, m->p_fref, m->i_stride[0], omx+1, omy, bw, bh );
        h->pixf.fpelcmp_x4[i_pixel]( m->p_fenc[0], pix, pix+16, pix+32, pix+48, 64, costs );
        costs[0] += p_cost_mvx[omx  ] + p_cost_mvy[omy-1];
        costs[1] += p_cost_mvx[omx  ] + p_cost_mvy[omy+1];
        costs[2] += p_cost_mvx[omx-1] + p_cost_mvy[omy  ];
        costs[3] += p_cost_mvx[omx+1] + p_cost_mvy[omy  ];
        bcost <<= 4;
        COPY1_IF_LT( bcost, (costs[0]<<4)+1 );
        COPY1_IF_LT( bcost, (costs[1]<<4)+3 );
        COPY1_IF_LT( bcost, (costs[2]<<4)+4 );
        COPY1_IF_LT( bcost, (costs[3]<<4)+12 );
        bmx -= (bcost<<28)>>30;
        bmy -= (bcost<<30)>>30;
        bcost >>= 4;
    }

    m->cost = bcost;
    m->mv[0] = bmx;
    m->mv[1] = bmy;
    m->cost_mv = p_cost_mvx[bmx] + p_cost_mvy[bmy];
}

#define BIME_CACHE( dx, dy, list )\
{\
    x264_me_t *m = m##list;\
    int i = 4 + 3*dx + dy;\
    int mvx = bm##list##x+dx;\
    int mvy = bm##list##y+dy;\
    stride[0][list][i] = bw;\
    src[0][list][i] = h->mc.get_ref( pixy_buf[list][i], &stride[0][list][i], &m->p_fref[0],\
                                     m->i_stride[0], mvx, mvy, bw, bh );\
    if( rd )\
    {\
            h->mc.mc_chroma( pixu_buf[list][i], pixv_buf[list][i], 8, m->p_fref[4], m->i_stride[1],\
                             mvx, 2*(mvy+mv##list##y_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift );\
    }\
}

#define SATD_THRESH(cost) (cost+(cost>>4))

/* Don't unroll the BIME_CACHE loop. I couldn't find any way to force this
 * other than making its iteration count not a compile-time constant. */
int x264_iter_kludge = 0;

#undef COST_MV_SATD
#define COST_MV_SATD( mx, my, dst, avoid_mvp ) \
{ \
    if( !avoid_mvp || !(mx == pmx && my == pmy) ) \
    { \
        h->mc.mc_luma( pix, FDEC_STRIDE, m->p_fref, m->i_stride[0], mx, my, bw, bh ); \
        dst = h->pixf.mbcmp[i_pixel]( m->p_fenc[0], FENC_STRIDE, pix, FDEC_STRIDE ) \
            + p_cost_mvx[mx] + p_cost_mvy[my]; \
        COPY1_IF_LT( bsatd, dst ); \
    } \
    else \
        dst = COST_MAX; \
}

#define COST_MV_RD( mx, my, satd, do_dir, mdir ) \
{ \
    if( satd <= SATD_THRESH(bsatd) ) \
    { \
        uint64_t cost; \
        M32( cache_mv ) = pack16to32_mask(mx,my); \
        if( m->i_pixel <= PIXEL_8x8 ) \
        { \
            h->mc.mc_chroma( pixu, pixv, FDEC_STRIDE, m->p_fref[4], m->i_stride[1], \
                             mx, 2*(my+mvy_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift ); \
        } \
        cost = x264_rd_cost_part( h, i_lambda2, i4, m->i_pixel ); \
        COPY4_IF_LT( bcost, cost, bmx, mx, bmy, my, dir, do_dir?mdir:dir ); \
    } \
}

void x264_me_refine_qpel_rd( x264_t *h, x264_me_t *m, int i_lambda2, int i4/*, int i_list */)
{
    int16_t *cache_mv = h->mb.cache.mv[x264_scan8[i4]];
    const uint16_t *p_cost_mvx, *p_cost_mvy;
    const int bw = x264_pixel_size[m->i_pixel].w;
    const int bh = x264_pixel_size[m->i_pixel].h;
    const int i_pixel = m->i_pixel;
    int chroma_v_shift = CHROMA_V_SHIFT;
    int mvy_offset = chroma_v_shift & 0 & m->i_ref ? (h->mb.i_mb_y & 1)*4 - 2 : 0;

    uint64_t bcost = COST_MAX64;
    int bmx = m->mv[0];
    int bmy = m->mv[1];
    int omx, omy, pmx, pmy;
    int satd, bsatd;
    int dir = -2;
    int i8 = i4>>2;
    uint16_t amvd;

    pixel *pix  = &h->mb.pic.p_fdec[0][block_idx_xy_fdec[i4]];
    pixel *pixu, *pixv;
    {
        pixu = &h->mb.pic.p_fdec[1][(i8>>1)*(8*FDEC_STRIDE>>chroma_v_shift)+(i8&1)*4];
        pixv = &h->mb.pic.p_fdec[2][(i8>>1)*(8*FDEC_STRIDE>>chroma_v_shift)+(i8&1)*4];
    }

    h->mb.b_skip_mc = 1;

    if( m->i_pixel != PIXEL_16x16 && i4 != 0 )
        x264_mb_predict_mv( h, i4, bw>>2, m->mvp );
    pmx = m->mvp[0];
    pmy = m->mvp[1];
    p_cost_mvx = m->p_cost_mv - pmx;
    p_cost_mvy = m->p_cost_mv - pmy;
    COST_MV_SATD( bmx, bmy, bsatd, 0 );
    if( m->i_pixel != PIXEL_16x16 )
        COST_MV_RD( bmx, bmy, 0, 0, 0 )
    else
        bcost = m->cost;

    /* check the predicted mv */
    if( (bmx != pmx || bmy != pmy)
        && pmx >= h->mb.mv_min_spel[0] && pmx <= h->mb.mv_max_spel[0]
        && pmy >= h->mb.mv_min_spel[1] && pmy <= h->mb.mv_max_spel[1] )
    {
        COST_MV_SATD( pmx, pmy, satd, 0 );
        COST_MV_RD  ( pmx, pmy, satd, 0, 0 );
        /* The hex motion search is guaranteed to not repeat the center candidate,
         * so if pmv is chosen, set the "MV to avoid checking" to bmv instead. */
        if( bmx == pmx && bmy == pmy )
        {
            pmx = m->mv[0];
            pmy = m->mv[1];
        }
    }

    if( bmy < h->mb.mv_min_spel[1] + 3 || bmy > h->mb.mv_max_spel[1] - 3 ||
        bmx < h->mb.mv_min_spel[0] + 3 || bmx > h->mb.mv_max_spel[0] - 3 )
    {
        h->mb.b_skip_mc = 0;
        return;
    }

    /* subpel hex search, same pattern as ME HEX. */
    dir = -2;
    omx = bmx;
    omy = bmy;
    for( int j = 0; j < 6; j++ )
    {
        COST_MV_SATD( omx + hex2[j+1][0], omy + hex2[j+1][1], satd, 1 );
        COST_MV_RD  ( omx + hex2[j+1][0], omy + hex2[j+1][1], satd, 1, j );
    }

    if( dir != -2 )
    {
        /* half hexagon, not overlapping the previous iteration */
        for( int i = 1; i < 10; i++ )
        {
            const int odir = mod6m1[dir+1];
            if( bmy < h->mb.mv_min_spel[1] + 3 ||
                bmy > h->mb.mv_max_spel[1] - 3 )
                break;
            dir = -2;
            omx = bmx;
            omy = bmy;
            for( int j = 0; j < 3; j++ )
            {
                COST_MV_SATD( omx + hex2[odir+j][0], omy + hex2[odir+j][1], satd, 1 );
                COST_MV_RD  ( omx + hex2[odir+j][0], omy + hex2[odir+j][1], satd, 1, odir-1+j );
            }
            if( dir == -2 )
                break;
        }
    }

    /* square refine, same pattern as ME HEX. */
    omx = bmx;
    omy = bmy;
    for( int i = 0; i < 8; i++ )
    {
        COST_MV_SATD( omx + square1[i+1][0], omy + square1[i+1][1], satd, 1 );
        COST_MV_RD  ( omx + square1[i+1][0], omy + square1[i+1][1], satd, 0, 0 );
    }

    m->cost = bcost;
    m->mv[0] = bmx;
    m->mv[1] = bmy;
    x264_macroblock_cache_mv ( h, block_idx_x[i4], block_idx_y[i4], bw>>2, bh>>2, pack16to32_mask(bmx, bmy) );
    amvd = pack8to16( X264_MIN(abs(bmx - m->mvp[0]),66), X264_MIN(abs(bmy - m->mvp[1]),66) );
    x264_macroblock_cache_mvd( h, block_idx_x[i4], block_idx_y[i4], bw>>2, bh>>2, amvd );
    h->mb.b_skip_mc = 0;
}
