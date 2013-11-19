/*****************************************************************************
 * mvpred.c: motion vector prediction
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *          Laurent Aimar <fenrir@via.ecp.fr>
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

void x264_mb_predict_mv( x264_t *h, int idx, int i_width, int16_t mvp[2] )
{
    const int i8 = x264_scan8[idx];
    const int i_ref= h->mb.cache.ref[i8];
    int     i_refa = h->mb.cache.ref[i8 - 1];
    int16_t *mv_a  = h->mb.cache.mv[i8 - 1];
    int     i_refb = h->mb.cache.ref[i8 - 8];
    int16_t *mv_b  = h->mb.cache.mv[i8 - 8];
    int     i_refc = h->mb.cache.ref[i8 - 8 + i_width];
    int16_t *mv_c  = h->mb.cache.mv[i8 - 8 + i_width];

    // Partitions not yet reached in scan order are unavailable.
    if( (idx&3) >= 2 + (i_width&1) || i_refc == -2 )
    {
        i_refc = h->mb.cache.ref[i8 - 8 - 1];
        mv_c   = h->mb.cache.mv[i8 - 8 - 1];
    }
    if( h->mb.i_partition == D_16x8 )
    {
        if( idx == 0 )
        {
            if( i_refb == i_ref )
            {
                CP32( mvp, mv_b );
                return;
            }
        }
        else
        {
            if( i_refa == i_ref )
            {
                CP32( mvp, mv_a );
                return;
            }
        }
    }
    else if( h->mb.i_partition == D_8x16 )
    {
        if( idx == 0 )
        {
            if( i_refa == i_ref )
            {
                CP32( mvp, mv_a );
                return;
            }
        }
        else
        {
            if( i_refc == i_ref )
            {
                CP32( mvp, mv_c );
                return;
            }
        }
    }

    int i_count = (i_refa == i_ref) + (i_refb == i_ref) + (i_refc == i_ref);

    if( i_count > 1 )
    {
median:
        x264_median_mv( mvp, mv_a, mv_b, mv_c );
    }
    else if( i_count == 1 )
    {
        if( i_refa == i_ref )
            CP32( mvp, mv_a );
        else if( i_refb == i_ref )
            CP32( mvp, mv_b );
        else
            CP32( mvp, mv_c );
    }
    else if( i_refb == -2 && i_refc == -2 && i_refa != -2 )
        CP32( mvp, mv_a );
    else
        goto median;
}

void x264_mb_predict_mv_16x16( x264_t *h, int i_ref, int16_t mvp[2] )
{
    int     i_refa = h->mb.cache.ref[X264_SCAN8_0 - 1];
    int16_t *mv_a  = h->mb.cache.mv[X264_SCAN8_0 - 1];
    int     i_refb = h->mb.cache.ref[X264_SCAN8_0 - 8];
    int16_t *mv_b  = h->mb.cache.mv[X264_SCAN8_0 - 8];
    int     i_refc = h->mb.cache.ref[X264_SCAN8_0 - 8 + 4];
    int16_t *mv_c  = h->mb.cache.mv[X264_SCAN8_0 - 8 + 4];
    if( i_refc == -2 )
    {
        i_refc = h->mb.cache.ref[X264_SCAN8_0 - 8 - 1];
        mv_c   = h->mb.cache.mv[X264_SCAN8_0 - 8 - 1];
    }

    int i_count = (i_refa == i_ref) + (i_refb == i_ref) + (i_refc == i_ref);

    if( i_count > 1 )
    {
median:
        x264_median_mv( mvp, mv_a, mv_b, mv_c );
    }
    else if( i_count == 1 )
    {
        if( i_refa == i_ref )
            CP32( mvp, mv_a );
        else if( i_refb == i_ref )
            CP32( mvp, mv_b );
        else
            CP32( mvp, mv_c );
    }
    else if( i_refb == -2 && i_refc == -2 && i_refa != -2 )
        CP32( mvp, mv_a );
    else
        goto median;
}


void x264_mb_predict_mv_pskip( x264_t *h, int16_t mv[2] )
{
    int     i_refa = h->mb.cache.ref[X264_SCAN8_0 - 1];
    int     i_refb = h->mb.cache.ref[X264_SCAN8_0 - 8];
    int16_t *mv_a  = h->mb.cache.mv[X264_SCAN8_0 - 1];
    int16_t *mv_b  = h->mb.cache.mv[X264_SCAN8_0 - 8];

    if( i_refa == -2 || i_refb == -2 ||
        !( i_refa | M32( mv_a ) ) ||
        !( i_refb | M32( mv_b ) ) )
    {
        M32( mv ) = 0;
    }
    else
        x264_mb_predict_mv_16x16( h, 0, mv );
}

/* This just improves encoder performance, it's not part of the spec */
void x264_mb_predict_mv_ref16x16( x264_t *h, int i_ref, int16_t mvc[9][2], int *i_mvc )
{
    int16_t (*mvr)[2] = h->mb.mvr[i_ref];
    int i = 0;

#define SET_MVP(mvp) \
    { \
        CP32( mvc[i], mvp ); \
        i++; \
    }

#define SET_IMVP(xy) \
    if( xy >= 0 ) \
    { \
        int shift = 1 - h->mb.field[xy]; \
        int16_t *mvp = h->mb.mvr[i_ref<<1>>shift][xy]; \
        mvc[i][0] = mvp[0]; \
        mvc[i][1] = mvp[1]<<1>>shift; \
        i++; \
    }

    if( i_ref == 0 && h->frames.b_have_lowres )
    {
        int idx = h->fenc->i_frame-h->fref[0]->i_frame-1;
        if( idx <= 0 )
        {
            int16_t (*lowres_mv)[2] = h->fenc->lowres_mvs[0][idx];
            if( lowres_mv[0][0] != 0x7fff )
            {
                M32( mvc[i] ) = (M32( lowres_mv[h->mb.i_mb_xy] )*2)&0xfffeffff;
                i++;
            }
        }
    }

    /* spatial predictors */
    {
        SET_MVP( mvr[h->mb.i_mb_left_xy[0]] );
        SET_MVP( mvr[h->mb.i_mb_top_xy] );
        SET_MVP( mvr[h->mb.i_mb_topleft_xy] );
        SET_MVP( mvr[h->mb.i_mb_topright_xy] );
    }
#undef SET_IMVP
#undef SET_MVP

    /* temporal predictors */
    if( h->fref[0]->i_ref > 0 )
    {
        x264_frame_t *l0 = h->fref[0];
        int field = h->mb.i_mb_y&1;
        int curpoc = h->fdec->i_poc + h->fdec->i_delta_poc[field];
        int refpoc = h->fref[i_ref]->i_poc;
        refpoc += l0->i_delta_poc[field^(i_ref&1)];

#define SET_TMVP( dx, dy ) \
        { \
            int mb_index = h->mb.i_mb_xy + dx + dy*h->mb.i_mb_stride; \
            int scale = (curpoc - refpoc) * l0->inv_ref_poc; \
            mvc[i][0] = (l0->mv16x16[mb_index][0]*scale + 128) >> 8; \
            mvc[i][1] = (l0->mv16x16[mb_index][1]*scale + 128) >> 8; \
            i++; \
        }

        SET_TMVP(0,0);
        if( h->mb.i_mb_x < h->mb.i_mb_width-1 )
            SET_TMVP(1,0);
        if( h->mb.i_mb_y < h->mb.i_mb_height-1 )
            SET_TMVP(0,1);
#undef SET_TMVP
    }

    *i_mvc = i;
}
