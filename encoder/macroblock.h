/*****************************************************************************
 * macroblock.h: macroblock encoding
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
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

#ifndef X264_ENCODER_MACROBLOCK_H
#define X264_ENCODER_MACROBLOCK_H

#include "common/macroblock.h"

extern const int x264_lambda2_tab[QP_MAX_MAX+1];
extern const uint16_t x264_lambda_tab[QP_MAX_MAX+1];

void x264_rdo_init( void );

int x264_macroblock_probe_skip( x264_t *h, int b_bidir );

#define x264_macroblock_probe_pskip( h )\
    x264_macroblock_probe_skip( h, 0 )
#define x264_macroblock_probe_bskip( h )\
    x264_macroblock_probe_skip( h, 1 )

void x264_predict_lossless_4x4( x264_t *h, pixel *p_dst, int p, int idx, int i_mode );
void x264_predict_lossless_16x16( x264_t *h, int p, int i_mode );
void x264_predict_lossless_chroma( x264_t *h, int i_mode );

void x264_macroblock_encode      ( x264_t *h );
void x264_macroblock_write_cavlc ( x264_t *h );

void x264_macroblock_encode_p8x8( x264_t *h, int i8 );
void x264_macroblock_encode_p4x4( x264_t *h, int i4 );
void x264_mb_encode_chroma( x264_t *h, int b_inter, int i_qp );

void x264_noise_reduction_update( x264_t *h );

static  int x264_quant_4x4( x264_t *h, dctcoef dct[16], int i_qp, int ctx_block_cat, int b_intra, int p, int idx )
{
    int i_quant_cat = b_intra ? (p?CQM_4IC:CQM_4IY) : (p?CQM_4PC:CQM_4PY);
    if( h->mb.b_noise_reduction )
        h->quantf.denoise_dct( dct, h->nr_residual_sum[0+!!p*2], h->nr_offset[0+!!p*2], 16 );
    return h->quantf.quant_4x4( dct, h->quant4_mf[i_quant_cat][i_qp], h->quant4_bias[i_quant_cat][i_qp] );
}


#define STORE_8x8_NNZ( p, idx, nz )\
do\
{\
    M16( &h->mb.cache.non_zero_count[x264_scan8[p*16+idx*4]+0] ) = (nz) * 0x0101;\
    M16( &h->mb.cache.non_zero_count[x264_scan8[p*16+idx*4]+8] ) = (nz) * 0x0101;\
} while(0)

#define CLEAR_16x16_NNZ( p ) \
do\
{\
    M32( &h->mb.cache.non_zero_count[x264_scan8[16*p+ 0]] ) = 0;\
    M32( &h->mb.cache.non_zero_count[x264_scan8[16*p+ 2]] ) = 0;\
    M32( &h->mb.cache.non_zero_count[x264_scan8[16*p+ 8]] ) = 0;\
    M32( &h->mb.cache.non_zero_count[x264_scan8[16*p+10]] ) = 0;\
} while(0)

/* A special for loop that iterates branchlessly over each set
 * bit in a 4-bit input. */
#define FOREACH_BIT(idx,start,mask) for( int idx = start, msk = mask, skip; msk && (skip = x264_ctz_4bit(msk), idx += skip, msk >>= skip+1, 1); idx++ )

static inline void x264_mb_encode_i4x4( x264_t *h, int p, int idx, int i_qp, int i_mode, int b_predict )
{
    int nz;
    pixel *p_src = &h->mb.pic.p_fenc[p][block_idx_xy_fenc[idx]];
    pixel *p_dst = &h->mb.pic.p_fdec[p][block_idx_xy_fdec[idx]];
    ALIGNED_ARRAY_N( dctcoef, dct4x4,[16] );

    if( b_predict )
    {
        if( h->mb.b_lossless )
            x264_predict_lossless_4x4( h, p_dst, p, idx, i_mode );
        else
            h->predict_4x4[i_mode]( p_dst );
    }

    if( h->mb.b_lossless )
    {
        nz = h->zigzagf.sub_4x4( h->dct.luma4x4[p*16+idx], p_src, p_dst );
        h->mb.cache.non_zero_count[x264_scan8[p*16+idx]] = nz;
        h->mb.i_cbp_luma |= nz<<(idx>>2);
        return;
    }

    h->dctf.sub4x4_dct( dct4x4, p_src, p_dst );

    nz = x264_quant_4x4( h, dct4x4, i_qp, ctx_cat_plane[DCT_LUMA_4x4][p], 1, p, idx );
    h->mb.cache.non_zero_count[x264_scan8[p*16+idx]] = nz;
    if( nz )
    {
        h->mb.i_cbp_luma |= 1<<(idx>>2);
        h->zigzagf.scan_4x4( h->dct.luma4x4[p*16+idx], dct4x4 );
        h->quantf.dequant_4x4( dct4x4, h->dequant4_mf[p?CQM_4IC:CQM_4IY], i_qp );
        h->dctf.add4x4_idct( p_dst, dct4x4 );
    }
}

#endif

