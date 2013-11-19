/*****************************************************************************
 * rdo.c: rate-distortion optimization
 *****************************************************************************
 * Copyright (C) 2005-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
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

/* duplicate all the writer functions, just calculating bit cost
 * instead of writing the bitstream.
 * TODO: use these for fast 1st pass too. */

#define RDO_SKIP_BS 1

/* CAVLC: produces exactly the same bit count as a normal encode */
/* this probably still leaves some unnecessary computations */
#define bs_write1(s,v)     ((s)->i_bits_encoded += 1)
#define bs_write(s,n,v)    ((s)->i_bits_encoded += (n))
#define bs_write_ue(s,v)   ((s)->i_bits_encoded += bs_size_ue(v))
#define bs_write_se(s,v)   ((s)->i_bits_encoded += bs_size_se(v))
#define bs_write_te(s,v,l) ((s)->i_bits_encoded += bs_size_te(v,l))
#define x264_macroblock_write_cavlc  static x264_macroblock_size_cavlc
#include "cavlc.c"

static  uint64_t cached_hadamard( x264_t *h, int size, int x, int y )
{
    static const uint8_t hadamard_shift_x[4] = {4,   4,   3,   3};
    static const uint8_t hadamard_shift_y[4] = {4-0, 3-0, 4-1, 3-1};
    static const uint8_t  hadamard_offset[4] = {0,   1,   3,   5};
    int cache_index = (x >> hadamard_shift_x[size]) + (y >> hadamard_shift_y[size])
                    + hadamard_offset[size];
    uint64_t res = h->mb.pic.fenc_hadamard_cache[cache_index];
    if( res )
        return res - 1;
    else
    {
        pixel *fenc = h->mb.pic.p_fenc[0] + x + y*FENC_STRIDE;
        res = h->pixf.hadamard_ac[size]( fenc, FENC_STRIDE );
        h->mb.pic.fenc_hadamard_cache[cache_index] = res + 1;
        return res;
    }
}

static  int cached_satd( x264_t *h, int size, int x, int y )
{
    static const uint8_t satd_shift_x[3] = {3,   2,   2};
    static const uint8_t satd_shift_y[3] = {2-1, 3-2, 2-2};
    static const uint8_t  satd_offset[3] = {0,   8,   16};
    ALIGNED_16( static pixel zero[16] ) = {0};
    int cache_index = (x >> satd_shift_x[size - PIXEL_8x4]) + (y >> satd_shift_y[size - PIXEL_8x4])
                    + satd_offset[size - PIXEL_8x4];
    int res = h->mb.pic.fenc_satd_cache[cache_index];
    if( res )
        return res - 1;
    else
    {
        pixel *fenc = h->mb.pic.p_fenc[0] + x + y*FENC_STRIDE;
        int dc = h->pixf.sad[size]( fenc, FENC_STRIDE, zero, 0 ) >> 1;
        res = h->pixf.satd[size]( fenc, FENC_STRIDE, zero, 0 ) - dc;
        h->mb.pic.fenc_satd_cache[cache_index] = res + 1;
        return res;
    }
}

/* Psy RD distortion metric: SSD plus "Absolute Difference of Complexities" */
/* SATD and SA8D are used to measure block complexity. */
/* The difference between SATD and SA8D scores are both used to avoid bias from the DCT size.  Using SATD */
/* only, for example, results in overusage of 8x8dct, while the opposite occurs when using SA8D. */

/* FIXME:  Is there a better metric than averaged SATD/SA8D difference for complexity difference? */
/* Hadamard transform is recursive, so a SATD+SA8D can be done faster by taking advantage of this fact. */
/* This optimization can also be used in non-RD transform decision. */

static  int ssd_plane( x264_t *h, int size, int p, int x, int y )
{
    ALIGNED_16( static pixel zero[16] ) = {0};
    int satd = 0;
    pixel *fdec = h->mb.pic.p_fdec[p] + x + y*FDEC_STRIDE;
    pixel *fenc = h->mb.pic.p_fenc[p] + x + y*FENC_STRIDE;
    if( p == 0 && h->mb.i_psy_rd )
    {
        /* If the plane is smaller than 8x8, we can't do an SA8D; this probably isn't a big problem. */
        if( size <= PIXEL_8x8 )
        {
            uint64_t fdec_acs = h->pixf.hadamard_ac[size]( fdec, FDEC_STRIDE );
            uint64_t fenc_acs = cached_hadamard( h, size, x, y );
            satd = abs((int32_t)fdec_acs - (int32_t)fenc_acs)
                 + abs((int32_t)(fdec_acs>>32) - (int32_t)(fenc_acs>>32));
            satd >>= 1;
        }
        else
        {
            int dc = h->pixf.sad[size]( fdec, FDEC_STRIDE, zero, 0 ) >> 1;
            satd = abs(h->pixf.satd[size]( fdec, FDEC_STRIDE, zero, 0 ) - dc - cached_satd( h, size, x, y ));
        }
        satd = (satd * h->mb.i_psy_rd * h->mb.i_psy_rd_lambda + 128) >> 8;
    }
    return h->pixf.ssd[size](fenc, FENC_STRIDE, fdec, FDEC_STRIDE) + satd;
}

static  int ssd_mb( x264_t *h )
{
    int chroma_size = h->luma2chroma_pixel[PIXEL_16x16];
    int chroma_ssd = ssd_plane(h, chroma_size, 1, 0, 0) + ssd_plane(h, chroma_size, 2, 0, 0);
    chroma_ssd = ((uint64_t)chroma_ssd * h->mb.i_chroma_lambda2_offset + 128) >> 8;
    return ssd_plane(h, PIXEL_16x16, 0, 0, 0) + chroma_ssd;
}

static int x264_rd_cost_mb( x264_t *h, int i_lambda2 )
{
    int i_ssd;
    int i_bits;
    int type_bak = h->mb.i_type;

    x264_macroblock_encode( h );

    if( h->mb.b_deblock_rdo )
        x264_macroblock_deblock( h );

    i_ssd = ssd_mb( h );

    if( IS_SKIP( h->mb.i_type ) )
    {
        i_bits = (1 * i_lambda2 + 128) >> 8;
    }
    else
    {
        x264_macroblock_size_cavlc( h );
        i_bits = ( h->out.bs.i_bits_encoded * i_lambda2 + 128 ) >> 8;
    }

    h->mb.i_type = type_bak;

    return i_ssd + i_bits;
}

/* partition RD functions use 8 bits more precision to avoid large rounding errors at low QPs */

static uint64_t x264_rd_cost_subpart( x264_t *h, int i_lambda2, int i4, int i_pixel )
{
    uint64_t i_ssd, i_bits;

    x264_macroblock_encode_p4x4( h, i4 );
    if( i_pixel == PIXEL_8x4 )
        x264_macroblock_encode_p4x4( h, i4+1 );
    if( i_pixel == PIXEL_4x8 )
        x264_macroblock_encode_p4x4( h, i4+2 );

    i_ssd = ssd_plane( h, i_pixel, 0, block_idx_x[i4]*4, block_idx_y[i4]*4 );

    i_bits = x264_subpartition_size_cavlc( h, i4, i_pixel );

    return (i_ssd<<8) + i_bits;
}

uint64_t x264_rd_cost_part( x264_t *h, int i_lambda2, int i4, int i_pixel )
{
    uint64_t i_ssd, i_bits;
    int i8 = i4 >> 2;

    if( i_pixel == PIXEL_16x16 )
    {
        int i_cost = x264_rd_cost_mb( h, i_lambda2 );
        return i_cost;
    }

    if( i_pixel > PIXEL_8x8 )
        return x264_rd_cost_subpart( h, i_lambda2, i4, i_pixel );

    h->mb.i_cbp_luma = 0;

    x264_macroblock_encode_p8x8( h, i8 );
    if( i_pixel == PIXEL_16x8 )
        x264_macroblock_encode_p8x8( h, i8+1 );
    if( i_pixel == PIXEL_8x16 )
        x264_macroblock_encode_p8x8( h, i8+2 );

    int ssd_x = 8*(i8&1);
    int ssd_y = 8*(i8>>1);
    i_ssd = ssd_plane( h, i_pixel, 0, ssd_x, ssd_y );
    int chromapix = h->luma2chroma_pixel[i_pixel];
    int chromassd = ssd_plane( h, chromapix, 1, ssd_x>>CHROMA_H_SHIFT, ssd_y>>CHROMA_V_SHIFT )
                  + ssd_plane( h, chromapix, 2, ssd_x>>CHROMA_H_SHIFT, ssd_y>>CHROMA_V_SHIFT );
    i_ssd += ((uint64_t)chromassd * h->mb.i_chroma_lambda2_offset + 128) >> 8;

    i_bits = x264_partition_size_cavlc( h, i8, i_pixel ) * i_lambda2;

    return (i_ssd<<8) + i_bits;
}

static uint64_t x264_rd_cost_i4x4( x264_t *h, int i_lambda2, int i4, int i_mode )
{
    uint64_t i_ssd, i_bits;
    int plane_count = 1;
    int i_qp = h->mb.i_qp;

    for( int p = 0; p < plane_count; p++ )
    {
        x264_mb_encode_i4x4( h, p, i4, i_qp, i_mode, 1 );
        i_qp = h->mb.i_chroma_qp;
    }

    i_ssd = ssd_plane( h, PIXEL_4x4, 0, block_idx_x[i4]*4, block_idx_y[i4]*4 );

    i_bits = x264_partition_i4x4_size_cavlc( h, i4, i_mode ) * i_lambda2;

    return (i_ssd<<8) + i_bits;
}

static uint64_t x264_rd_cost_chroma( x264_t *h, int i_lambda2, int i_mode, int b_dct )
{
    uint64_t i_ssd, i_bits;

    if( b_dct )
        x264_mb_encode_chroma( h, 0, h->mb.i_chroma_qp );

    int chromapix = h->luma2chroma_pixel[PIXEL_16x16];
    i_ssd = ssd_plane( h, chromapix, 1, 0, 0 )
          + ssd_plane( h, chromapix, 2, 0, 0 );

    h->mb.i_chroma_pred_mode = i_mode;

    i_bits = x264_chroma_size_cavlc( h ) * i_lambda2;

    return (i_ssd<<8) + i_bits;
}

/* precalculate the cost of coding various combinations of bits in a single context */
void x264_rdo_init( void )
{

}

typedef struct
{
    uint64_t score;
    int level_idx; // index into level_tree[]
} trellis_node_t;

typedef struct
{
    uint16_t next;
    uint16_t abs_level;
} trellis_level_t;

// TODO:
// save cabac state between blocks?
// use trellis' RD score instead of x264_mb_decimate_score?
// code 8x8 sig/last flags forwards with deadzone and save the contexts at
//   each position?
// change weights when using CQMs?

// possible optimizations:
// make scores fit in 32bit
// save quantized coefs during rd, to avoid a duplicate trellis in the final encode
// if trellissing all MBRD modes, finish SSD calculation so we can skip all of
//   the normal dequant/idct/ssd/cabac

// the unquant_mf here is not the same as dequant_mf:
// in normal operation (dct->quant->dequant->idct) the dct and idct are not
// normalized. quant/dequant absorb those scaling factors.
// in this function, we just do (quant->unquant) and want the output to be
// comparable to the input. so unquant is the direct inverse of quant,
// and uses the dct scaling factors, not the idct ones.

#define SIGN(x,y) ((x^(y >> 31))-(y >> 31))

#define SET_LEVEL(ndst, nsrc, l) {\
    if( sizeof(trellis_level_t) == sizeof(uint32_t) )\
        M32( &level_tree[levels_used] ) = pack16to32( nsrc.level_idx, l );\
    else\
        level_tree[levels_used] = (trellis_level_t){ nsrc.level_idx, l };\
    ndst.level_idx = levels_used;\
    levels_used++;\
}

