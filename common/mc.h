/*****************************************************************************
 * mc.h: motion compensation
 *****************************************************************************
 * Copyright (C) 2004-2013 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
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

#ifndef X264_MC_H
#define X264_MC_H


/* Do the MC
 * XXX: Only width = 4, 8 or 16 are valid
 * width == 4 -> height == 4 or 8
 * width == 8 -> height == 4 or 8 or 16
 * width == 16-> height == 8 or 16
 * */

typedef struct
{
    void (*mc_luma)( pixel *dst, intptr_t i_dst, pixel **src, intptr_t i_src,
                     int mvx, int mvy, int i_width, int i_height );

    /* may round up the dimensions if they're not a power of 2 */
    pixel* (*get_ref)( pixel *dst, intptr_t *i_dst, pixel **src, intptr_t i_src,
                       int mvx, int mvy, int i_width, int i_height );

    /* mc_chroma may write up to 2 bytes of garbage to the right of dst,
     * so it must be run from left to right. */
    void (*mc_chroma)( pixel *dstu, pixel *dstv, intptr_t i_dst, pixel *src, intptr_t i_src,
                       int mvx, int mvy, int i_width, int i_height );

    void (*avg[12])( pixel *dst,  intptr_t dst_stride, pixel *src1, intptr_t src1_stride,
                     pixel *src2, intptr_t src2_stride );

    /* only 16x16, 8x8, and 4x4 defined */
    void (*copy[7])( pixel *dst, intptr_t dst_stride, pixel *src, intptr_t src_stride, int i_height );
    void (*copy_16x16_unaligned)( pixel *dst, intptr_t dst_stride, pixel *src, intptr_t src_stride, int i_height );

    void (*store_interleave_chroma)( pixel *dst, intptr_t i_dst, pixel *srcu, pixel *srcv, int height );
    void (*load_deinterleave_chroma_fenc)( pixel *dst, pixel *src, intptr_t i_src, int height );
    void (*load_deinterleave_chroma_fdec)( pixel *dst, pixel *src, intptr_t i_src, int height );

    void (*plane_copy)( pixel *dst, intptr_t i_dst, pixel *src, intptr_t i_src, int w, int h );
    void (*plane_copy_interleave)( pixel *dst,  intptr_t i_dst, pixel *srcu, intptr_t i_srcu,
                                   pixel *srcv, intptr_t i_srcv, int w, int h );
    /* may write up to 15 pixels off the end of each plane */
    void (*plane_copy_deinterleave)( pixel *dstu, intptr_t i_dstu, pixel *dstv, intptr_t i_dstv,
                                     pixel *src,  intptr_t i_src, int w, int h );
    void (*plane_copy_deinterleave_rgb)( pixel *dsta, intptr_t i_dsta, pixel *dstb, intptr_t i_dstb,
                                         pixel *dstc, intptr_t i_dstc, pixel *src,  intptr_t i_src, int pw, int w, int h );
    void (*hpel_filter)( pixel *dsth, pixel *dstv, pixel *dstc, pixel *src,
                         intptr_t i_stride, int i_width, int i_height, int16_t *buf );

    void *(*memcpy_aligned)( void *dst, const void *src, size_t n );
    void (*memzero_aligned)( void *dst, size_t n );

    /* successive elimination prefilter */
    void (*integral_init4h)( uint16_t *sum, pixel *pix, intptr_t stride );
    void (*integral_init8h)( uint16_t *sum, pixel *pix, intptr_t stride );
    void (*integral_init4v)( uint16_t *sum8, uint16_t *sum4, intptr_t stride );
    void (*integral_init8v)( uint16_t *sum8, intptr_t stride );

    void (*frame_init_lowres_core)( pixel *src0, pixel *dst0, pixel *dsth, pixel *dstv, pixel *dstc,
                                    intptr_t src_stride, intptr_t dst_stride, int width, int height );
} x264_mc_functions_t;

void x264_mc_init( x264_mc_functions_t *pf );

#endif
