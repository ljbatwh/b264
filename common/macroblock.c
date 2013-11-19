/*****************************************************************************
 * macroblock.c: macroblock common functions
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Jason Garrett-Glaser <darkshikari@gmail.com>
 *          Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
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
#include "encoder/me.h"

static NOINLINE void x264_mb_mc_0xywh( x264_t *h, int x, int y, int width, int height )
{
    int i8    = x264_scan8[0]+x+8*y;
    int i_ref = h->mb.cache.ref[i8];
    int mvx   = x264_clip3( h->mb.cache.mv[i8][0], h->mb.mv_min[0], h->mb.mv_max[0] ) + 4*4*x;
    int mvy   = x264_clip3( h->mb.cache.mv[i8][1], h->mb.mv_min[1], h->mb.mv_max[1] ) + 4*4*y;

    h->mc.mc_luma( &h->mb.pic.p_fdec[0][4*y*FDEC_STRIDE+4*x], FDEC_STRIDE,
                   &h->mb.pic.p_fref[i_ref][0], h->mb.pic.i_stride[0],
                   mvx, mvy, 4*width, 4*height );

    {
        int v_shift = CHROMA_V_SHIFT;
        int offset = (4*FDEC_STRIDE>>v_shift)*y + 2*x;
        height = 4*height >> v_shift;

        h->mc.mc_chroma( &h->mb.pic.p_fdec[1][offset],
                         &h->mb.pic.p_fdec[2][offset], FDEC_STRIDE,
                         h->mb.pic.p_fref[i_ref][4], h->mb.pic.i_stride[1],
                         mvx, 2*mvy>>v_shift, 2*width, height );
    }
}



void x264_mb_mc_8x8( x264_t *h, int i8 )
{
    int x = 2*(i8&1);
    int y = 2*(i8>>1);

    if( h->sh.i_type == SLICE_TYPE_P )
    {
        switch( h->mb.i_sub_partition[i8] )
        {
            case D_L0_8x8:
                x264_mb_mc_0xywh( h, x, y, 2, 2 );
                break;
            case D_L0_8x4:
                x264_mb_mc_0xywh( h, x, y+0, 2, 1 );
                x264_mb_mc_0xywh( h, x, y+1, 2, 1 );
                break;
            case D_L0_4x8:
                x264_mb_mc_0xywh( h, x+0, y, 1, 2 );
                x264_mb_mc_0xywh( h, x+1, y, 1, 2 );
                break;
            case D_L0_4x4:
                x264_mb_mc_0xywh( h, x+0, y+0, 1, 1 );
                x264_mb_mc_0xywh( h, x+1, y+0, 1, 1 );
                x264_mb_mc_0xywh( h, x+0, y+1, 1, 1 );
                x264_mb_mc_0xywh( h, x+1, y+1, 1, 1 );
                break;
        }
    }
    else
    {
        int scan8 = x264_scan8[0] + x + 8*y;

        if( h->mb.cache.ref[scan8] >= 0 )
            x264_mb_mc_0xywh( h, x, y, 2, 2 );
    }
}

void x264_mb_mc( x264_t *h )
{
    if( h->mb.i_partition == D_8x8 )
    {
        for( int i = 0; i < 4; i++ )
            x264_mb_mc_8x8( h, i );
    }
    else
    {
        int ref0a = h->mb.cache.ref[x264_scan8[ 0]];
        int ref0b = h->mb.cache.ref[x264_scan8[12]];
        //int ref1a = h->mb.cache.ref[1][x264_scan8[ 0]];
        //int ref1b = h->mb.cache.ref[1][x264_scan8[12]];

        if( h->mb.i_partition == D_16x16 )
        {
            if( ref0a >= 0 )
                x264_mb_mc_0xywh ( h, 0, 0, 4, 4 );
        }
        else if( h->mb.i_partition == D_16x8 )
        {
            if( ref0a >= 0 )
                x264_mb_mc_0xywh ( h, 0, 0, 4, 2 );

            if( ref0b >= 0 )
                x264_mb_mc_0xywh ( h, 0, 2, 4, 2 );
        }
        else if( h->mb.i_partition == D_8x16 )
        {
            if( ref0a >= 0 )
                x264_mb_mc_0xywh ( h, 0, 0, 2, 4 );

            if( ref0b >= 0 )
                x264_mb_mc_0xywh ( h, 2, 0, 2, 4 );
        }
    }
}

int x264_macroblock_cache_allocate( x264_t *h )
{
    int i_mb_count = h->mb.i_mb_count;

    h->mb.i_mb_stride = h->mb.i_mb_width;
    h->mb.i_b8_stride = h->mb.i_mb_width * 2;
    h->mb.i_b4_stride = h->mb.i_mb_width * 4;

    PREALLOC_INIT

    PREALLOC( h->mb.qp, i_mb_count * sizeof(int8_t) );
    PREALLOC( h->mb.cbp, i_mb_count * sizeof(int16_t) );
    PREALLOC( h->mb.slice_table, i_mb_count * sizeof(uint16_t) );

    /* 0 -> 3 top(4), 4 -> 6 : left(3) */
    PREALLOC( h->mb.intra4x4_pred_mode, i_mb_count * 8 * sizeof(int8_t) );

    /* all coeffs */
    PREALLOC( h->mb.non_zero_count, i_mb_count * 48 * sizeof(uint8_t) );

    PREALLOC_END( h->mb.base );

    memset( h->mb.slice_table, -1, i_mb_count * sizeof(uint16_t) );

    return 0;
fail:
    return -1;
}
void x264_macroblock_cache_free( x264_t *h )
{
    x264_free( h->mb.base );
}

int x264_macroblock_thread_allocate( x264_t *h )
{

    for( int i = 0; i < 2; i++ )
        for( int j = 0; j < 2; j++ )
        {
            CHECKED_MALLOC( h->intra_border_backup[i][j], (h->sps->i_mb_width*16+32) * sizeof(pixel) );
            h->intra_border_backup[i][j] += 16;
        }


    /* Allocate scratch buffer */
    int scratch_size = 0;
    {
        int buf_hpel = (h->fdec->i_width[0]+48+32) * sizeof(int16_t);
        int buf_ssim = h->param.analyse.b_ssim * 8 * (h->param.i_width/4+3) * sizeof(int);
        int me_range = X264_MIN(h->param.analyse.i_me_range, h->param.analyse.i_mv_range);
        int buf_tesa = (h->param.analyse.i_me_method >= X264_ME_ESA) *
            ((me_range*2+24) * sizeof(int16_t) + (me_range+4) * (me_range+1) * 4 * sizeof(mvsad_t));
        scratch_size = X264_MAX3( buf_hpel, buf_ssim, buf_tesa );
    }
    scratch_size = scratch_size;
    if( scratch_size )
        CHECKED_MALLOC( h->scratch_buffer, scratch_size );
    else
        h->scratch_buffer = NULL;

    int buf_lookahead_threads = (h->mb.i_mb_height) * sizeof(int) * 2;
    CHECKED_MALLOC( h->scratch_buffer2, buf_lookahead_threads );

    return 0;
fail:
    return -1;
}

void x264_macroblock_thread_free( x264_t *h )
{

    for( int i = 0; i < 2; i++ )
        for( int j = 0; j < ( 2); j++ )
            x264_free( h->intra_border_backup[i][j] - 16 );

    x264_free( h->scratch_buffer );
    x264_free( h->scratch_buffer2 );
}

void x264_macroblock_slice_init( x264_t *h )
{
    h->mb.mv = h->fdec->mv;
    h->mb.mvr[0] = h->fdec->mv16x16;
    h->mb.ref = h->fdec->ref;
    h->mb.type = h->fdec->mb_type;
    h->mb.partition = h->fdec->mb_partition;
    h->mb.field = h->fdec->field;

    h->fdec->i_ref = h->i_ref;
    for( int i = 0; i < h->i_ref; i++ )
        h->fdec->ref_poc[i] = h->fref[i]->i_poc;

    /* init with not available (for top right idx=7,15) */
    memset( h->mb.cache.ref, -2, sizeof( h->mb.cache.ref ) );

    if( h->i_ref > 0 ){
        int curpoc = h->fdec->i_poc + h->fdec->i_delta_poc[0];
        int refpoc = h->fref[0]->i_poc + h->fref[0]->i_delta_poc[0];
        int delta = curpoc - refpoc;

        h->fdec->inv_ref_poc = (256 + delta/2) / delta;
    }

    h->mb.i_neighbour4[6] =
    h->mb.i_neighbour4[9] =
    h->mb.i_neighbour4[12] =
    h->mb.i_neighbour4[14] = MB_LEFT|MB_TOP|MB_TOPLEFT|MB_TOPRIGHT;
    h->mb.i_neighbour4[3] =
    h->mb.i_neighbour4[7] =
    h->mb.i_neighbour4[11] =
    h->mb.i_neighbour4[13] =
    h->mb.i_neighbour4[15] = MB_LEFT|MB_TOP|MB_TOPLEFT;
}

void x264_macroblock_thread_init( x264_t *h )
{
    h->mb.i_me_method = h->param.analyse.i_me_method;
    h->mb.i_subpel_refine = h->param.analyse.i_subpel_refine;
    h->mb.b_chroma_me = h->param.analyse.b_chroma_me &&
                        ((h->sh.i_type == SLICE_TYPE_P && h->mb.i_subpel_refine >= 5) );
    h->mb.b_dct_decimate = (h->param.analyse.b_dct_decimate && h->sh.i_type != SLICE_TYPE_I);
    h->mb.i_mb_prev_xy = -1;

    /*          4:2:0                      4:2:2                      4:4:4
     * fdec            fenc       fdec            fenc       fdec            fenc
     * y y y y y y y   Y Y Y Y    y y y y y y y   Y Y Y Y    y y y y y y y   Y Y Y Y
     * y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y
     * y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y
     * y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y    y Y Y Y Y       Y Y Y Y
     * y Y Y Y Y       U U V V    y Y Y Y Y       U U V V    y Y Y Y Y       U U U U
     * u u u   v v v   U U V V    u u u   v v v   U U V V    u u u u u u u   U U U U
     * u U U   v V V              u U U   v V V   U U V V    u U U U U       U U U U
     * u U U   v V V              u U U   v V V   U U V V    u U U U U       U U U U
     *                            u U U   v V V              u U U U U       V V V V
     *                            u U U   v V V              u U U U U       V V V V
     *                                                       v v v v v v v   V V V V
     *                                                       v V V V V       V V V V
     *                                                       v V V V V
     *                                                       v V V V V
     *                                                       v V V V V
     */
    h->mb.pic.p_fenc[0] = h->mb.pic.fenc_buf;
    h->mb.pic.p_fdec[0] = h->mb.pic.fdec_buf + 2*FDEC_STRIDE;
    h->mb.pic.p_fenc[1] = h->mb.pic.fenc_buf + 16*FENC_STRIDE;
    h->mb.pic.p_fdec[1] = h->mb.pic.fdec_buf + 19*FDEC_STRIDE;
    h->mb.pic.p_fenc[2] = h->mb.pic.fenc_buf + 16*FENC_STRIDE + 8;
    h->mb.pic.p_fdec[2] = h->mb.pic.fdec_buf + 19*FDEC_STRIDE + 16;
}


void x264_copy_column8( pixel *dst, pixel *src )
{
    // input pointers are offset by 4 rows because that's faster (smaller instruction size on x86)
    for( int i = -4; i < 4; i++ )
        dst[i*FDEC_STRIDE] = src[i*FDEC_STRIDE];
}

static void  x264_macroblock_load_pic_pointers( x264_t *h, int mb_x, int mb_y, int i, int b_chroma )
{
    int height = b_chroma ? 16 >> CHROMA_V_SHIFT : 16;
    int i_stride = h->fdec->i_stride[i];
    int i_stride2 = i_stride << 0;
    int i_pix_offset = 16 * mb_x + height * mb_y * i_stride;
    pixel *plane_fdec = &h->fdec->plane[i][i_pix_offset];
    int fdec_idx = !(mb_y&1);
    pixel *intra_fdec = &h->intra_border_backup[fdec_idx][i][mb_x*16];
    int ref_pix_offset[2] = { i_pix_offset, i_pix_offset };
    /* ref_pix_offset[0] references the current field and [1] the opposite field. */
    h->mb.pic.i_stride[i] = i_stride2;
    h->mb.pic.p_fenc_plane[i] = &h->fenc->plane[i][i_pix_offset];
    if( b_chroma )
    {
        h->mc.load_deinterleave_chroma_fenc( h->mb.pic.p_fenc[1], h->mb.pic.p_fenc_plane[1], i_stride2, height );
        memcpy( h->mb.pic.p_fdec[1]-FDEC_STRIDE, intra_fdec, 8*sizeof(pixel) );
        memcpy( h->mb.pic.p_fdec[2]-FDEC_STRIDE, intra_fdec+8, 8*sizeof(pixel) );
        h->mb.pic.p_fdec[1][-FDEC_STRIDE-1] = intra_fdec[-1-8];
        h->mb.pic.p_fdec[2][-FDEC_STRIDE-1] = intra_fdec[-1];
    }
    else
    {
        h->mc.copy[PIXEL_16x16]( h->mb.pic.p_fenc[i], FENC_STRIDE, h->mb.pic.p_fenc_plane[i], i_stride2, 16 );
        memcpy( h->mb.pic.p_fdec[i]-FDEC_STRIDE, intra_fdec, 24*sizeof(pixel) );
        h->mb.pic.p_fdec[i][-FDEC_STRIDE-1] = intra_fdec[-1];
    }
    if( h->mb.b_reencode_mb )
    {
        for( int j = 0; j < height; j++ )
            if( b_chroma )
            {
                h->mb.pic.p_fdec[1][-1+j*FDEC_STRIDE] = plane_fdec[-2+j*i_stride2];
                h->mb.pic.p_fdec[2][-1+j*FDEC_STRIDE] = plane_fdec[-1+j*i_stride2];
            }
            else
                h->mb.pic.p_fdec[i][-1+j*FDEC_STRIDE] = plane_fdec[-1+j*i_stride2];
    }
    pixel *plane_src, **filtered_src;
    for( int j = 0; j < h->mb.pic.i_fref; j++ )
    {
        {
            plane_src = h->fref[j]->plane[i];
            filtered_src = h->fref[j]->filtered[i];
        }
        h->mb.pic.p_fref[j][i*4] = plane_src + ref_pix_offset[j&1];

        if( !b_chroma )
        {
            for( int k = 1; k < 4; k++ )
                h->mb.pic.p_fref[j][i*4+k] = filtered_src[k] + ref_pix_offset[j&1];
        }
    }
}

static const x264_left_table_t left_indices[4] =
{
    /* Current is progressive */
    {{ 4, 4, 5, 5}, { 3,  3,  7,  7}, {16+1, 16+1, 32+1, 32+1}, {0, 0, 1, 1}, {0, 0, 0, 0}},
    {{ 6, 6, 3, 3}, {11, 11, 15, 15}, {16+5, 16+5, 32+5, 32+5}, {2, 2, 3, 3}, {1, 1, 1, 1}},
    /* Current is interlaced */
    {{ 4, 6, 4, 6}, { 3, 11,  3, 11}, {16+1, 16+1, 32+1, 32+1}, {0, 2, 0, 2}, {0, 1, 0, 1}},
    /* Both same */
    {{ 4, 5, 6, 3}, { 3,  7, 11, 15}, {16+1, 16+5, 32+1, 32+5}, {0, 1, 2, 3}, {0, 0, 1, 1}}
};

static void  x264_macroblock_cache_load_neighbours( x264_t *h, int mb_x, int mb_y)
{
    int top_y = mb_y - 1;
    int top = top_y * h->mb.i_mb_stride + mb_x;

    h->mb.i_mb_x = mb_x;
    h->mb.i_mb_y = mb_y;
    h->mb.i_mb_xy = mb_y * h->mb.i_mb_stride + mb_x;
    h->mb.i_b8_xy = 2*(mb_y * h->mb.i_b8_stride + mb_x);
    h->mb.i_b4_xy = 4*(mb_y * h->mb.i_b4_stride + mb_x);
    h->mb.left_b8[0] =
    h->mb.left_b8[1] = -1;
    h->mb.left_b4[0] =
    h->mb.left_b4[1] = -1;
    h->mb.i_neighbour = 0;
    h->mb.i_neighbour_intra = 0;
    h->mb.i_neighbour_frame = 0;
    h->mb.i_mb_top_xy = -1;
    h->mb.i_mb_top_y = -1;
    h->mb.i_mb_left_xy[0] = h->mb.i_mb_left_xy[1] = -1;
    h->mb.i_mb_topleft_xy = -1;
    h->mb.i_mb_topright_xy = -1;
    h->mb.i_mb_type_top = -1;
    h->mb.i_mb_type_left[0] = h->mb.i_mb_type_left[1] = -1;
    h->mb.i_mb_type_topleft = -1;
    h->mb.i_mb_type_topright = -1;
    h->mb.left_index_table = &left_indices[3];
    h->mb.topleft_partition = 0;

    int topleft_y = top_y;
    int topright_y = top_y;
    int left[2];

    left[0] = left[1] = h->mb.i_mb_xy - 1;
    h->mb.left_b8[0] = h->mb.left_b8[1] = h->mb.i_b8_xy - 2;
    h->mb.left_b4[0] = h->mb.left_b4[1] = h->mb.i_b4_xy - 4;

    if( mb_x > 0 )
    {
        h->mb.i_neighbour_frame |= MB_LEFT;
        h->mb.i_mb_left_xy[0] = left[0];
        h->mb.i_mb_left_xy[1] = left[1];
        h->mb.i_mb_type_left[0] = h->mb.type[h->mb.i_mb_left_xy[0]];
        h->mb.i_mb_type_left[1] = h->mb.type[h->mb.i_mb_left_xy[1]];
        if( h->mb.slice_table[left[0]] == h->sh.i_first_mb )
        {
            h->mb.i_neighbour |= MB_LEFT;

            // FIXME: We don't currently support constrained intra + mbaff.
            if( !h->param.b_constrained_intra || IS_INTRA( h->mb.i_mb_type_left[0] ) )
                h->mb.i_neighbour_intra |= MB_LEFT;
        }
    }

    /* We can't predict from the previous threadslice since it hasn't been encoded yet. */
    if( (h->i_threadslice_start) != (mb_y) )
    {
        if( top >= 0 )
        {
            h->mb.i_neighbour_frame |= MB_TOP;
            h->mb.i_mb_top_xy = top;
            h->mb.i_mb_top_y = top_y;
            h->mb.i_mb_type_top = h->mb.type[h->mb.i_mb_top_xy];
            if( h->mb.slice_table[top] == h->sh.i_first_mb )
            {
                h->mb.i_neighbour |= MB_TOP;

                if( !h->param.b_constrained_intra || IS_INTRA( h->mb.i_mb_type_top ) )
                    h->mb.i_neighbour_intra |= MB_TOP;
            }
        }

        if( mb_x > 0 && topleft_y >= 0  )
        {
            h->mb.i_neighbour_frame |= MB_TOPLEFT;
            h->mb.i_mb_topleft_xy = h->mb.i_mb_stride*topleft_y + mb_x - 1;
            h->mb.i_mb_topleft_y = topleft_y;
            h->mb.i_mb_type_topleft = h->mb.type[h->mb.i_mb_topleft_xy];
            if( h->mb.slice_table[h->mb.i_mb_topleft_xy] == h->sh.i_first_mb )
            {
                h->mb.i_neighbour |= MB_TOPLEFT;

                if( !h->param.b_constrained_intra || IS_INTRA( h->mb.i_mb_type_topleft ) )
                    h->mb.i_neighbour_intra |= MB_TOPLEFT;
            }
        }

        if( mb_x < h->mb.i_mb_width - 1 && topright_y >= 0 )
        {
            h->mb.i_neighbour_frame |= MB_TOPRIGHT;
            h->mb.i_mb_topright_xy = h->mb.i_mb_stride*topright_y + mb_x + 1;
            h->mb.i_mb_topright_y = topright_y;
            h->mb.i_mb_type_topright = h->mb.type[h->mb.i_mb_topright_xy];
            if( h->mb.slice_table[h->mb.i_mb_topright_xy] == h->sh.i_first_mb )
            {
                h->mb.i_neighbour |= MB_TOPRIGHT;

                if( !h->param.b_constrained_intra || IS_INTRA( h->mb.i_mb_type_topright ) )
                    h->mb.i_neighbour_intra |= MB_TOPRIGHT;
            }
        }
    }
}

#define LTOP 0
#define LBOT 0

static void  x264_macroblock_cache_load( x264_t *h, int mb_x, int mb_y)
{
    x264_macroblock_cache_load_neighbours( h, mb_x, mb_y );

    int *left = h->mb.i_mb_left_xy;
    int top  = h->mb.i_mb_top_xy;
    int top_y = h->mb.i_mb_top_y;
    int s8x8 = h->mb.i_b8_stride;
    int s4x4 = h->mb.i_b4_stride;
    int top_8x8 = (2*top_y+1) * s8x8 + 2*mb_x;
    int top_4x4 = (4*top_y+3) * s4x4 + 4*mb_x;

    /* GCC pessimizes direct loads from heap-allocated arrays due to aliasing. */
    /* By only dereferencing them once, we avoid this issue. */
    int8_t (*i4x4)[8] = h->mb.intra4x4_pred_mode;
    uint8_t (*nnz)[48] = h->mb.non_zero_count;
    int16_t *cbp = h->mb.cbp;

    const x264_left_table_t *left_index_table = h->mb.left_index_table;

    h->mb.cache.deblock_strength = h->deblock_strength[mb_y&1][mb_x];

    /* load cache */
    if( h->mb.i_neighbour & MB_TOP )
    {
        h->mb.cache.i_cbp_top = cbp[top];
        /* load intra4x4 */
        CP32( &h->mb.cache.intra4x4_pred_mode[x264_scan8[0] - 8], &i4x4[top][0] );

        /* load non_zero_count */
        CP32( &h->mb.cache.non_zero_count[x264_scan8[ 0] - 8], &nnz[top][12] );
        CP32( &h->mb.cache.non_zero_count[x264_scan8[16] - 8], &nnz[top][16-4 + (16>>CHROMA_V_SHIFT)] );
        CP32( &h->mb.cache.non_zero_count[x264_scan8[32] - 8], &nnz[top][32-4 + (16>>CHROMA_V_SHIFT)] );
    }
    else
    {
        h->mb.cache.i_cbp_top = -1;

        /* load intra4x4 */
        M32( &h->mb.cache.intra4x4_pred_mode[x264_scan8[0] - 8] ) = 0xFFFFFFFFU;

        /* load non_zero_count */
        M32( &h->mb.cache.non_zero_count[x264_scan8[ 0] - 8] ) = 0x80808080U;
        M32( &h->mb.cache.non_zero_count[x264_scan8[16] - 8] ) = 0x80808080U;
        M32( &h->mb.cache.non_zero_count[x264_scan8[32] - 8] ) = 0x80808080U;
    }

    if( h->mb.i_neighbour & MB_LEFT )
    {
        int ltop = left[LTOP];
        int lbot = ltop;
        h->mb.cache.i_cbp_left = cbp[ltop];

        /* load intra4x4 */
        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 0] - 1] = i4x4[ltop][left_index_table->intra[0]];
        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 2] - 1] = i4x4[ltop][left_index_table->intra[1]];
        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 8] - 1] = i4x4[lbot][left_index_table->intra[2]];
        h->mb.cache.intra4x4_pred_mode[x264_scan8[10] - 1] = i4x4[lbot][left_index_table->intra[3]];

        /* load non_zero_count */
        h->mb.cache.non_zero_count[x264_scan8[ 0] - 1] = nnz[ltop][left_index_table->nnz[0]];
        h->mb.cache.non_zero_count[x264_scan8[ 2] - 1] = nnz[ltop][left_index_table->nnz[1]];
        h->mb.cache.non_zero_count[x264_scan8[ 8] - 1] = nnz[lbot][left_index_table->nnz[2]];
        h->mb.cache.non_zero_count[x264_scan8[10] - 1] = nnz[lbot][left_index_table->nnz[3]];

        {
            h->mb.cache.non_zero_count[x264_scan8[16+ 0] - 1] = nnz[ltop][left_index_table->nnz_chroma[0]];
            h->mb.cache.non_zero_count[x264_scan8[16+ 2] - 1] = nnz[lbot][left_index_table->nnz_chroma[1]];
            h->mb.cache.non_zero_count[x264_scan8[32+ 0] - 1] = nnz[ltop][left_index_table->nnz_chroma[2]];
            h->mb.cache.non_zero_count[x264_scan8[32+ 2] - 1] = nnz[lbot][left_index_table->nnz_chroma[3]];
        }
    }
    else
    {
        h->mb.cache.i_cbp_left = -1;

        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 0] - 1] =
        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 2] - 1] =
        h->mb.cache.intra4x4_pred_mode[x264_scan8[ 8] - 1] =
        h->mb.cache.intra4x4_pred_mode[x264_scan8[10] - 1] = -1;

        /* load non_zero_count */
        h->mb.cache.non_zero_count[x264_scan8[ 0] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[ 2] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[ 8] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[10] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[16+ 0] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[16+ 2] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[32+ 0] - 1] =
        h->mb.cache.non_zero_count[x264_scan8[32+ 2] - 1] = 0x80;
    }

    {
        x264_copy_column8( h->mb.pic.p_fdec[0]-1+ 4*FDEC_STRIDE, h->mb.pic.p_fdec[0]+15+ 4*FDEC_STRIDE );
        x264_copy_column8( h->mb.pic.p_fdec[0]-1+12*FDEC_STRIDE, h->mb.pic.p_fdec[0]+15+12*FDEC_STRIDE );
        x264_macroblock_load_pic_pointers( h, mb_x, mb_y, 0, 0 );
        {
            x264_copy_column8( h->mb.pic.p_fdec[1]-1+ 4*FDEC_STRIDE, h->mb.pic.p_fdec[1]+ 7+ 4*FDEC_STRIDE );
            x264_copy_column8( h->mb.pic.p_fdec[2]-1+ 4*FDEC_STRIDE, h->mb.pic.p_fdec[2]+ 7+ 4*FDEC_STRIDE );
            x264_macroblock_load_pic_pointers( h, mb_x, mb_y, 1, 1 );
        }
    }

    if( h->fdec->integral )
    {
        int offset = 16 * (mb_x + mb_y * h->fdec->i_stride[0]);
        for( int i = 0; i < h->mb.pic.i_fref; i++ )
            h->mb.pic.p_integral[i] = &h->fref[i]->integral[offset];
    }

    /* load ref/mv/mvd */
    {
        int16_t (*mv)[2] = h->mb.mv;
        int8_t *ref = h->mb.ref;

        int i8 = x264_scan8[0] - 1 - 1*8;
        if( h->mb.i_neighbour & MB_TOPLEFT )
        {
            int ir = top_8x8 - 1;
            int iv = top_4x4 - 1;
            h->mb.cache.ref[i8] = ref[ir];
            CP32( h->mb.cache.mv[i8], mv[iv] );
        }
        else
        {
            h->mb.cache.ref[i8] = -2;
            M32( h->mb.cache.mv[i8] ) = 0;
        }

        i8 = x264_scan8[0] - 8;
        if( h->mb.i_neighbour & MB_TOP )
        {
            h->mb.cache.ref[i8+0] =
            h->mb.cache.ref[i8+1] = ref[top_8x8 + 0];
            h->mb.cache.ref[i8+2] =
            h->mb.cache.ref[i8+3] = ref[top_8x8 + 1];
            CP128( h->mb.cache.mv[i8], mv[top_4x4] );
        }
        else
        {
            M128( h->mb.cache.mv[i8] ) = M128_ZERO;
            M32( &h->mb.cache.ref[i8] ) = (uint8_t)(-2) * 0x01010101U;
        }

        i8 = x264_scan8[0] + 4 - 1*8;
        if( h->mb.i_neighbour & MB_TOPRIGHT )
        {
            int ir = top_8x8 + 2;
            int iv = top_4x4 + 4;
            h->mb.cache.ref[i8] = ref[ir];
            CP32( h->mb.cache.mv[i8], mv[iv] );
        }
        else
             h->mb.cache.ref[i8] = -2;

        i8 = x264_scan8[0] - 1;
        if( h->mb.i_neighbour & MB_LEFT )
        {
            {
                const int ir = h->mb.i_b8_xy - 1;
                const int iv = h->mb.i_b4_xy - 1;
                h->mb.cache.ref[i8+0*8] =
                h->mb.cache.ref[i8+1*8] = ref[ir + 0*s8x8];
                h->mb.cache.ref[i8+2*8] =
                h->mb.cache.ref[i8+3*8] = ref[ir + 1*s8x8];

                CP32( h->mb.cache.mv[i8+0*8], mv[iv + 0*s4x4] );
                CP32( h->mb.cache.mv[i8+1*8], mv[iv + 1*s4x4] );
                CP32( h->mb.cache.mv[i8+2*8], mv[iv + 2*s4x4] );
                CP32( h->mb.cache.mv[i8+3*8], mv[iv + 3*s4x4] );
            }
        }
        else
        {
            for( int i = 0; i < 4; i++ )
            {
                h->mb.cache.ref[i8+i*8] = -2;
                M32( h->mb.cache.mv[i8+i*8] ) = 0;
            }
        }
    }

    /* Check whether skip here would cause decoder to predict interlace mode incorrectly.
     * FIXME: It might be better to change the interlace type rather than forcing a skip to be non-skip. */
    h->mb.b_allow_skip = 1;

    if( h->sh.i_type == SLICE_TYPE_P )
        x264_mb_predict_mv_pskip( h, h->mb.cache.pskip_mv );

    h->mb.i_neighbour4[0] = (h->mb.i_neighbour_intra & (MB_TOP|MB_LEFT|MB_TOPLEFT))
                            | ((h->mb.i_neighbour_intra & MB_TOP) ? MB_TOPRIGHT : 0);
    h->mb.i_neighbour4[4] =
    h->mb.i_neighbour4[1] = MB_LEFT | ((h->mb.i_neighbour_intra & MB_TOP) ? (MB_TOP|MB_TOPLEFT|MB_TOPRIGHT) : 0);
    h->mb.i_neighbour4[2] =
    h->mb.i_neighbour4[8] =
    h->mb.i_neighbour4[10] = MB_TOP|MB_TOPRIGHT | ((h->mb.i_neighbour_intra & MB_LEFT) ? (MB_LEFT|MB_TOPLEFT) : 0);
    h->mb.i_neighbour4[5] = MB_LEFT | (h->mb.i_neighbour_intra & MB_TOPRIGHT)
                            | ((h->mb.i_neighbour_intra & MB_TOP) ? MB_TOP|MB_TOPLEFT : 0);
}

void x264_macroblock_cache_load_progressive( x264_t *h, int mb_x, int mb_y )
{
    x264_macroblock_cache_load( h, mb_x, mb_y );
}

void x264_macroblock_deblock_strength( x264_t *h )
{
    uint8_t (*bs)[8][4] = h->mb.cache.deblock_strength;
    if( IS_INTRA( h->mb.i_type ) )
    {
        memset( bs[0][1], 3, 3*4*sizeof(uint8_t) );
        memset( bs[1][1], 3, 3*4*sizeof(uint8_t) );
        return;
    }

    int neighbour_changed = 0;
    if( h->sh.i_disable_deblocking_filter_idc != 2 )
    {
        neighbour_changed = h->mb.i_neighbour_frame&~h->mb.i_neighbour;
        h->mb.i_neighbour = h->mb.i_neighbour_frame;
    }

    /* If we have multiple slices and we're deblocking on slice edges, we
     * have to reload neighbour data. */
    if( neighbour_changed )
    {
        int top_y = h->mb.i_mb_top_y;
        int top_8x8 = (2*top_y+1) * h->mb.i_b8_stride + 2*h->mb.i_mb_x;
        int top_4x4 = (4*top_y+3) * h->mb.i_b4_stride + 4*h->mb.i_mb_x;
        int s8x8 = h->mb.i_b8_stride;
        int s4x4 = h->mb.i_b4_stride;

        uint8_t (*nnz)[48] = h->mb.non_zero_count;
        const x264_left_table_t *left_index_table = &left_indices[3];

        if( neighbour_changed & MB_TOP )
            CP32( &h->mb.cache.non_zero_count[x264_scan8[0] - 8], &nnz[h->mb.i_mb_top_xy][12] );

        if( neighbour_changed & MB_LEFT )
        {
            int *left = h->mb.i_mb_left_xy;
            h->mb.cache.non_zero_count[x264_scan8[0 ] - 1] = nnz[left[0]][left_index_table->nnz[0]];
            h->mb.cache.non_zero_count[x264_scan8[2 ] - 1] = nnz[left[0]][left_index_table->nnz[1]];
            h->mb.cache.non_zero_count[x264_scan8[8 ] - 1] = nnz[left[1]][left_index_table->nnz[2]];
            h->mb.cache.non_zero_count[x264_scan8[10] - 1] = nnz[left[1]][left_index_table->nnz[3]];
        }

        {
            int16_t (*mv)[2] = h->mb.mv;
            int8_t *ref = h->mb.ref;

            int i8 = x264_scan8[0] - 8;
            if( neighbour_changed & MB_TOP )
            {
                h->mb.cache.ref[i8+0] =
                h->mb.cache.ref[i8+1] = ref[top_8x8 + 0];
                h->mb.cache.ref[i8+2] =
                h->mb.cache.ref[i8+3] = ref[top_8x8 + 1];
                CP128( h->mb.cache.mv[i8], mv[top_4x4] );
            }

            i8 = x264_scan8[0] - 1;
            if( neighbour_changed & MB_LEFT )
            {
                h->mb.cache.ref[i8+0*8] =
                h->mb.cache.ref[i8+1*8] = ref[h->mb.left_b8[0] + 1 + s8x8*left_index_table->ref[0]];
                h->mb.cache.ref[i8+2*8] =
                h->mb.cache.ref[i8+3*8] = ref[h->mb.left_b8[1] + 1 + s8x8*left_index_table->ref[2]];

                CP32( h->mb.cache.mv[i8+0*8], mv[h->mb.left_b4[0] + 3 + s4x4*left_index_table->mv[0]] );
                CP32( h->mb.cache.mv[i8+1*8], mv[h->mb.left_b4[0] + 3 + s4x4*left_index_table->mv[1]] );
                CP32( h->mb.cache.mv[i8+2*8], mv[h->mb.left_b4[1] + 3 + s4x4*left_index_table->mv[2]] );
                CP32( h->mb.cache.mv[i8+3*8], mv[h->mb.left_b4[1] + 3 + s4x4*left_index_table->mv[3]] );
            }
        }
    }

    h->loopf.deblock_strength( h->mb.cache.non_zero_count, h->mb.cache.ref, h->mb.cache.mv,
                               bs, 4 );
}

static void  x264_macroblock_store_pic( x264_t *h, int mb_x, int mb_y, int i, int b_chroma )
{
    int height = b_chroma ? 16>>CHROMA_V_SHIFT : 16;
    int i_stride = h->fdec->i_stride[i];
    int i_stride2 = i_stride ;
    int i_pix_offset = 16 * mb_x + height * mb_y * i_stride;
    if( b_chroma )
        h->mc.store_interleave_chroma( &h->fdec->plane[1][i_pix_offset], i_stride2, h->mb.pic.p_fdec[1], h->mb.pic.p_fdec[2], height );
    else
        h->mc.copy[PIXEL_16x16]( &h->fdec->plane[i][i_pix_offset], i_stride2, h->mb.pic.p_fdec[i], FDEC_STRIDE, 16 );
}

static void  x264_macroblock_backup_intra( x264_t *h, int mb_x, int mb_y)
{
    /* In MBAFF we store the last two rows in intra_border_backup[0] and [1].
     * For progressive mbs this is the bottom two rows, and for interlaced the
     * bottom row of each field. We also store samples needed for the next
     * mbpair in intra_border_backup[2]. */
    int backup_dst =  (mb_y&1) ;
    memcpy( &h->intra_border_backup[backup_dst][0][mb_x*16  ], h->mb.pic.p_fdec[0]+FDEC_STRIDE*15, 16*sizeof(pixel) );

    {
        int backup_src = (15>>CHROMA_V_SHIFT) * FDEC_STRIDE;
        memcpy( &h->intra_border_backup[backup_dst][1][mb_x*16  ], h->mb.pic.p_fdec[1]+backup_src, 8*sizeof(pixel) );
        memcpy( &h->intra_border_backup[backup_dst][1][mb_x*16+8], h->mb.pic.p_fdec[2]+backup_src, 8*sizeof(pixel) );
    }
}

void x264_macroblock_cache_save( x264_t *h )
{
    const int i_mb_xy = h->mb.i_mb_xy;
    const int i_mb_type = x264_mb_type_fix[h->mb.i_type];
    const int s8x8 = h->mb.i_b8_stride;
    const int s4x4 = h->mb.i_b4_stride;
    const int i_mb_4x4 = h->mb.i_b4_xy;
    const int i_mb_8x8 = h->mb.i_b8_xy;

    /* GCC pessimizes direct stores to heap-allocated arrays due to aliasing. */
    /* By only dereferencing them once, we avoid this issue. */
    int8_t *i4x4 = h->mb.intra4x4_pred_mode[i_mb_xy];
    uint8_t *nnz = h->mb.non_zero_count[i_mb_xy];

    {
        x264_macroblock_backup_intra( h, h->mb.i_mb_x, h->mb.i_mb_y );
        x264_macroblock_store_pic( h, h->mb.i_mb_x, h->mb.i_mb_y, 0, 0 );
        x264_macroblock_store_pic( h, h->mb.i_mb_x, h->mb.i_mb_y, 1, 1 );
    }

    h->mb.type[i_mb_xy] = i_mb_type;
    h->mb.slice_table[i_mb_xy] = h->sh.i_first_mb;
    h->mb.partition[i_mb_xy] = IS_INTRA( i_mb_type ) ? D_16x16 : h->mb.i_partition;
    h->mb.i_mb_prev_xy = i_mb_xy;

    /* save intra4x4 */
    if( i_mb_type == I_4x4 )
    {
        CP32( &i4x4[0], &h->mb.cache.intra4x4_pred_mode[x264_scan8[10]] );
        M32( &i4x4[4] ) = pack8to32( h->mb.cache.intra4x4_pred_mode[x264_scan8[5] ],
                                     h->mb.cache.intra4x4_pred_mode[x264_scan8[7] ],
                                     h->mb.cache.intra4x4_pred_mode[x264_scan8[13] ], 0);
    }
    else if( !h->param.b_constrained_intra || IS_INTRA(i_mb_type) )
        M64( i4x4 ) = I_PRED_4x4_DC * 0x0101010101010101ULL;
    else
        M64( i4x4 ) = (uint8_t)(-1) * 0x0101010101010101ULL;


    if( i_mb_type == I_PCM )
    {
        h->mb.qp[i_mb_xy] = 0;
        h->mb.i_last_dqp = 0;
        h->mb.i_cbp_chroma = 2;
        h->mb.i_cbp_luma = 0xf;
        h->mb.cbp[i_mb_xy] = (h->mb.i_cbp_chroma << 4) | h->mb.i_cbp_luma | 0x700;
        for( int i = 0; i < 48; i++ )
            h->mb.cache.non_zero_count[x264_scan8[i]] = 16;
    }
    else
    {
        if( h->mb.i_type != I_16x16 && h->mb.i_cbp_luma == 0 && h->mb.i_cbp_chroma == 0 )
            h->mb.i_qp = h->mb.i_last_qp;
        h->mb.qp[i_mb_xy] = h->mb.i_qp;
        h->mb.i_last_dqp = h->mb.i_qp - h->mb.i_last_qp;
        h->mb.i_last_qp = h->mb.i_qp;
    }

    /* save non zero count */
    CP32( &nnz[ 0+0*4], &h->mb.cache.non_zero_count[x264_scan8[ 0]] );
    CP32( &nnz[ 0+1*4], &h->mb.cache.non_zero_count[x264_scan8[ 2]] );
    CP32( &nnz[ 0+2*4], &h->mb.cache.non_zero_count[x264_scan8[ 8]] );
    CP32( &nnz[ 0+3*4], &h->mb.cache.non_zero_count[x264_scan8[10]] );
    CP32( &nnz[16+0*4], &h->mb.cache.non_zero_count[x264_scan8[16+0]] );
    CP32( &nnz[16+1*4], &h->mb.cache.non_zero_count[x264_scan8[16+2]] );
    CP32( &nnz[32+0*4], &h->mb.cache.non_zero_count[x264_scan8[32+0]] );
    CP32( &nnz[32+1*4], &h->mb.cache.non_zero_count[x264_scan8[32+2]] );

    if( h->sh.i_type != SLICE_TYPE_I )
    {
        int16_t (*mv0)[2] = &h->mb.mv[i_mb_4x4];
        int8_t *ref0 = &h->mb.ref[i_mb_8x8];
        if( !IS_INTRA( i_mb_type ) )
        {
            ref0[0+0*s8x8] = h->mb.cache.ref[x264_scan8[0]];
            ref0[1+0*s8x8] = h->mb.cache.ref[x264_scan8[4]];
            ref0[0+1*s8x8] = h->mb.cache.ref[x264_scan8[8]];
            ref0[1+1*s8x8] = h->mb.cache.ref[x264_scan8[12]];
            CP128( &mv0[0*s4x4], h->mb.cache.mv[x264_scan8[0]+8*0] );
            CP128( &mv0[1*s4x4], h->mb.cache.mv[x264_scan8[0]+8*1] );
            CP128( &mv0[2*s4x4], h->mb.cache.mv[x264_scan8[0]+8*2] );
            CP128( &mv0[3*s4x4], h->mb.cache.mv[x264_scan8[0]+8*3] );
        }
        else
        {
            M16( &ref0[0*s8x8] ) = (uint8_t)(-1) * 0x0101;
            M16( &ref0[1*s8x8] ) = (uint8_t)(-1) * 0x0101;
            M128( &mv0[0*s4x4] ) = M128_ZERO;
            M128( &mv0[1*s4x4] ) = M128_ZERO;
            M128( &mv0[2*s4x4] ) = M128_ZERO;
            M128( &mv0[3*s4x4] ) = M128_ZERO;
        }
    }

}

