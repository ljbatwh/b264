/*****************************************************************************
 * frame.h: frame handling
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

#ifndef X264_FRAME_H
#define X264_FRAME_H

/* number of pixels past the edge of the frame, for motion estimation/compensation */
#define PADH 32
#define PADV 32

typedef struct x264_frame
{
    /* */
    uint8_t *base;       /* Base pointer for all malloced data in this frame. */
    int     i_poc;
    int     i_delta_poc[2];
    int     i_type;
    int     i_qpplus1;
    int64_t i_duration;  /* in SPS time_scale units (i.e 2 * timebase units) used for vfr */
    float   f_duration;  /* in seconds */
    int64_t i_cpb_duration;
    int64_t i_cpb_delay; /* in SPS time_scale units (i.e 2 * timebase units) */
    int64_t i_dpb_output_delay;
    x264_param_t *param;

    int     i_frame;     /* Presentation frame number */
    int     i_coded;     /* Coded frame number */
    int64_t i_field_cnt; /* Presentation field count */
    int     i_frame_num; /* 7.4.3 frame_num */
    int     b_kept_as_ref;
    int     b_keyframe;
    uint8_t b_fdec;
    uint8_t b_last_minigop_bframe; /* this frame is the last b in a sequence of bframes */
    float   f_qp_avg_rc; /* QPs as decided by ratecontrol */
    float   f_qp_avg_aq; /* QPs as decided by AQ in addition to ratecontrol */
    float   f_crf_avg;   /* Average effective CRF for this frame */
    int     i_poc_l0ref0; /* poc of first refframe in L0, used to check if direct temporal is possible */

    /* YUV buffer */
    int     i_csp; /* Internal csp */
    int     i_plane;
    int     i_stride[3];
    int     i_width[3];
    int     i_lines[3];
    int     i_stride_lowres;
    int     i_width_lowres;
    int     i_lines_lowres;
    pixel *plane[3];
    pixel *plane_fld[3];
    pixel *filtered[3][4]; /* plane[0], H, V, HV */
    pixel *filtered_fld[3][4];
    pixel *lowres[4]; /* half-size copy of input frame: Orig, H, V, HV */
    uint16_t *integral;

    /* for unrestricted mv we allocate more data than needed
     * allocated data are stored in buffer */
    pixel *buffer[4];
    pixel *buffer_fld[4];
    pixel *buffer_lowres[4];

    int b_duplicate;
    struct x264_frame *orig;

    /* motion data */
    int8_t  *mb_type;
    uint8_t *mb_partition;
    int16_t (*mv)[2]; //just for p ref. [2] x,y
    int16_t (*mv16x16)[2];
    int16_t (*lowres_mvs[2][X264_BFRAME_MAX+1])[2];
    uint8_t *field;
    uint8_t *effective_qp;

    /* Stored as (lists_used << LOWRES_COST_SHIFT) + (cost).
     * Doesn't need special addressing for intra cost because
     * lists_used is guaranteed to be zero in that cast. */
    uint16_t (*lowres_costs[X264_BFRAME_MAX+2][X264_BFRAME_MAX+2]);
    #define LOWRES_COST_MASK ((1<<14)-1)
    #define LOWRES_COST_SHIFT 14

    int     *lowres_mv_costs[2][X264_BFRAME_MAX+1];
    int8_t  *ref;
    int     i_ref;
    int     ref_poc[X264_REF_MAX];
    int16_t inv_ref_poc; // inverse values of ref0 poc to avoid divisions in temporal MV prediction

    /* for adaptive B-frame decision.
     * contains the SATD cost of the lowres frame encoded in various modes
     * FIXME: how big an array do we need? */
    int     i_cost_est[2]; //for intra and inter
    int     i_cost_est_aq[2];
    int     i_satd; // the i_cost_est of the selected frametype
    int     i_intra_mbs[2];
    int     *i_row_satds[2]; //every elements is a array, length is mbs
    int     *i_row_satd;
    int     *i_row_bits;
    float   *f_row_qp;
    float   *f_row_qscale;
    float   *f_qp_offset;
    float   *f_qp_offset_aq;
    int     b_intra_calculated;
    uint16_t *i_intra_cost;
    uint16_t *i_propagate_cost;
    uint16_t *i_inv_qscale_factor;
    int     b_scenecut; /* Set to zero if the frame cannot possibly be part of a real scenecut. */
    uint32_t i_pixel_sum[3];
    uint64_t i_pixel_ssd[3];

    /* hrd */
    x264_hrd_t hrd_timing;

    /* vbv */
    double f_planned_cpb_duration;

    /* threading */
    int     i_lines_completed; /* in pixels */
    int     i_reference_count; /* number of threads using this frame (not necessarily the number of pointers) */
    int     i_slice_count; /* Atomically written to/read from with slice threads */

    /* periodic intra refresh */
    float   f_pir_position;
    int     i_pir_start_col;
    int     i_pir_end_col;
    int     i_frames_since_pir;

    /* user sei */
    x264_sei_t extra_sei;

    /* user data */
    void *opaque;

} x264_frame_t;

/* synchronized frame list */
typedef struct
{
   x264_frame_t **list;
   int i_max_size;
   int i_size;
} x264_sync_frame_list_t;

typedef void (*x264_deblock_inter_t)( pixel *pix, intptr_t stride, int alpha, int beta, int8_t *tc0 );
typedef void (*x264_deblock_intra_t)( pixel *pix, intptr_t stride, int alpha, int beta );
typedef struct
{
    x264_deblock_inter_t deblock_luma[2];
    x264_deblock_inter_t deblock_chroma[2];
    x264_deblock_inter_t deblock_h_chroma_420;
    x264_deblock_intra_t deblock_luma_intra[2];
    x264_deblock_intra_t deblock_chroma_intra[2];
    x264_deblock_intra_t deblock_h_chroma_420_intra;
    x264_deblock_inter_t deblock_luma_mbaff;
    x264_deblock_inter_t deblock_chroma_mbaff;
    x264_deblock_inter_t deblock_chroma_420_mbaff;
    x264_deblock_intra_t deblock_luma_intra_mbaff;
    x264_deblock_intra_t deblock_chroma_intra_mbaff;
    x264_deblock_intra_t deblock_chroma_420_intra_mbaff;
    void (*deblock_strength) ( uint8_t nnz[X264_SCAN8_SIZE], int8_t ref[X264_SCAN8_LUMA_SIZE],
                               int16_t mv[X264_SCAN8_LUMA_SIZE][2], uint8_t bs[8][4], int mvy_limit);
} x264_deblock_function_t;

void          x264_frame_delete( x264_frame_t *frame );

int           x264_frame_copy_picture( x264_t *h, x264_frame_t *dst, x264_picture_t *src );

void          x264_frame_expand_border( x264_t *h, x264_frame_t *frame, int mb_y );
void          x264_frame_expand_border_filtered( x264_t *h, x264_frame_t *frame, int mb_y, int b_end );
void          x264_frame_expand_border_lowres( x264_frame_t *frame );
void          x264_frame_expand_border_chroma( x264_t *h, x264_frame_t *frame, int plane );
void          x264_frame_expand_border_mod16( x264_t *h, x264_frame_t *frame );

void          x264_frame_deblock_row( x264_t *h, int mb_y );
void          x264_macroblock_deblock( x264_t *h );

void          x264_frame_filter( x264_t *h, x264_frame_t *frame, int mb_y, int b_end );
void          x264_frame_init_lowres( x264_t *h, x264_frame_t *frame );

void          x264_deblock_init( x264_deblock_function_t *pf);

void          x264_frame_cond_wait( x264_frame_t *frame, int i_lines_completed );
int           x264_frame_new_slice( x264_t *h, x264_frame_t *frame );

void          x264_threadslice_cond_broadcast( x264_t *h, int pass );
void          x264_threadslice_cond_wait( x264_t *h, int pass );

void          x264_frame_push( x264_frame_t **list, x264_frame_t *frame );
x264_frame_t *x264_frame_pop( x264_frame_t **list );
void          x264_frame_unshift( x264_frame_t **list, x264_frame_t *frame );
x264_frame_t *x264_frame_shift( x264_frame_t **list );
void          x264_frame_push_unused( x264_t *h, x264_frame_t *frame );
void          x264_frame_push_blank_unused( x264_t *h, x264_frame_t *frame );
x264_frame_t *x264_frame_pop_blank_unused( x264_t *h );
x264_frame_t *x264_frame_pop_unused( x264_t *h, int b_fdec );
void          x264_frame_delete_list( x264_frame_t **list );


#endif
