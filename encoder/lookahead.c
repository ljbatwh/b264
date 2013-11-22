/*****************************************************************************
 * lookahead.c: high-level lookahead functions
 *****************************************************************************
 * Copyright (C) 2010-2013 Avail Media and x264 project
 *
 * Authors: Michael Kazmier <mkazmier@availmedia.com>
 *          Alex Giladi <agiladi@availmedia.com>
 *          Steven Walters <kemuri9@gmail.com>
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

/* LOOKAHEAD (threaded and non-threaded mode)
 *
 * Lookahead types:
 *     [1] Slice type / scene cut;
 *
 * In non-threaded mode, we run the existing slicetype decision code as it was.
 * In threaded mode, we run in a separate thread, that lives between the calls
 * to x264_encoder_open() and x264_encoder_close(), and performs lookahead for
 * the number of frames specified in rc_lookahead.  Recommended setting is
 * # of bframes + # of threads.
 */
#include "common/common.h"
#include "analyse.h"

static void x264_lookahead_update_last_nonb( x264_t *h, x264_frame_t *new_nonb )
{
    if( h->lookahead->last_nonb )
        x264_frame_push_unused( h, h->lookahead->last_nonb );
    h->lookahead->last_nonb = new_nonb;
    new_nonb->i_reference_count++;
}

int x264_lookahead_init( x264_t *h )
{
    x264_lookahead_t *look;
    CHECKED_MALLOCZERO( look, sizeof(x264_lookahead_t) );
    h->lookahead = look;

    look->i_last_keyframe = - h->param.i_keyint_max;

    return 0;

fail:
    x264_free( look );
    return -1;
}

void x264_lookahead_delete( x264_t *h )
{
    if( h->lookahead->last_nonb )
        x264_frame_push_unused( h, h->lookahead->last_nonb );
}

void x264_lookahead_get_frame( x264_t *h )
{
    if( h->frames.current || !h->lookahead->current_frame )
        return ;

    x264_slicetype_decide( h );

    if( h->lookahead->last_nonb )
        x264_frame_push_unused( h, h->lookahead->last_nonb );
    h->lookahead->last_nonb = h->lookahead->current_frame;
    h->lookahead->current_frame->i_reference_count++;

    h->frames.current = h->lookahead->current_frame;
    h->lookahead->current_frame = 0;
}
