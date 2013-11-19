/*****************************************************************************
 * predict.c: intra prediction
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

/* predict4x4 are inspired from ffmpeg h264 decoder */


#include "common.h"


/****************************************************************************
 * 16x16 prediction for intra luma block
 ****************************************************************************/

#define PREDICT_16x16_DC(v)\
    for( int i = 0; i < 16; i++ )\
    {\
        MPIXEL_X4( src+ 0 ) = v;\
        MPIXEL_X4( src+ 4 ) = v;\
        MPIXEL_X4( src+ 8 ) = v;\
        MPIXEL_X4( src+12 ) = v;\
        src += FDEC_STRIDE;\
    }

void x264_predict_16x16_dc_c( pixel *src )
{
    int dc = 0;

    for( int i = 0; i < 16; i++ )
    {
        dc += src[-1 + i * FDEC_STRIDE];
        dc += src[i - FDEC_STRIDE];
    }
    pixel4 dcsplat = PIXEL_SPLAT_X4( ( dc + 16 ) >> 5 );

    PREDICT_16x16_DC( dcsplat );
}
static void x264_predict_16x16_dc_left_c( pixel *src )
{
    int dc = 0;

    for( int i = 0; i < 16; i++ )
        dc += src[-1 + i * FDEC_STRIDE];
    pixel4 dcsplat = PIXEL_SPLAT_X4( ( dc + 8 ) >> 4 );

    PREDICT_16x16_DC( dcsplat );
}
static void x264_predict_16x16_dc_top_c( pixel *src )
{
    int dc = 0;

    for( int i = 0; i < 16; i++ )
        dc += src[i - FDEC_STRIDE];
    pixel4 dcsplat = PIXEL_SPLAT_X4( ( dc + 8 ) >> 4 );

    PREDICT_16x16_DC( dcsplat );
}
static void x264_predict_16x16_dc_128_c( pixel *src )
{
    PREDICT_16x16_DC( PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) ) );
}
void x264_predict_16x16_h_c( pixel *src )
{
    for( int i = 0; i < 16; i++ )
    {
        const pixel4 v = PIXEL_SPLAT_X4( src[-1] );
        MPIXEL_X4( src+ 0 ) = v;
        MPIXEL_X4( src+ 4 ) = v;
        MPIXEL_X4( src+ 8 ) = v;
        MPIXEL_X4( src+12 ) = v;
        src += FDEC_STRIDE;
    }
}
void x264_predict_16x16_v_c( pixel *src )
{
    pixel4 v0 = MPIXEL_X4( &src[ 0-FDEC_STRIDE] );
    pixel4 v1 = MPIXEL_X4( &src[ 4-FDEC_STRIDE] );
    pixel4 v2 = MPIXEL_X4( &src[ 8-FDEC_STRIDE] );
    pixel4 v3 = MPIXEL_X4( &src[12-FDEC_STRIDE] );

    for( int i = 0; i < 16; i++ )
    {
        MPIXEL_X4( src+ 0 ) = v0;
        MPIXEL_X4( src+ 4 ) = v1;
        MPIXEL_X4( src+ 8 ) = v2;
        MPIXEL_X4( src+12 ) = v3;
        src += FDEC_STRIDE;
    }
}
void x264_predict_16x16_p_c( pixel *src )
{
    int H = 0, V = 0;

    /* calculate H and V */
    for( int i = 0; i <= 7; i++ )
    {
        H += ( i + 1 ) * ( src[ 8 + i - FDEC_STRIDE ] - src[6 -i -FDEC_STRIDE] );
        V += ( i + 1 ) * ( src[-1 + (8+i)*FDEC_STRIDE] - src[-1 + (6-i)*FDEC_STRIDE] );
    }

    int a = 16 * ( src[-1 + 15*FDEC_STRIDE] + src[15 - FDEC_STRIDE] );
    int b = ( 5 * H + 32 ) >> 6;
    int c = ( 5 * V + 32 ) >> 6;

    int i00 = a - b * 7 - c * 7 + 16;

    for( int y = 0; y < 16; y++ )
    {
        int pix = i00;
        for( int x = 0; x < 16; x++ )
        {
            src[x] = x264_clip_pixel( pix>>5 );
            pix += b;
        }
        src += FDEC_STRIDE;
        i00 += c;
    }
}


/****************************************************************************
 * 8x8 prediction for intra chroma block (4:2:0)
 ****************************************************************************/

static void x264_predict_8x8c_dc_128_c( pixel *src )
{
    for( int y = 0; y < 8; y++ )
    {
        MPIXEL_X4( src+0 ) = PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) );
        MPIXEL_X4( src+4 ) = PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) );
        src += FDEC_STRIDE;
    }
}
static void x264_predict_8x8c_dc_left_c( pixel *src )
{
    int dc0 = 0, dc1 = 0;

    for( int y = 0; y < 4; y++ )
    {
        dc0 += src[y * FDEC_STRIDE     - 1];
        dc1 += src[(y+4) * FDEC_STRIDE - 1];
    }
    pixel4 dc0splat = PIXEL_SPLAT_X4( ( dc0 + 2 ) >> 2 );
    pixel4 dc1splat = PIXEL_SPLAT_X4( ( dc1 + 2 ) >> 2 );

    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc0splat;
        MPIXEL_X4( src+4 ) = dc0splat;
        src += FDEC_STRIDE;
    }
    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc1splat;
        MPIXEL_X4( src+4 ) = dc1splat;
        src += FDEC_STRIDE;
    }

}
static void x264_predict_8x8c_dc_top_c( pixel *src )
{
    int dc0 = 0, dc1 = 0;

    for( int x = 0; x < 4; x++ )
    {
        dc0 += src[x     - FDEC_STRIDE];
        dc1 += src[x + 4 - FDEC_STRIDE];
    }
    pixel4 dc0splat = PIXEL_SPLAT_X4( ( dc0 + 2 ) >> 2 );
    pixel4 dc1splat = PIXEL_SPLAT_X4( ( dc1 + 2 ) >> 2 );

    for( int y = 0; y < 8; y++ )
    {
        MPIXEL_X4( src+0 ) = dc0splat;
        MPIXEL_X4( src+4 ) = dc1splat;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x8c_dc_c( pixel *src )
{
    int s0 = 0, s1 = 0, s2 = 0, s3 = 0;

    /*
          s0 s1
       s2
       s3
    */
    for( int i = 0; i < 4; i++ )
    {
        s0 += src[i - FDEC_STRIDE];
        s1 += src[i + 4 - FDEC_STRIDE];
        s2 += src[-1 + i * FDEC_STRIDE];
        s3 += src[-1 + (i+4)*FDEC_STRIDE];
    }
    /*
       dc0 dc1
       dc2 dc3
     */
    pixel4 dc0 = PIXEL_SPLAT_X4( ( s0 + s2 + 4 ) >> 3 );
    pixel4 dc1 = PIXEL_SPLAT_X4( ( s1 + 2 ) >> 2 );
    pixel4 dc2 = PIXEL_SPLAT_X4( ( s3 + 2 ) >> 2 );
    pixel4 dc3 = PIXEL_SPLAT_X4( ( s1 + s3 + 4 ) >> 3 );

    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc0;
        MPIXEL_X4( src+4 ) = dc1;
        src += FDEC_STRIDE;
    }

    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc2;
        MPIXEL_X4( src+4 ) = dc3;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x8c_h_c( pixel *src )
{
    for( int i = 0; i < 8; i++ )
    {
        pixel4 v = PIXEL_SPLAT_X4( src[-1] );
        MPIXEL_X4( src+0 ) = v;
        MPIXEL_X4( src+4 ) = v;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x8c_v_c( pixel *src )
{
    pixel4 v0 = MPIXEL_X4( src+0-FDEC_STRIDE );
    pixel4 v1 = MPIXEL_X4( src+4-FDEC_STRIDE );

    for( int i = 0; i < 8; i++ )
    {
        MPIXEL_X4( src+0 ) = v0;
        MPIXEL_X4( src+4 ) = v1;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x8c_p_c( pixel *src )
{
    int H = 0, V = 0;

    for( int i = 0; i < 4; i++ )
    {
        H += ( i + 1 ) * ( src[4+i - FDEC_STRIDE] - src[2 - i -FDEC_STRIDE] );
        V += ( i + 1 ) * ( src[-1 +(i+4)*FDEC_STRIDE] - src[-1+(2-i)*FDEC_STRIDE] );
    }

    int a = 16 * ( src[-1+7*FDEC_STRIDE] + src[7 - FDEC_STRIDE] );
    int b = ( 17 * H + 16 ) >> 5;
    int c = ( 17 * V + 16 ) >> 5;
    int i00 = a -3*b -3*c + 16;

    for( int y = 0; y < 8; y++ )
    {
        int pix = i00;
        for( int x = 0; x < 8; x++ )
        {
            src[x] = x264_clip_pixel( pix>>5 );
            pix += b;
        }
        src += FDEC_STRIDE;
        i00 += c;
    }
}

/****************************************************************************
 * 8x16 prediction for intra chroma block (4:2:2)
 ****************************************************************************/

static void x264_predict_8x16c_dc_128_c( pixel *src )
{
    for( int y = 0; y < 16; y++ )
    {
        MPIXEL_X4( src+0 ) = PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) );
        MPIXEL_X4( src+4 ) = PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) );
        src += FDEC_STRIDE;
    }
}
static void x264_predict_8x16c_dc_left_c( pixel *src )
{
    for( int i = 0; i < 4; i++ )
    {
        int dc = 0;

        for( int y = 0; y < 4; y++ )
            dc += src[y*FDEC_STRIDE - 1];

        pixel4 dcsplat = PIXEL_SPLAT_X4( (dc + 2) >> 2 );

        for( int y = 0; y < 4; y++ )
        {
            MPIXEL_X4( src+0 ) = dcsplat;
            MPIXEL_X4( src+4 ) = dcsplat;
            src += FDEC_STRIDE;
        }
    }
}
static void x264_predict_8x16c_dc_top_c( pixel *src )
{
    int dc0 = 0, dc1 = 0;

    for(int  x = 0; x < 4; x++ )
    {
        dc0 += src[x     - FDEC_STRIDE];
        dc1 += src[x + 4 - FDEC_STRIDE];
    }
    pixel4 dc0splat = PIXEL_SPLAT_X4( ( dc0 + 2 ) >> 2 );
    pixel4 dc1splat = PIXEL_SPLAT_X4( ( dc1 + 2 ) >> 2 );

    for( int y = 0; y < 16; y++ )
    {
        MPIXEL_X4( src+0 ) = dc0splat;
        MPIXEL_X4( src+4 ) = dc1splat;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x16c_dc_c( pixel *src )
{
    int s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0;

    /*
          s0 s1
       s2
       s3
       s4
       s5
    */
    for( int i = 0; i < 4; i++ )
    {
        s0 += src[i+0 - FDEC_STRIDE];
        s1 += src[i+4 - FDEC_STRIDE];
        s2 += src[-1 + (i+0)  * FDEC_STRIDE];
        s3 += src[-1 + (i+4)  * FDEC_STRIDE];
        s4 += src[-1 + (i+8)  * FDEC_STRIDE];
        s5 += src[-1 + (i+12) * FDEC_STRIDE];
    }
    /*
       dc0 dc1
       dc2 dc3
       dc4 dc5
       dc6 dc7
    */
    pixel4 dc0 = PIXEL_SPLAT_X4( ( s0 + s2 + 4 ) >> 3 );
    pixel4 dc1 = PIXEL_SPLAT_X4( ( s1 + 2 ) >> 2 );
    pixel4 dc2 = PIXEL_SPLAT_X4( ( s3 + 2 ) >> 2 );
    pixel4 dc3 = PIXEL_SPLAT_X4( ( s1 + s3 + 4 ) >> 3 );
    pixel4 dc4 = PIXEL_SPLAT_X4( ( s4 + 2 ) >> 2 );
    pixel4 dc5 = PIXEL_SPLAT_X4( ( s1 + s4 + 4 ) >> 3 );
    pixel4 dc6 = PIXEL_SPLAT_X4( ( s5 + 2 ) >> 2 );
    pixel4 dc7 = PIXEL_SPLAT_X4( ( s1 + s5 + 4 ) >> 3 );

    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc0;
        MPIXEL_X4( src+4 ) = dc1;
        src += FDEC_STRIDE;
    }
    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc2;
        MPIXEL_X4( src+4 ) = dc3;
        src += FDEC_STRIDE;
    }
    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc4;
        MPIXEL_X4( src+4 ) = dc5;
        src += FDEC_STRIDE;
    }
    for( int y = 0; y < 4; y++ )
    {
        MPIXEL_X4( src+0 ) = dc6;
        MPIXEL_X4( src+4 ) = dc7;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x16c_h_c( pixel *src )
{
    for( int i = 0; i < 16; i++ )
    {
        pixel4 v = PIXEL_SPLAT_X4( src[-1] );
        MPIXEL_X4( src+0 ) = v;
        MPIXEL_X4( src+4 ) = v;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x16c_v_c( pixel *src )
{
    pixel4 v0 = MPIXEL_X4( src+0-FDEC_STRIDE );
    pixel4 v1 = MPIXEL_X4( src+4-FDEC_STRIDE );

    for( int i = 0; i < 16; i++ )
    {
        MPIXEL_X4( src+0 ) = v0;
        MPIXEL_X4( src+4 ) = v1;
        src += FDEC_STRIDE;
    }
}
void x264_predict_8x16c_p_c( pixel *src )
{
    int H = 0;
    int V = 0;

    for( int i = 0; i < 4; i++ )
        H += ( i + 1 ) * ( src[4 + i - FDEC_STRIDE] - src[2 - i - FDEC_STRIDE] );
    for( int i = 0; i < 8; i++ )
        V += ( i + 1 ) * ( src[-1 + (i+8)*FDEC_STRIDE] - src[-1 + (6-i)*FDEC_STRIDE] );

    int a = 16 * ( src[-1 + 15*FDEC_STRIDE] + src[7 - FDEC_STRIDE] );
    int b = ( 17 * H + 16 ) >> 5;
    int c = ( 5 * V + 32 ) >> 6;
    int i00 = a -3*b -7*c + 16;

    for( int y = 0; y < 16; y++ )
    {
        int pix = i00;
        for( int x = 0; x < 8; x++ )
        {
            src[x] = x264_clip_pixel( pix>>5 );
            pix += b;
        }
        src += FDEC_STRIDE;
        i00 += c;
    }
}

/****************************************************************************
 * 4x4 prediction for intra luma block
 ****************************************************************************/

#define SRC(x,y) src[(x)+(y)*FDEC_STRIDE]
#define SRC_X4(x,y) MPIXEL_X4( &SRC(x,y) )

#define PREDICT_4x4_DC(v)\
    SRC_X4(0,0) = SRC_X4(0,1) = SRC_X4(0,2) = SRC_X4(0,3) = v;

static void x264_predict_4x4_dc_128_c( pixel *src )
{
    PREDICT_4x4_DC( PIXEL_SPLAT_X4( 1 << (BIT_DEPTH-1) ) );
}
static void x264_predict_4x4_dc_left_c( pixel *src )
{
    pixel4 dc = PIXEL_SPLAT_X4( (SRC(-1,0) + SRC(-1,1) + SRC(-1,2) + SRC(-1,3) + 2) >> 2 );
    PREDICT_4x4_DC( dc );
}
static void x264_predict_4x4_dc_top_c( pixel *src )
{
    pixel4 dc = PIXEL_SPLAT_X4( (SRC(0,-1) + SRC(1,-1) + SRC(2,-1) + SRC(3,-1) + 2) >> 2 );
    PREDICT_4x4_DC( dc );
}
void x264_predict_4x4_dc_c( pixel *src )
{
    pixel4 dc = PIXEL_SPLAT_X4( (SRC(-1,0) + SRC(-1,1) + SRC(-1,2) + SRC(-1,3) +
                                 SRC(0,-1) + SRC(1,-1) + SRC(2,-1) + SRC(3,-1) + 4) >> 3 );
    PREDICT_4x4_DC( dc );
}
void x264_predict_4x4_h_c( pixel *src )
{
    SRC_X4(0,0) = PIXEL_SPLAT_X4( SRC(-1,0) );
    SRC_X4(0,1) = PIXEL_SPLAT_X4( SRC(-1,1) );
    SRC_X4(0,2) = PIXEL_SPLAT_X4( SRC(-1,2) );
    SRC_X4(0,3) = PIXEL_SPLAT_X4( SRC(-1,3) );
}
void x264_predict_4x4_v_c( pixel *src )
{
    PREDICT_4x4_DC(SRC_X4(0,-1));
}

#define PREDICT_4x4_LOAD_LEFT\
    int l0 = SRC(-1,0);\
    int l1 = SRC(-1,1);\
    int l2 = SRC(-1,2);\
    UNUSED int l3 = SRC(-1,3);

#define PREDICT_4x4_LOAD_TOP\
    int t0 = SRC(0,-1);\
    int t1 = SRC(1,-1);\
    int t2 = SRC(2,-1);\
    UNUSED int t3 = SRC(3,-1);

#define PREDICT_4x4_LOAD_TOP_RIGHT\
    int t4 = SRC(4,-1);\
    int t5 = SRC(5,-1);\
    int t6 = SRC(6,-1);\
    UNUSED int t7 = SRC(7,-1);

#define F1(a,b)   (((a)+(b)+1)>>1)
#define F2(a,b,c) (((a)+2*(b)+(c)+2)>>2)

static void x264_predict_4x4_ddl_c( pixel *src )
{
    PREDICT_4x4_LOAD_TOP
    PREDICT_4x4_LOAD_TOP_RIGHT
    SRC(0,0)= F2(t0,t1,t2);
    SRC(1,0)=SRC(0,1)= F2(t1,t2,t3);
    SRC(2,0)=SRC(1,1)=SRC(0,2)= F2(t2,t3,t4);
    SRC(3,0)=SRC(2,1)=SRC(1,2)=SRC(0,3)= F2(t3,t4,t5);
    SRC(3,1)=SRC(2,2)=SRC(1,3)= F2(t4,t5,t6);
    SRC(3,2)=SRC(2,3)= F2(t5,t6,t7);
    SRC(3,3)= F2(t6,t7,t7);
}
static void x264_predict_4x4_ddr_c( pixel *src )
{
    int lt = SRC(-1,-1);
    PREDICT_4x4_LOAD_LEFT
    PREDICT_4x4_LOAD_TOP
    SRC(3,0)= F2(t3,t2,t1);
    SRC(2,0)=SRC(3,1)= F2(t2,t1,t0);
    SRC(1,0)=SRC(2,1)=SRC(3,2)= F2(t1,t0,lt);
    SRC(0,0)=SRC(1,1)=SRC(2,2)=SRC(3,3)= F2(t0,lt,l0);
    SRC(0,1)=SRC(1,2)=SRC(2,3)= F2(lt,l0,l1);
    SRC(0,2)=SRC(1,3)= F2(l0,l1,l2);
    SRC(0,3)= F2(l1,l2,l3);
}

static void x264_predict_4x4_vr_c( pixel *src )
{
    int lt = SRC(-1,-1);
    PREDICT_4x4_LOAD_LEFT
    PREDICT_4x4_LOAD_TOP
    SRC(0,3)= F2(l2,l1,l0);
    SRC(0,2)= F2(l1,l0,lt);
    SRC(0,1)=SRC(1,3)= F2(l0,lt,t0);
    SRC(0,0)=SRC(1,2)= F1(lt,t0);
    SRC(1,1)=SRC(2,3)= F2(lt,t0,t1);
    SRC(1,0)=SRC(2,2)= F1(t0,t1);
    SRC(2,1)=SRC(3,3)= F2(t0,t1,t2);
    SRC(2,0)=SRC(3,2)= F1(t1,t2);
    SRC(3,1)= F2(t1,t2,t3);
    SRC(3,0)= F1(t2,t3);
}

static void x264_predict_4x4_hd_c( pixel *src )
{
    int lt= SRC(-1,-1);
    PREDICT_4x4_LOAD_LEFT
    PREDICT_4x4_LOAD_TOP
    SRC(0,3)= F1(l2,l3);
    SRC(1,3)= F2(l1,l2,l3);
    SRC(0,2)=SRC(2,3)= F1(l1,l2);
    SRC(1,2)=SRC(3,3)= F2(l0,l1,l2);
    SRC(0,1)=SRC(2,2)= F1(l0,l1);
    SRC(1,1)=SRC(3,2)= F2(lt,l0,l1);
    SRC(0,0)=SRC(2,1)= F1(lt,l0);
    SRC(1,0)=SRC(3,1)= F2(t0,lt,l0);
    SRC(2,0)= F2(t1,t0,lt);
    SRC(3,0)= F2(t2,t1,t0);
}

static void x264_predict_4x4_vl_c( pixel *src )
{
    PREDICT_4x4_LOAD_TOP
    PREDICT_4x4_LOAD_TOP_RIGHT
    SRC(0,0)= F1(t0,t1);
    SRC(0,1)= F2(t0,t1,t2);
    SRC(1,0)=SRC(0,2)= F1(t1,t2);
    SRC(1,1)=SRC(0,3)= F2(t1,t2,t3);
    SRC(2,0)=SRC(1,2)= F1(t2,t3);
    SRC(2,1)=SRC(1,3)= F2(t2,t3,t4);
    SRC(3,0)=SRC(2,2)= F1(t3,t4);
    SRC(3,1)=SRC(2,3)= F2(t3,t4,t5);
    SRC(3,2)= F1(t4,t5);
    SRC(3,3)= F2(t4,t5,t6);
}

static void x264_predict_4x4_hu_c( pixel *src )
{
    PREDICT_4x4_LOAD_LEFT
    SRC(0,0)= F1(l0,l1);
    SRC(1,0)= F2(l0,l1,l2);
    SRC(2,0)=SRC(0,1)= F1(l1,l2);
    SRC(3,0)=SRC(1,1)= F2(l1,l2,l3);
    SRC(2,1)=SRC(0,2)= F1(l2,l3);
    SRC(3,1)=SRC(1,2)= F2(l2,l3,l3);
    SRC(3,2)=SRC(1,3)=SRC(0,3)=
    SRC(2,2)=SRC(2,3)=SRC(3,3)= l3;
}

/****************************************************************************
 * Exported functions:
 ****************************************************************************/
void x264_predict_16x16_init( x264_predict_t pf[7] )
{
    pf[I_PRED_16x16_V ]     = x264_predict_16x16_v_c;
    pf[I_PRED_16x16_H ]     = x264_predict_16x16_h_c;
    pf[I_PRED_16x16_DC]     = x264_predict_16x16_dc_c;
    pf[I_PRED_16x16_P ]     = x264_predict_16x16_p_c;
    pf[I_PRED_16x16_DC_LEFT]= x264_predict_16x16_dc_left_c;
    pf[I_PRED_16x16_DC_TOP ]= x264_predict_16x16_dc_top_c;
    pf[I_PRED_16x16_DC_128 ]= x264_predict_16x16_dc_128_c;

}

void x264_predict_8x8c_init( x264_predict_t pf[7] )
{
    pf[I_PRED_CHROMA_V ]     = x264_predict_8x8c_v_c;
    pf[I_PRED_CHROMA_H ]     = x264_predict_8x8c_h_c;
    pf[I_PRED_CHROMA_DC]     = x264_predict_8x8c_dc_c;
    pf[I_PRED_CHROMA_P ]     = x264_predict_8x8c_p_c;
    pf[I_PRED_CHROMA_DC_LEFT]= x264_predict_8x8c_dc_left_c;
    pf[I_PRED_CHROMA_DC_TOP ]= x264_predict_8x8c_dc_top_c;
    pf[I_PRED_CHROMA_DC_128 ]= x264_predict_8x8c_dc_128_c;

}

void x264_predict_8x16c_init( x264_predict_t pf[7] )
{
    pf[I_PRED_CHROMA_V ]     = x264_predict_8x16c_v_c;
    pf[I_PRED_CHROMA_H ]     = x264_predict_8x16c_h_c;
    pf[I_PRED_CHROMA_DC]     = x264_predict_8x16c_dc_c;
    pf[I_PRED_CHROMA_P ]     = x264_predict_8x16c_p_c;
    pf[I_PRED_CHROMA_DC_LEFT]= x264_predict_8x16c_dc_left_c;
    pf[I_PRED_CHROMA_DC_TOP ]= x264_predict_8x16c_dc_top_c;
    pf[I_PRED_CHROMA_DC_128 ]= x264_predict_8x16c_dc_128_c;

}


void x264_predict_4x4_init( x264_predict_t pf[12] )
{
    pf[I_PRED_4x4_V]      = x264_predict_4x4_v_c;
    pf[I_PRED_4x4_H]      = x264_predict_4x4_h_c;
    pf[I_PRED_4x4_DC]     = x264_predict_4x4_dc_c;
    pf[I_PRED_4x4_DDL]    = x264_predict_4x4_ddl_c;
    pf[I_PRED_4x4_DDR]    = x264_predict_4x4_ddr_c;
    pf[I_PRED_4x4_VR]     = x264_predict_4x4_vr_c;
    pf[I_PRED_4x4_HD]     = x264_predict_4x4_hd_c;
    pf[I_PRED_4x4_VL]     = x264_predict_4x4_vl_c;
    pf[I_PRED_4x4_HU]     = x264_predict_4x4_hu_c;
    pf[I_PRED_4x4_DC_LEFT]= x264_predict_4x4_dc_left_c;
    pf[I_PRED_4x4_DC_TOP] = x264_predict_4x4_dc_top_c;
    pf[I_PRED_4x4_DC_128] = x264_predict_4x4_dc_128_c;

}

