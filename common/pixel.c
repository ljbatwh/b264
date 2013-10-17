/*****************************************************************************
 * pixel.c: pixel metrics
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

#include "common.h"

/****************************************************************************
 * pixel_sad_WxH
 ****************************************************************************/
/*
#define PIXEL_SAD_C( name, lx, ly ) \
static int name( pixel *pix1, intptr_t i_stride_pix1,  \
                 pixel *pix2, intptr_t i_stride_pix2 ) \
{                                                   \
    int i_sum = 0;                                  \
    for( int y = 0; y < ly; y++ )                   \
    {                                               \
        for( int x = 0; x < lx; x++ )               \
        {                                           \
            i_sum += abs( pix1[x] - pix2[x] );      \
        }                                           \
        pix1 += i_stride_pix1;                      \
        pix2 += i_stride_pix2;                      \
    }                                               \
    return i_sum;                                   \
}


PIXEL_SAD_C( x264_pixel_sad_16x16, 16, 16 )
PIXEL_SAD_C( x264_pixel_sad_16x8,  16,  8 )
PIXEL_SAD_C( x264_pixel_sad_8x16,   8, 16 )
PIXEL_SAD_C( x264_pixel_sad_8x8,    8,  8 )
PIXEL_SAD_C( x264_pixel_sad_8x4,    8,  4 )
PIXEL_SAD_C( x264_pixel_sad_4x16,   4, 16 )
PIXEL_SAD_C( x264_pixel_sad_4x8,    4,  8 )
PIXEL_SAD_C( x264_pixel_sad_4x4,    4,  4 )
*/
static int x264_pixel_sad_16x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 16; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_16x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 16; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_8x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_8x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_8x4( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 4; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_4x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_4x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_sad_4x4( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 4; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            i_sum += abs( pix1[x] - pix2[x] );
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}

/****************************************************************************
 * pixel_ssd_WxH
 ****************************************************************************/
/*
#define PIXEL_SSD_C( name, lx, ly ) \
static int name( pixel *pix1, intptr_t i_stride_pix1,  \
                 pixel *pix2, intptr_t i_stride_pix2 ) \
{                                                   \
    int i_sum = 0;                                  \
    for( int y = 0; y < ly; y++ )                   \
    {                                               \
        for( int x = 0; x < lx; x++ )               \
        {                                           \
            int d = pix1[x] - pix2[x];              \
            i_sum += d*d;                           \
        }                                           \
        pix1 += i_stride_pix1;                      \
        pix2 += i_stride_pix2;                      \
    }                                               \
    return i_sum;                                   \
}

PIXEL_SSD_C( x264_pixel_ssd_16x16, 16, 16 )
PIXEL_SSD_C( x264_pixel_ssd_16x8,  16,  8 )
PIXEL_SSD_C( x264_pixel_ssd_8x16,   8, 16 )
PIXEL_SSD_C( x264_pixel_ssd_8x8,    8,  8 )
PIXEL_SSD_C( x264_pixel_ssd_8x4,    8,  4 )
PIXEL_SSD_C( x264_pixel_ssd_4x16,   4, 16 )
PIXEL_SSD_C( x264_pixel_ssd_4x8,    4,  8 )
PIXEL_SSD_C( x264_pixel_ssd_4x4,    4,  4 )
*/

static int x264_pixel_ssd_16x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 16; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_16x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 16; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_8x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_8x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_8x4( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 4; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_4x16( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_4x8( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}
static int x264_pixel_ssd_4x4( pixel *pix1, intptr_t i_stride_pix1, pixel *pix2, intptr_t i_stride_pix2 ) {
    int i_sum = 0;
    for ( int y = 0; y < 4; y++ ) {
        for ( int x = 0; x < 4; x++ ) {
            int d = pix1[x] - pix2[x];
            i_sum += d*d;
        }
        pix1 += i_stride_pix1;
        pix2 += i_stride_pix2;
    }
    return i_sum;
}

uint64_t x264_pixel_ssd_wxh( x264_pixel_function_t *pf, pixel *pix1, intptr_t i_pix1,
                             pixel *pix2, intptr_t i_pix2, int i_width, int i_height )
{
    uint64_t i_ssd = 0;
    int y;
    int align = !(((intptr_t)pix1 | (intptr_t)pix2 | i_pix1 | i_pix2) & 15);

#define SSD(size) i_ssd += pf->ssd[size]( pix1 + y*i_pix1 + x, i_pix1, \
                                          pix2 + y*i_pix2 + x, i_pix2 );
    for( y = 0; y < i_height-15; y += 16 )
    {
        int x = 0;
        if( align )
            for( ; x < i_width-15; x += 16 )
                SSD(PIXEL_16x16);
        for( ; x < i_width-7; x += 8 )
            SSD(PIXEL_8x16);
    }
    if( y < i_height-7 )
        for( int x = 0; x < i_width-7; x += 8 )
            SSD(PIXEL_8x8);
#undef SSD

#define SSD1 { int d = pix1[y*i_pix1+x] - pix2[y*i_pix2+x]; i_ssd += d*d; }
    if( i_width & 7 )
    {
        for( y = 0; y < (i_height & ~7); y++ )
            for( int x = i_width & ~7; x < i_width; x++ )
                SSD1;
    }
    if( i_height & 7 )
    {
        for( y = i_height & ~7; y < i_height; y++ )
            for( int x = 0; x < i_width; x++ )
                SSD1;
    }
#undef SSD1

    return i_ssd;
}

static void pixel_ssd_nv12_core( pixel *pixuv1, intptr_t stride1, pixel *pixuv2, intptr_t stride2,
                                 int width, int height, uint64_t *ssd_u, uint64_t *ssd_v )
{
    *ssd_u = 0, *ssd_v = 0;
    for( int y = 0; y < height; y++, pixuv1+=stride1, pixuv2+=stride2 )
        for( int x = 0; x < width; x++ )
        {
            int du = pixuv1[2*x]   - pixuv2[2*x];
            int dv = pixuv1[2*x+1] - pixuv2[2*x+1];
            *ssd_u += du*du;
            *ssd_v += dv*dv;
        }
}

void x264_pixel_ssd_nv12( x264_pixel_function_t *pf, pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2,
                          int i_width, int i_height, uint64_t *ssd_u, uint64_t *ssd_v )
{
    pf->ssd_nv12_core( pix1, i_pix1, pix2, i_pix2, i_width&~7, i_height, ssd_u, ssd_v );
    if( i_width&7 )
    {
        uint64_t tmp[2];
        pixel_ssd_nv12_core( pix1+(i_width&~7), i_pix1, pix2+(i_width&~7), i_pix2, i_width&7, i_height, &tmp[0], &tmp[1] );
        *ssd_u += tmp[0];
        *ssd_v += tmp[1];
    }
}

/****************************************************************************
 * pixel_var_wxh
 ****************************************************************************/
/*
#define PIXEL_VAR_C( name, w, h ) \
static uint64_t name( pixel *pix, intptr_t i_stride ) \
{                                             \
    uint32_t sum = 0, sqr = 0;                \
    for( int y = 0; y < h; y++ )              \
    {                                         \
        for( int x = 0; x < w; x++ )          \
        {                                     \
            sum += pix[x];                    \
            sqr += pix[x] * pix[x];           \
        }                                     \
        pix += i_stride;                      \
    }                                         \
    return sum + ((uint64_t)sqr << 32);       \
}

PIXEL_VAR_C( x264_pixel_var_16x16, 16, 16 )
PIXEL_VAR_C( x264_pixel_var_8x16,   8, 16 )
PIXEL_VAR_C( x264_pixel_var_8x8,    8,  8 )
*/
static uint64_t x264_pixel_var_16x16( pixel *pix, intptr_t i_stride ) {
    uint32_t sum = 0, sqr = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 16; x++ ) {
            sum += pix[x];
            sqr += pix[x] * pix[x];
        }
        pix += i_stride;
    }
    return sum + ((uint64_t)sqr << 32);
}
static uint64_t x264_pixel_var_8x16( pixel *pix, intptr_t i_stride ) {
    uint32_t sum = 0, sqr = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            sum += pix[x];
            sqr += pix[x] * pix[x];
        }
        pix += i_stride;
    }
    return sum + ((uint64_t)sqr << 32);
}
static uint64_t x264_pixel_var_8x8( pixel *pix, intptr_t i_stride ) {
    uint32_t sum = 0, sqr = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            sum += pix[x];
            sqr += pix[x] * pix[x];
        }
        pix += i_stride;
    }
    return sum + ((uint64_t)sqr << 32);
}

/****************************************************************************
 * pixel_var2_wxh
 ****************************************************************************/
/*
#define PIXEL_VAR2_C( name, w, h, shift ) \
static int name( pixel *pix1, intptr_t i_stride1, pixel *pix2, intptr_t i_stride2, int *ssd ) \
{ \
    uint32_t var = 0, sum = 0, sqr = 0; \
    for( int y = 0; y < h; y++ ) \
    { \
        for( int x = 0; x < w; x++ ) \
        { \
            int diff = pix1[x] - pix2[x]; \
            sum += diff; \
            sqr += diff * diff; \
        } \
        pix1 += i_stride1; \
        pix2 += i_stride2; \
    } \
    sum = abs(sum); \
    var = sqr - ((uint64_t)sum * sum >> shift); \
    *ssd = sqr; \
    return var; \
}

PIXEL_VAR2_C( x264_pixel_var2_8x16, 8, 16, 7 )
PIXEL_VAR2_C( x264_pixel_var2_8x8,  8,  8, 6 )
*/
static int x264_pixel_var2_8x16( pixel *pix1, intptr_t i_stride1, pixel *pix2, intptr_t i_stride2, int *ssd ) {
    uint32_t var = 0, sum = 0, sqr = 0;
    for ( int y = 0; y < 16; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            int diff = pix1[x] - pix2[x];
            sum += diff;
            sqr += diff * diff;
        }
        pix1 += i_stride1;
        pix2 += i_stride2;
    }
    sum = abs(sum);
    var = sqr - ((uint64_t)sum * sum >> 7);
    *ssd = sqr;
    return var;
}
static int x264_pixel_var2_8x8( pixel *pix1, intptr_t i_stride1, pixel *pix2, intptr_t i_stride2, int *ssd ) {
    uint32_t var = 0, sum = 0, sqr = 0;
    for ( int y = 0; y < 8; y++ ) {
        for ( int x = 0; x < 8; x++ ) {
            int diff = pix1[x] - pix2[x];
            sum += diff;
            sqr += diff * diff;
        }
        pix1 += i_stride1;
        pix2 += i_stride2;
    }
    sum = abs(sum);
    var = sqr - ((uint64_t)sum * sum >> 6);
    *ssd = sqr;
    return var;
}

typedef uint16_t sum_t;
typedef uint32_t sum2_t;
#define BITS_PER_SUM (8 * sizeof(sum_t))

#define HADAMARD4(d0, d1, d2, d3, s0, s1, s2, s3) {\
    sum2_t t0 = s0 + s1;\
    sum2_t t1 = s0 - s1;\
    sum2_t t2 = s2 + s3;\
    sum2_t t3 = s2 - s3;\
    d0 = t0 + t2;\
    d2 = t0 - t2;\
    d1 = t1 + t3;\
    d3 = t1 - t3;\
}

// in: a pseudo-simd number of the form x+(y<<16)
// return: abs(x)+(abs(y)<<16)
static ALWAYS_INLINE sum2_t abs2( sum2_t a )
{
    sum2_t s = ((a>>(BITS_PER_SUM-1))&(((sum2_t)1<<BITS_PER_SUM)+1))*((sum_t)-1);
    return (a+s)^s;
}

/****************************************************************************
 * pixel_satd_WxH: sum of 4x4 Hadamard transformed differences
 ****************************************************************************/

static NOINLINE int x264_pixel_satd_4x4( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )
{
    sum2_t tmp[4][2];
    sum2_t a0, a1, a2, a3, b0, b1;
    sum2_t sum = 0;
    for( int i = 0; i < 4; i++, pix1 += i_pix1, pix2 += i_pix2 )
    {
        a0 = pix1[0] - pix2[0];
        a1 = pix1[1] - pix2[1];
        b0 = (a0+a1) + ((a0-a1)<<BITS_PER_SUM);
        a2 = pix1[2] - pix2[2];
        a3 = pix1[3] - pix2[3];
        b1 = (a2+a3) + ((a2-a3)<<BITS_PER_SUM);
        tmp[i][0] = b0 + b1;
        tmp[i][1] = b0 - b1;
    }
    for( int i = 0; i < 2; i++ )
    {
        HADAMARD4( a0, a1, a2, a3, tmp[0][i], tmp[1][i], tmp[2][i], tmp[3][i] );
        a0 = abs2(a0) + abs2(a1) + abs2(a2) + abs2(a3);
        sum += ((sum_t)a0) + (a0>>BITS_PER_SUM);
    }
    return sum >> 1;
}

static NOINLINE int x264_pixel_satd_8x4( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )
{
    sum2_t tmp[4][4];
    sum2_t a0, a1, a2, a3;
    sum2_t sum = 0;
    for( int i = 0; i < 4; i++, pix1 += i_pix1, pix2 += i_pix2 )
    {
        a0 = (pix1[0] - pix2[0]) + ((sum2_t)(pix1[4] - pix2[4]) << BITS_PER_SUM);
        a1 = (pix1[1] - pix2[1]) + ((sum2_t)(pix1[5] - pix2[5]) << BITS_PER_SUM);
        a2 = (pix1[2] - pix2[2]) + ((sum2_t)(pix1[6] - pix2[6]) << BITS_PER_SUM);
        a3 = (pix1[3] - pix2[3]) + ((sum2_t)(pix1[7] - pix2[7]) << BITS_PER_SUM);
        HADAMARD4( tmp[i][0], tmp[i][1], tmp[i][2], tmp[i][3], a0,a1,a2,a3 );
    }
    for( int i = 0; i < 4; i++ )
    {
        HADAMARD4( a0, a1, a2, a3, tmp[0][i], tmp[1][i], tmp[2][i], tmp[3][i] );
        sum += abs2(a0) + abs2(a1) + abs2(a2) + abs2(a3);
    }
    return (((sum_t)sum) + (sum>>BITS_PER_SUM)) >> 1;
}
/*
#define PIXEL_SATD_C( w, h, sub )\
static int x264_pixel_satd_##w##x##h( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )\
{\
    int sum = sub( pix1, i_pix1, pix2, i_pix2 )\
            + sub( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );\
    if( w==16 )\
        sum+= sub( pix1+8, i_pix1, pix2+8, i_pix2 )\
            + sub( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );\
    if( h==16 )\
        sum+= sub( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 )\
            + sub( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );\
    if( w==16 && h==16 )\
        sum+= sub( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 )\
            + sub( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );\
    return sum;\
}
PIXEL_SATD_C( 16, 16, x264_pixel_satd_8x4 )
PIXEL_SATD_C( 16, 8,  x264_pixel_satd_8x4 )
PIXEL_SATD_C( 8,  16, x264_pixel_satd_8x4 )
PIXEL_SATD_C( 8,  8,  x264_pixel_satd_8x4 )
PIXEL_SATD_C( 4,  16, x264_pixel_satd_4x4 )
PIXEL_SATD_C( 4,  8,  x264_pixel_satd_4x4 )
*/
static int x264_pixel_satd_16x16( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_8x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 16==16 && 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}
static int x264_pixel_satd_16x8( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_8x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 16==16 && 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}
static int x264_pixel_satd_8x16( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_8x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 8==16 && 16==16 ) sum+= x264_pixel_satd_8x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}
static int x264_pixel_satd_8x8( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_8x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 8==16 && 8==16 ) sum+= x264_pixel_satd_8x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_8x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}
static int x264_pixel_satd_4x16( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_4x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 4==16 ) sum+= x264_pixel_satd_4x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_4x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 16==16 ) sum+= x264_pixel_satd_4x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 4==16 && 16==16 ) sum+= x264_pixel_satd_4x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}
static int x264_pixel_satd_4x8( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 ) {
    int sum = x264_pixel_satd_4x4( pix1, i_pix1, pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+4*i_pix1, i_pix1, pix2+4*i_pix2, i_pix2 );
    if ( 4==16 ) sum+= x264_pixel_satd_4x4( pix1+8, i_pix1, pix2+8, i_pix2 ) + x264_pixel_satd_4x4( pix1+8+4*i_pix1, i_pix1, pix2+8+4*i_pix2, i_pix2 );
    if ( 8==16 ) sum+= x264_pixel_satd_4x4( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+12*i_pix1, i_pix1, pix2+12*i_pix2, i_pix2 );
    if ( 4==16 && 8==16 ) sum+= x264_pixel_satd_4x4( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 ) + x264_pixel_satd_4x4( pix1+8+12*i_pix1, i_pix1, pix2+8+12*i_pix2, i_pix2 );
    return sum;
}

static NOINLINE int sa8d_8x8( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )
{
    sum2_t tmp[8][4];
    sum2_t a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3;
    sum2_t sum = 0;
    for( int i = 0; i < 8; i++, pix1 += i_pix1, pix2 += i_pix2 )
    {
        a0 = pix1[0] - pix2[0];
        a1 = pix1[1] - pix2[1];
        b0 = (a0+a1) + ((a0-a1)<<BITS_PER_SUM);
        a2 = pix1[2] - pix2[2];
        a3 = pix1[3] - pix2[3];
        b1 = (a2+a3) + ((a2-a3)<<BITS_PER_SUM);
        a4 = pix1[4] - pix2[4];
        a5 = pix1[5] - pix2[5];
        b2 = (a4+a5) + ((a4-a5)<<BITS_PER_SUM);
        a6 = pix1[6] - pix2[6];
        a7 = pix1[7] - pix2[7];
        b3 = (a6+a7) + ((a6-a7)<<BITS_PER_SUM);
        HADAMARD4( tmp[i][0], tmp[i][1], tmp[i][2], tmp[i][3], b0,b1,b2,b3 );
    }
    for( int i = 0; i < 4; i++ )
    {
        HADAMARD4( a0, a1, a2, a3, tmp[0][i], tmp[1][i], tmp[2][i], tmp[3][i] );
        HADAMARD4( a4, a5, a6, a7, tmp[4][i], tmp[5][i], tmp[6][i], tmp[7][i] );
        b0  = abs2(a0+a4) + abs2(a0-a4);
        b0 += abs2(a1+a5) + abs2(a1-a5);
        b0 += abs2(a2+a6) + abs2(a2-a6);
        b0 += abs2(a3+a7) + abs2(a3-a7);
        sum += (sum_t)b0 + (b0>>BITS_PER_SUM);
    }
    return sum;
}

static int x264_pixel_sa8d_8x8( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )
{
    int sum = sa8d_8x8( pix1, i_pix1, pix2, i_pix2 );
    return (sum+2)>>2;
}

static int x264_pixel_sa8d_16x16( pixel *pix1, intptr_t i_pix1, pixel *pix2, intptr_t i_pix2 )
{
    int sum = sa8d_8x8( pix1, i_pix1, pix2, i_pix2 )
            + sa8d_8x8( pix1+8, i_pix1, pix2+8, i_pix2 )
            + sa8d_8x8( pix1+8*i_pix1, i_pix1, pix2+8*i_pix2, i_pix2 )
            + sa8d_8x8( pix1+8+8*i_pix1, i_pix1, pix2+8+8*i_pix2, i_pix2 );
    return (sum+2)>>2;
}
//8x8 hadamard and remove dc
static NOINLINE uint64_t pixel_hadamard_ac( pixel *pix, intptr_t stride )
{
    sum2_t tmp[32];
    sum2_t a0, a1, a2, a3, dc;
    sum2_t sum4 = 0, sum8 = 0;
    for( int i = 0; i < 8; i++, pix+=stride )
    {
        sum2_t *t = tmp + (i&3) + (i&4)*4;
        a0 = (pix[0]+pix[1]) + ((sum2_t)(pix[0]-pix[1])<<BITS_PER_SUM);
        a1 = (pix[2]+pix[3]) + ((sum2_t)(pix[2]-pix[3])<<BITS_PER_SUM);
        t[0] = a0 + a1;
        t[4] = a0 - a1;
        a2 = (pix[4]+pix[5]) + ((sum2_t)(pix[4]-pix[5])<<BITS_PER_SUM);
        a3 = (pix[6]+pix[7]) + ((sum2_t)(pix[6]-pix[7])<<BITS_PER_SUM);
        t[8] = a2 + a3;
        t[12] = a2 - a3;
    }
    for( int i = 0; i < 8; i++ )
    {
        HADAMARD4( a0, a1, a2, a3, tmp[i*4+0], tmp[i*4+1], tmp[i*4+2], tmp[i*4+3] );
        tmp[i*4+0] = a0;
        tmp[i*4+1] = a1;
        tmp[i*4+2] = a2;
        tmp[i*4+3] = a3;
        sum4 += abs2(a0) + abs2(a1) + abs2(a2) + abs2(a3);
    }
    for( int i = 0; i < 8; i++ )
    {
        HADAMARD4( a0,a1,a2,a3, tmp[i], tmp[8+i], tmp[16+i], tmp[24+i] );
        sum8 += abs2(a0) + abs2(a1) + abs2(a2) + abs2(a3);
    }
    dc = (sum_t)(tmp[0] + tmp[8] + tmp[16] + tmp[24]);
    sum4 = (sum_t)sum4 + (sum4>>BITS_PER_SUM) - dc;
    sum8 = (sum_t)sum8 + (sum8>>BITS_PER_SUM) - dc;
    return ((uint64_t)sum8<<32) + sum4;
}
/*
#define HADAMARD_AC(w,h) \
static uint64_t x264_pixel_hadamard_ac_##w##x##h( pixel *pix, intptr_t stride )\
{\
    uint64_t sum = pixel_hadamard_ac( pix, stride );\
    if( w==16 )\
        sum += pixel_hadamard_ac( pix+8, stride );\
    if( h==16 )\
        sum += pixel_hadamard_ac( pix+8*stride, stride );\
    if( w==16 && h==16 )\
        sum += pixel_hadamard_ac( pix+8*stride+8, stride );\
    return ((sum>>34)<<32) + ((uint32_t)sum>>1);\
}
HADAMARD_AC( 16, 16 )
HADAMARD_AC( 16, 8 )
HADAMARD_AC( 8, 16 )
HADAMARD_AC( 8, 8 )
*/
//using divide and conquer algorithm
//16x16 hadamard divide into four 8x8 hadamard
static uint64_t x264_pixel_hadamard_ac_16x16( pixel *pix, intptr_t stride ) {
    uint64_t sum = pixel_hadamard_ac( pix, stride );
    sum += pixel_hadamard_ac( pix+8, stride );
    sum += pixel_hadamard_ac( pix+8*stride, stride );
    sum += pixel_hadamard_ac( pix+8*stride+8, stride );
    return ((sum>>34)<<32) + ((uint32_t)sum>>1);
}
//16x8 hadamard divide into two 8x8 hadamard
static uint64_t x264_pixel_hadamard_ac_16x8( pixel *pix, intptr_t stride ) {
    uint64_t sum = pixel_hadamard_ac( pix, stride );
    sum += pixel_hadamard_ac( pix+8, stride );
    return ((sum>>34)<<32) + ((uint32_t)sum>>1);
}
static uint64_t x264_pixel_hadamard_ac_8x16( pixel *pix, intptr_t stride ) {
    uint64_t sum = pixel_hadamard_ac( pix, stride );
    sum += pixel_hadamard_ac( pix+8*stride, stride );
    return ((sum>>34)<<32) + ((uint32_t)sum>>1);
}
static uint64_t x264_pixel_hadamard_ac_8x8( pixel *pix, intptr_t stride ) {
    uint64_t sum = pixel_hadamard_ac( pix, stride );
    return ((sum>>34)<<32) + ((uint32_t)sum>>1);
}


/****************************************************************************
 * pixel_sad_x4
 ****************************************************************************/
/*
#define SAD_X( size ) \
static void x264_pixel_sad_x3_##size( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2,\
                                      intptr_t i_stride, int scores[3] )\
{\
    scores[0] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix0, i_stride );\
    scores[1] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix1, i_stride );\
    scores[2] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix2, i_stride );\
}\
static void x264_pixel_sad_x4_##size( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3,\
                                      intptr_t i_stride, int scores[4] )\
{\
    scores[0] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix0, i_stride );\
    scores[1] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix1, i_stride );\
    scores[2] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix2, i_stride );\
    scores[3] = x264_pixel_sad_##size( fenc, FENC_STRIDE, pix3, i_stride );\
}

SAD_X( 16x16 )
SAD_X( 16x8 )
SAD_X( 8x16 )
SAD_X( 8x8 )
SAD_X( 8x4 )
SAD_X( 4x8 )
SAD_X( 4x4 )

#if !HIGH_BIT_DEPTH
#if ARCH_UltraSPARC
SAD_X( 16x16_vis )
SAD_X( 16x8_vis )
SAD_X( 8x16_vis )
SAD_X( 8x8_vis )
#endif
#endif // !HIGH_BIT_DEPTH
*/
static void x264_pixel_sad_x3_16x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_16x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_16x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_16x16( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_16x16( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_16x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_16x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_16x16( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_16x16( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_16x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_16x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_16x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_16x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_16x8( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_16x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_16x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_16x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_16x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_8x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_8x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x16( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_8x16( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_8x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x16( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_8x16( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_8x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_8x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_8x8( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_8x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_8x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_8x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_8x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x4( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_8x4( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_8x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_8x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_8x4( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_8x4( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_4x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_4x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_4x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_4x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_4x8( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_4x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_4x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_4x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_4x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_sad_x3_4x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_sad_4x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_4x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_4x4( fenc, 16, pix2, i_stride );
}
static void x264_pixel_sad_x4_4x4( pixel *fenc, pixel *pix0, pixel *pix1,pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_sad_4x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_sad_4x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_sad_4x4( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_sad_4x4( fenc, 16, pix3, i_stride );
}

/****************************************************************************
 * pixel_satd_x4
 * no faster than single satd, but needed for satd to be a drop-in replacement for sad
 ****************************************************************************/
/*
#define SATD_X( size, cpu ) \
static void x264_pixel_satd_x3_##size##cpu( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2,\
                                            intptr_t i_stride, int scores[3] )\
{\
    scores[0] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix0, i_stride );\
    scores[1] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix1, i_stride );\
    scores[2] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix2, i_stride );\
}\
static void x264_pixel_satd_x4_##size##cpu( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3,\
                                            intptr_t i_stride, int scores[4] )\
{\
    scores[0] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix0, i_stride );\
    scores[1] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix1, i_stride );\
    scores[2] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix2, i_stride );\
    scores[3] = x264_pixel_satd_##size##cpu( fenc, FENC_STRIDE, pix3, i_stride );\
}
#define SATD_X_DECL6( cpu )\
SATD_X( 16x16, cpu )\
SATD_X( 16x8, cpu )\
SATD_X( 8x16, cpu )\
SATD_X( 8x8, cpu )\
SATD_X( 8x4, cpu )\
SATD_X( 4x8, cpu )
#define SATD_X_DECL7( cpu )\
SATD_X_DECL6( cpu )\
SATD_X( 4x4, cpu )

SATD_X_DECL7()

    

*/

static void x264_pixel_satd_x3_16x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_16x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_16x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_16x16( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_16x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_16x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_16x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_16x16( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_16x16( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_16x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_16x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_16x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_16x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_16x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_16x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_16x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_16x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_16x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_8x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_8x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x16( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_8x16( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_8x16( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x16( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x16( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_8x16( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_8x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_8x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_8x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_8x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_8x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_8x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_8x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x4( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_8x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_8x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_8x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_8x4( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_8x4( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_4x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_4x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_4x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_4x8( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_4x8( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_4x8( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_4x8( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_4x8( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_4x8( fenc, 16, pix3, i_stride );
}
static void x264_pixel_satd_x3_4x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, intptr_t i_stride, int scores[3] ) {
    scores[0] = x264_pixel_satd_4x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_4x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_4x4( fenc, 16, pix2, i_stride );
}
static void x264_pixel_satd_x4_4x4( pixel *fenc, pixel *pix0, pixel *pix1, pixel *pix2, pixel *pix3, intptr_t i_stride, int scores[4] ) {
    scores[0] = x264_pixel_satd_4x4( fenc, 16, pix0, i_stride );
    scores[1] = x264_pixel_satd_4x4( fenc, 16, pix1, i_stride );
    scores[2] = x264_pixel_satd_4x4( fenc, 16, pix2, i_stride );
    scores[3] = x264_pixel_satd_4x4( fenc, 16, pix3, i_stride );
}
/*
#define INTRA_MBCMP_8x8( mbcmp, cpu, cpu2 )\
void x264_intra_##mbcmp##_x3_8x8##cpu( pixel *fenc, pixel edge[36], int res[3] )\
{\
    ALIGNED_ARRAY_16( pixel, pix, [8*FDEC_STRIDE] );\
    x264_predict_8x8_v##cpu2( pix, edge );\
    res[0] = x264_pixel_##mbcmp##_8x8##cpu( pix, FDEC_STRIDE, fenc, FENC_STRIDE );\
    x264_predict_8x8_h##cpu2( pix, edge );\
    res[1] = x264_pixel_##mbcmp##_8x8##cpu( pix, FDEC_STRIDE, fenc, FENC_STRIDE );\
    x264_predict_8x8_dc##cpu2( pix, edge );\
    res[2] = x264_pixel_##mbcmp##_8x8##cpu( pix, FDEC_STRIDE, fenc, FENC_STRIDE );\
}

INTRA_MBCMP_8x8( sad,, _c )
INTRA_MBCMP_8x8(sa8d,, _c )

*/

void x264_intra_sad_x3_8x8( pixel *fenc, pixel edge[36], int res[3] ) {
    pixel pix [8*32] __attribute__((aligned(16)));
    x264_predict_8x8_v_c( pix, edge );
    res[0] = x264_pixel_sad_8x8( pix, 32, fenc, 16 );
    x264_predict_8x8_h_c( pix, edge );
    res[1] = x264_pixel_sad_8x8( pix, 32, fenc, 16 );
    x264_predict_8x8_dc_c( pix, edge );
    res[2] = x264_pixel_sad_8x8( pix, 32, fenc, 16 );
}
void x264_intra_sa8d_x3_8x8( pixel *fenc, pixel edge[36], int res[3] ) {
    pixel pix [8*32] __attribute__((aligned(16)));
    x264_predict_8x8_v_c( pix, edge );
    res[0] = x264_pixel_sa8d_8x8( pix, 32, fenc, 16 );
    x264_predict_8x8_h_c( pix, edge );
    res[1] = x264_pixel_sa8d_8x8( pix, 32, fenc, 16 );
    x264_predict_8x8_dc_c( pix, edge );
    res[2] = x264_pixel_sa8d_8x8( pix, 32, fenc, 16 );
}
/*
#define INTRA_MBCMP( mbcmp, size, pred1, pred2, pred3, chroma, cpu, cpu2 )\
void x264_intra_##mbcmp##_x3_##size##chroma##cpu( pixel *fenc, pixel *fdec, int res[3] )\
{\
    x264_predict_##size##chroma##_##pred1##cpu2( fdec );\
    res[0] = x264_pixel_##mbcmp##_##size##cpu( fdec, FDEC_STRIDE, fenc, FENC_STRIDE );\
    x264_predict_##size##chroma##_##pred2##cpu2( fdec );\
    res[1] = x264_pixel_##mbcmp##_##size##cpu( fdec, FDEC_STRIDE, fenc, FENC_STRIDE );\
    x264_predict_##size##chroma##_##pred3##cpu2( fdec );\
    res[2] = x264_pixel_##mbcmp##_##size##cpu( fdec, FDEC_STRIDE, fenc, FENC_STRIDE );\
}

INTRA_MBCMP( sad,  4x4,   v, h, dc,  ,, _c )
INTRA_MBCMP(satd,  4x4,   v, h, dc,  ,, _c )
INTRA_MBCMP( sad,  8x8,  dc, h,  v, c,, _c )
INTRA_MBCMP(satd,  8x8,  dc, h,  v, c,, _c )
INTRA_MBCMP( sad,  8x16, dc, h,  v, c,, _c )
INTRA_MBCMP(satd,  8x16, dc, h,  v, c,, _c )
INTRA_MBCMP( sad, 16x16,  v, h, dc,  ,, _c )
INTRA_MBCMP(satd, 16x16,  v, h, dc,  ,, _c )


*/
void x264_intra_sad_x3_4x4( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_4x4_v_c( fdec );
    res[0] = x264_pixel_sad_4x4( fdec, 32, fenc, 16 );
    x264_predict_4x4_h_c( fdec );
    res[1] = x264_pixel_sad_4x4( fdec, 32, fenc, 16 );
    x264_predict_4x4_dc_c( fdec );
    res[2] = x264_pixel_sad_4x4( fdec, 32, fenc, 16 );
}
void x264_intra_satd_x3_4x4( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_4x4_v_c( fdec );
    res[0] = x264_pixel_satd_4x4( fdec, 32, fenc, 16 );
    x264_predict_4x4_h_c( fdec );
    res[1] = x264_pixel_satd_4x4( fdec, 32, fenc, 16 );
    x264_predict_4x4_dc_c( fdec );
    res[2] = x264_pixel_satd_4x4( fdec, 32, fenc, 16 );
}
void x264_intra_sad_x3_8x8c( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_8x8c_dc_c( fdec );
    res[0] = x264_pixel_sad_8x8( fdec, 32, fenc, 16 );
    x264_predict_8x8c_h_c( fdec );
    res[1] = x264_pixel_sad_8x8( fdec, 32, fenc, 16 );
    x264_predict_8x8c_v_c( fdec );
    res[2] = x264_pixel_sad_8x8( fdec, 32, fenc, 16 );
}
void x264_intra_satd_x3_8x8c( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_8x8c_dc_c( fdec );
    res[0] = x264_pixel_satd_8x8( fdec, 32, fenc, 16 );
    x264_predict_8x8c_h_c( fdec );
    res[1] = x264_pixel_satd_8x8( fdec, 32, fenc, 16 );
    x264_predict_8x8c_v_c( fdec );
    res[2] = x264_pixel_satd_8x8( fdec, 32, fenc, 16 );
}
void x264_intra_sad_x3_8x16c( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_8x16c_dc_c( fdec );
    res[0] = x264_pixel_sad_8x16( fdec, 32, fenc, 16 );
    x264_predict_8x16c_h_c( fdec );
    res[1] = x264_pixel_sad_8x16( fdec, 32, fenc, 16 );
    x264_predict_8x16c_v_c( fdec );
    res[2] = x264_pixel_sad_8x16( fdec, 32, fenc, 16 );
}
void x264_intra_satd_x3_8x16c( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_8x16c_dc_c( fdec );
    res[0] = x264_pixel_satd_8x16( fdec, 32, fenc, 16 );
    x264_predict_8x16c_h_c( fdec );
    res[1] = x264_pixel_satd_8x16( fdec, 32, fenc, 16 );
    x264_predict_8x16c_v_c( fdec );
    res[2] = x264_pixel_satd_8x16( fdec, 32, fenc, 16 );
}
void x264_intra_sad_x3_16x16( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_16x16_v_c( fdec );
    res[0] = x264_pixel_sad_16x16( fdec, 32, fenc, 16 );
    x264_predict_16x16_h_c( fdec );
    res[1] = x264_pixel_sad_16x16( fdec, 32, fenc, 16 );
    x264_predict_16x16_dc_c( fdec );
    res[2] = x264_pixel_sad_16x16( fdec, 32, fenc, 16 );
}
void x264_intra_satd_x3_16x16( pixel *fenc, pixel *fdec, int res[3] ) {
    x264_predict_16x16_v_c( fdec );
    res[0] = x264_pixel_satd_16x16( fdec, 32, fenc, 16 );
    x264_predict_16x16_h_c( fdec );
    res[1] = x264_pixel_satd_16x16( fdec, 32, fenc, 16 );
    x264_predict_16x16_dc_c( fdec );
    res[2] = x264_pixel_satd_16x16( fdec, 32, fenc, 16 );
}


// No C implementation of intra_satd_x9. See checkasm for its behavior,
// or see x264_mb_analyse_intra for the entirely different algorithm we
// use when lacking an asm implementation of it.



/****************************************************************************
 * structural similarity metric
 ****************************************************************************/
static void ssim_4x4x2_core( const pixel *pix1, intptr_t stride1,
                             const pixel *pix2, intptr_t stride2,
                             int sums[2][4] )
{
    for( int z = 0; z < 2; z++ )
    {
        uint32_t s1 = 0, s2 = 0, ss = 0, s12 = 0;
        for( int y = 0; y < 4; y++ )
            for( int x = 0; x < 4; x++ )
            {
                int a = pix1[x+y*stride1];
                int b = pix2[x+y*stride2];
                s1  += a;
                s2  += b;
                ss  += a*a;
                ss  += b*b;
                s12 += a*b;
            }
        sums[z][0] = s1;
        sums[z][1] = s2;
        sums[z][2] = ss;
        sums[z][3] = s12;
        pix1 += 4;
        pix2 += 4;
    }
}

static float ssim_end1( int s1, int s2, int ss, int s12 )
{
/* Maximum value for 10-bit is: ss*64 = (2^10-1)^2*16*4*64 = 4286582784, which will overflow in some cases.
 * s1*s1, s2*s2, and s1*s2 also obtain this value for edge cases: ((2^10-1)*16*4)^2 = 4286582784.
 * Maximum value for 9-bit is: ss*64 = (2^9-1)^2*16*4*64 = 1069551616, which will not overflow. */
 /*
#if BIT_DEPTH > 9
#define type float
    static const float ssim_c1 = .01*.01*PIXEL_MAX*PIXEL_MAX*64;
    static const float ssim_c2 = .03*.03*PIXEL_MAX*PIXEL_MAX*64*63;
#else
*/
#define type int
    static const int ssim_c1 = (int)(.01*.01*PIXEL_MAX*PIXEL_MAX*64 + .5);
    static const int ssim_c2 = (int)(.03*.03*PIXEL_MAX*PIXEL_MAX*64*63 + .5);
//#endif
    type fs1 = s1;
    type fs2 = s2;
    type fss = ss;
    type fs12 = s12;
    type vars = fss*64 - fs1*fs1 - fs2*fs2;
    type covar = fs12*64 - fs1*fs2;
    return (float)(2*fs1*fs2 + ssim_c1) * (float)(2*covar + ssim_c2)
         / ((float)(fs1*fs1 + fs2*fs2 + ssim_c1) * (float)(vars + ssim_c2));
#undef type
}

static float ssim_end4( int sum0[5][4], int sum1[5][4], int width )
{
    float ssim = 0.0;
    for( int i = 0; i < width; i++ )
        ssim += ssim_end1( sum0[i][0] + sum0[i+1][0] + sum1[i][0] + sum1[i+1][0],
                           sum0[i][1] + sum0[i+1][1] + sum1[i][1] + sum1[i+1][1],
                           sum0[i][2] + sum0[i+1][2] + sum1[i][2] + sum1[i+1][2],
                           sum0[i][3] + sum0[i+1][3] + sum1[i][3] + sum1[i+1][3] );
    return ssim;
}

float x264_pixel_ssim_wxh( x264_pixel_function_t *pf,
                           pixel *pix1, intptr_t stride1,
                           pixel *pix2, intptr_t stride2,
                           int width, int height, void *buf, int *cnt )
{
    int z = 0;
    float ssim = 0.0;
    int (*sum0)[4] = buf;
    int (*sum1)[4] = sum0 + (width >> 2) + 3;
    width >>= 2;
    height >>= 2;
    for( int y = 1; y < height; y++ )
    {
        for( ; z <= y; z++ )
        {
            XCHG( void*, sum0, sum1 );
            for( int x = 0; x < width; x+=2 )
                pf->ssim_4x4x2_core( &pix1[4*(x+z*stride1)], stride1, &pix2[4*(x+z*stride2)], stride2, &sum0[x] );
        }
        for( int x = 0; x < width-1; x += 4 )
            ssim += pf->ssim_end4( sum0+x, sum1+x, X264_MIN(4,width-x-1) );
    }
    *cnt = (height-1) * (width-1);
    return ssim;
}
//Absolute Sum Difference 8
//1. pix1 and pix2 differnce
//2. sum up in 8 byte width
//3. absolute
static int pixel_asd8( pixel *pix1, intptr_t stride1, pixel *pix2, intptr_t stride2, int height )
{
    int sum = 0;
    for( int y = 0; y < height; y++, pix1 += stride1, pix2 += stride2 )
        for( int x = 0; x < 8; x++ )
            sum += pix1[x] - pix2[x];
    return abs( sum );
}

/****************************************************************************
 * successive elimination
 ****************************************************************************/
static int x264_pixel_ads4( int enc_dc[4], uint16_t *sums, int delta,
                            uint16_t *cost_mvx, int16_t *mvs, int width, int thresh )
{
    int nmv = 0;
    for( int i = 0; i < width; i++, sums++ )
    {
        int ads = abs( enc_dc[0] - sums[0] )
                + abs( enc_dc[1] - sums[8] )
                + abs( enc_dc[2] - sums[delta] )
                + abs( enc_dc[3] - sums[delta+8] )
                + cost_mvx[i];
        if( ads < thresh )
            mvs[nmv++] = i;
    }
    return nmv;
}

static int x264_pixel_ads2( int enc_dc[2], uint16_t *sums, int delta,
                            uint16_t *cost_mvx, int16_t *mvs, int width, int thresh )
{
    int nmv = 0;
    for( int i = 0; i < width; i++, sums++ )
    {
        int ads = abs( enc_dc[0] - sums[0] )
                + abs( enc_dc[1] - sums[delta] )
                + cost_mvx[i];
        if( ads < thresh )
            mvs[nmv++] = i;
    }
    return nmv;
}

static int x264_pixel_ads1( int enc_dc[1], uint16_t *sums, int delta,
                            uint16_t *cost_mvx, int16_t *mvs, int width, int thresh )
{
    int nmv = 0;
    for( int i = 0; i<width; i++, sums++ )
    {
        int ads = abs( enc_dc[0] - sums[0] )
                + cost_mvx[i];
        if( ads < thresh )
            mvs[nmv++] = i;
    }
    return nmv;
}


/****************************************************************************
 * x264_pixel_init:
 ****************************************************************************/
void x264_pixel_init( int cpu, x264_pixel_function_t *pixf )
{
    memset( pixf, 0, sizeof(*pixf) );
/*
#define INIT2_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_16x16] = x264_pixel_##name2##_16x16##cpu;\
    pixf->name1[PIXEL_16x8]  = x264_pixel_##name2##_16x8##cpu;
#define INIT4_NAME( name1, name2, cpu ) \
    INIT2_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_8x16]  = x264_pixel_##name2##_8x16##cpu;\
    pixf->name1[PIXEL_8x8]   = x264_pixel_##name2##_8x8##cpu;
#define INIT5_NAME( name1, name2, cpu ) \
    INIT4_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_8x4]   = x264_pixel_##name2##_8x4##cpu;
#define INIT6_NAME( name1, name2, cpu ) \
    INIT5_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_4x8]   = x264_pixel_##name2##_4x8##cpu;
#define INIT7_NAME( name1, name2, cpu ) \
    INIT6_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_4x4]   = x264_pixel_##name2##_4x4##cpu;
#define INIT8_NAME( name1, name2, cpu ) \
    INIT7_NAME( name1, name2, cpu ) \
    pixf->name1[PIXEL_4x16]  = x264_pixel_##name2##_4x16##cpu;
#define INIT2( name, cpu ) INIT2_NAME( name, name, cpu )
#define INIT4( name, cpu ) INIT4_NAME( name, name, cpu )
#define INIT5( name, cpu ) INIT5_NAME( name, name, cpu )
#define INIT6( name, cpu ) INIT6_NAME( name, name, cpu )
#define INIT7( name, cpu ) INIT7_NAME( name, name, cpu )
#define INIT8( name, cpu ) INIT8_NAME( name, name, cpu )

#define INIT_ADS( cpu ) \
    pixf->ads[PIXEL_16x16] = x264_pixel_ads4##cpu;\
    pixf->ads[PIXEL_16x8] = x264_pixel_ads2##cpu;\
    pixf->ads[PIXEL_8x8] = x264_pixel_ads1##cpu;

    INIT8( sad, );
    INIT8_NAME( sad_aligned, sad, );
    INIT7( sad_x3, );
    INIT7( sad_x4, );
    INIT8( ssd, );
    INIT8( satd, );
    INIT7( satd_x3, );
    INIT7( satd_x4, );
    INIT4( hadamard_ac, );
    INIT_ADS( );
*/
    pixf->sad[PIXEL_16x16] = x264_pixel_sad_16x16;
    pixf->sad[PIXEL_16x8] = x264_pixel_sad_16x8;
    pixf->sad[PIXEL_8x16] = x264_pixel_sad_8x16;
    pixf->sad[PIXEL_8x8] = x264_pixel_sad_8x8;
    pixf->sad[PIXEL_8x4] = x264_pixel_sad_8x4;
    pixf->sad[PIXEL_4x8] = x264_pixel_sad_4x8;
    pixf->sad[PIXEL_4x4] = x264_pixel_sad_4x4;
    pixf->sad[PIXEL_4x16] = x264_pixel_sad_4x16;;
    pixf->sad_aligned[PIXEL_16x16] = x264_pixel_sad_16x16;
    pixf->sad_aligned[PIXEL_16x8] = x264_pixel_sad_16x8;
    pixf->sad_aligned[PIXEL_8x16] = x264_pixel_sad_8x16;
    pixf->sad_aligned[PIXEL_8x8] = x264_pixel_sad_8x8;
    pixf->sad_aligned[PIXEL_8x4] = x264_pixel_sad_8x4;
    pixf->sad_aligned[PIXEL_4x8] = x264_pixel_sad_4x8;
    pixf->sad_aligned[PIXEL_4x4] = x264_pixel_sad_4x4;
    pixf->sad_aligned[PIXEL_4x16] = x264_pixel_sad_4x16;;
    pixf->sad_x3[PIXEL_16x16] = x264_pixel_sad_x3_16x16;
    pixf->sad_x3[PIXEL_16x8] = x264_pixel_sad_x3_16x8;
    pixf->sad_x3[PIXEL_8x16] = x264_pixel_sad_x3_8x16;
    pixf->sad_x3[PIXEL_8x8] = x264_pixel_sad_x3_8x8;
    pixf->sad_x3[PIXEL_8x4] = x264_pixel_sad_x3_8x4;
    pixf->sad_x3[PIXEL_4x8] = x264_pixel_sad_x3_4x8;
    pixf->sad_x3[PIXEL_4x4] = x264_pixel_sad_x3_4x4;;
    pixf->sad_x4[PIXEL_16x16] = x264_pixel_sad_x4_16x16;
    pixf->sad_x4[PIXEL_16x8] = x264_pixel_sad_x4_16x8;
    pixf->sad_x4[PIXEL_8x16] = x264_pixel_sad_x4_8x16;
    pixf->sad_x4[PIXEL_8x8] = x264_pixel_sad_x4_8x8;
    pixf->sad_x4[PIXEL_8x4] = x264_pixel_sad_x4_8x4;
    pixf->sad_x4[PIXEL_4x8] = x264_pixel_sad_x4_4x8;
    pixf->sad_x4[PIXEL_4x4] = x264_pixel_sad_x4_4x4;;
    pixf->ssd[PIXEL_16x16] = x264_pixel_ssd_16x16;
    pixf->ssd[PIXEL_16x8] = x264_pixel_ssd_16x8;
    pixf->ssd[PIXEL_8x16] = x264_pixel_ssd_8x16;
    pixf->ssd[PIXEL_8x8] = x264_pixel_ssd_8x8;
    pixf->ssd[PIXEL_8x4] = x264_pixel_ssd_8x4;
    pixf->ssd[PIXEL_4x8] = x264_pixel_ssd_4x8;
    pixf->ssd[PIXEL_4x4] = x264_pixel_ssd_4x4;
    pixf->ssd[PIXEL_4x16] = x264_pixel_ssd_4x16;;
    pixf->satd[PIXEL_16x16] = x264_pixel_satd_16x16;
    pixf->satd[PIXEL_16x8] = x264_pixel_satd_16x8;
    pixf->satd[PIXEL_8x16] = x264_pixel_satd_8x16;
    pixf->satd[PIXEL_8x8] = x264_pixel_satd_8x8;
    pixf->satd[PIXEL_8x4] = x264_pixel_satd_8x4;
    pixf->satd[PIXEL_4x8] = x264_pixel_satd_4x8;
    pixf->satd[PIXEL_4x4] = x264_pixel_satd_4x4;
    pixf->satd[PIXEL_4x16] = x264_pixel_satd_4x16;;
    pixf->satd_x3[PIXEL_16x16] = x264_pixel_satd_x3_16x16;
    pixf->satd_x3[PIXEL_16x8] = x264_pixel_satd_x3_16x8;
    pixf->satd_x3[PIXEL_8x16] = x264_pixel_satd_x3_8x16;
    pixf->satd_x3[PIXEL_8x8] = x264_pixel_satd_x3_8x8;
    pixf->satd_x3[PIXEL_8x4] = x264_pixel_satd_x3_8x4;
    pixf->satd_x3[PIXEL_4x8] = x264_pixel_satd_x3_4x8;
    pixf->satd_x3[PIXEL_4x4] = x264_pixel_satd_x3_4x4;;
    pixf->satd_x4[PIXEL_16x16] = x264_pixel_satd_x4_16x16;
    pixf->satd_x4[PIXEL_16x8] = x264_pixel_satd_x4_16x8;
    pixf->satd_x4[PIXEL_8x16] = x264_pixel_satd_x4_8x16;
    pixf->satd_x4[PIXEL_8x8] = x264_pixel_satd_x4_8x8;
    pixf->satd_x4[PIXEL_8x4] = x264_pixel_satd_x4_8x4;
    pixf->satd_x4[PIXEL_4x8] = x264_pixel_satd_x4_4x8;
    pixf->satd_x4[PIXEL_4x4] = x264_pixel_satd_x4_4x4;;
    pixf->hadamard_ac[PIXEL_16x16] = x264_pixel_hadamard_ac_16x16;
    pixf->hadamard_ac[PIXEL_16x8] = x264_pixel_hadamard_ac_16x8;
    pixf->hadamard_ac[PIXEL_8x16] = x264_pixel_hadamard_ac_8x16;
    pixf->hadamard_ac[PIXEL_8x8] = x264_pixel_hadamard_ac_8x8;;
    pixf->ads[PIXEL_16x16] = x264_pixel_ads4;
    pixf->ads[PIXEL_16x8] = x264_pixel_ads2;
    pixf->ads[PIXEL_8x8] = x264_pixel_ads1;;

    pixf->sa8d[PIXEL_16x16] = x264_pixel_sa8d_16x16;
    pixf->sa8d[PIXEL_8x8]   = x264_pixel_sa8d_8x8;
    pixf->var[PIXEL_16x16] = x264_pixel_var_16x16;
    pixf->var[PIXEL_8x16]  = x264_pixel_var_8x16;
    pixf->var[PIXEL_8x8]   = x264_pixel_var_8x8;
    pixf->var2[PIXEL_8x16]  = x264_pixel_var2_8x16;
    pixf->var2[PIXEL_8x8]   = x264_pixel_var2_8x8;

    pixf->ssd_nv12_core = pixel_ssd_nv12_core;
    pixf->ssim_4x4x2_core = ssim_4x4x2_core;
    pixf->ssim_end4 = ssim_end4;
    pixf->asd8 = pixel_asd8;

    pixf->intra_sad_x3_4x4    = x264_intra_sad_x3_4x4;
    pixf->intra_satd_x3_4x4   = x264_intra_satd_x3_4x4;
    pixf->intra_sad_x3_8x8    = x264_intra_sad_x3_8x8;
    pixf->intra_sa8d_x3_8x8   = x264_intra_sa8d_x3_8x8;
    pixf->intra_sad_x3_8x8c   = x264_intra_sad_x3_8x8c;
    pixf->intra_satd_x3_8x8c  = x264_intra_satd_x3_8x8c;
    pixf->intra_sad_x3_8x16c  = x264_intra_sad_x3_8x16c;
    pixf->intra_satd_x3_8x16c = x264_intra_satd_x3_8x16c;
    pixf->intra_sad_x3_16x16  = x264_intra_sad_x3_16x16;
    pixf->intra_satd_x3_16x16 = x264_intra_satd_x3_16x16;

    pixf->ads[PIXEL_8x16] =
    pixf->ads[PIXEL_8x4] =
    pixf->ads[PIXEL_4x8] = pixf->ads[PIXEL_16x8];
    pixf->ads[PIXEL_4x4] = pixf->ads[PIXEL_8x8];
}

