/*****************************************************************************
 * vlc.c : vlc tables
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
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

#include "common.h"

/* [nC] */
const vlc_t x264_coeff0_token[6] =
{
    { 0x1, 1 }, /* str=1 */
    { 0x3, 2 }, /* str=11 */
    { 0xf, 4 }, /* str=1111 */
    { 0x3, 6 }, /* str=000011 */
    { 0x1, 2 }, /* str=01 */
    { 0x1, 1 }, /* str=1 */
};

/* [nC][i_total_coeff-1][i_trailing] */
const vlc_t x264_coeff_token[6][16][4] =
{
    { /* table 0 */
        { /* i_total 1 */
            { 0x5, 6 }, /* str=000101 */
            { 0x1, 2 }, /* str=01 */
        },
        { /* i_total 2 */
            { 0x7, 8 }, /* str=00000111 */
            { 0x4, 6 }, /* str=000100 */
            { 0x1, 3 }, /* str=001 */
        },
        { /* i_total 3 */
            { 0x7, 9 }, /* str=000000111 */
            { 0x6, 8 }, /* str=00000110 */
            { 0x5, 7 }, /* str=0000101 */
            { 0x3, 5 }, /* str=00011 */
        },
        { /* i_total 4 */
            { 0x7, 10 }, /* str=0000000111 */
            { 0x6, 9 },  /* str=000000110 */
            { 0x5, 8 },  /* str=00000101 */
            { 0x3, 6 },  /* str=000011 */
        },
        { /* i_total 5 */
            { 0x7, 11 }, /* str=00000000111 */
            { 0x6, 10 }, /* str=0000000110 */
            { 0x5, 9 },  /* str=000000101 */
            { 0x4, 7 },  /* str=0000100 */
        },
        { /* i_total 6 */
            { 0xf, 13 }, /* str=0000000001111 */
            { 0x6, 11 }, /* str=00000000110 */
            { 0x5, 10 }, /* str=0000000101 */
            { 0x4, 8 },  /* str=00000100 */
        },
        { /* i_total 7 */
            { 0xb, 13 }, /* str=0000000001011 */
            { 0xe, 13 }, /* str=0000000001110 */
            { 0x5, 11 }, /* str=00000000101 */
            { 0x4, 9 },  /* str=000000100 */
        },
        { /* i_total 8 */
            { 0x8, 13 }, /* str=0000000001000 */
            { 0xa, 13 }, /* str=0000000001010 */
            { 0xd, 13 }, /* str=0000000001101 */
            { 0x4, 10 }, /* str=0000000100 */
        },
        { /* i_total 9 */
            { 0xf, 14 }, /* str=00000000001111 */
            { 0xe, 14 }, /* str=00000000001110 */
            { 0x9, 13 }, /* str=0000000001001 */
            { 0x4, 11 }, /* str=00000000100 */
        },
        { /* i_total 10 */
            { 0xb, 14 }, /* str=00000000001011 */
            { 0xa, 14 }, /* str=00000000001010 */
            { 0xd, 14 }, /* str=00000000001101 */
            { 0xc, 13 }, /* str=0000000001100 */
        },
        { /* i_total 14 */
            { 0xf, 15 }, /* str=000000000001111 */
            { 0xe, 15 }, /* str=000000000001110 */
            { 0x9, 14 }, /* str=00000000001001 */
            { 0xc, 14 }, /* str=00000000001100 */
        },
        { /* i_total 12 */
            { 0xb, 15 }, /* str=000000000001011 */
            { 0xa, 15 }, /* str=000000000001010 */
            { 0xd, 15 }, /* str=000000000001101 */
            { 0x8, 14 }, /* str=00000000001000 */
        },
        { /* i_total 13 */
            { 0xf, 16 }, /* str=0000000000001111 */
            { 0x1, 15 }, /* str=000000000000001 */
            { 0x9, 15 }, /* str=000000000001001 */
            { 0xc, 15 }, /* str=000000000001100 */
        },
        { /* i_total 14 */
            { 0xb, 16 }, /* str=0000000000001011 */
            { 0xe, 16 }, /* str=0000000000001110 */
            { 0xd, 16 }, /* str=0000000000001101 */
            { 0x8, 15 }, /* str=000000000001000 */
        },
        { /* i_total 15 */
            { 0x7, 16 }, /* str=0000000000000111 */
            { 0xa, 16 }, /* str=0000000000001010 */
            { 0x9, 16 }, /* str=0000000000001001 */
            { 0xc, 16 }, /* str=0000000000001100 */
        },
        { /* i_total 16 */
            { 0x4, 16 }, /* str=0000000000000100 */
            { 0x6, 16 }, /* str=0000000000000110 */
            { 0x5, 16 }, /* str=0000000000000101 */
            { 0x8, 16 }, /* str=0000000000001000 */
        },
    },
    { /* table 1 */
        { /* i_total 1 */
            { 0xb, 6 }, /* str=001011 */
            { 0x2, 2 }, /* str=10 */
        },
        { /* i_total 2 */
            { 0x7, 6 }, /* str=000111 */
            { 0x7, 5 }, /* str=00111 */
            { 0x3, 3 }, /* str=011 */
        },
        { /* i_total 3 */
            { 0x7, 7 }, /* str=0000111 */
            { 0xa, 6 }, /* str=001010 */
            { 0x9, 6 }, /* str=001001 */
            { 0x5, 4 }, /* str=0101 */
        },
        { /* i_total 4 */
            { 0x7, 8 }, /* str=00000111 */
            { 0x6, 6 }, /* str=000110 */
            { 0x5, 6 }, /* str=000101 */
            { 0x4, 4 }, /* str=0100 */
        },
        { /* i_total 5 */
            { 0x4, 8 }, /* str=00000100 */
            { 0x6, 7 }, /* str=0000110 */
            { 0x5, 7 }, /* str=0000101 */
            { 0x6, 5 }, /* str=00110 */
        },
        { /* i_total 6 */
            { 0x7, 9 }, /* str=000000111 */
            { 0x6, 8 }, /* str=00000110 */
            { 0x5, 8 }, /* str=00000101 */
            { 0x8, 6 }, /* str=001000 */
        },
        { /* i_total 7 */
            { 0xf, 11 }, /* str=00000001111 */
            { 0x6, 9 },  /* str=000000110 */
            { 0x5, 9 },  /* str=000000101 */
            { 0x4, 6 },  /* str=000100 */
        },
        { /* i_total 8 */
            { 0xb, 11 }, /* str=00000001011 */
            { 0xe, 11 }, /* str=00000001110 */
            { 0xd, 11 }, /* str=00000001101 */
            { 0x4, 7 },  /* str=0000100 */
        },
        { /* i_total 9 */
            { 0xf, 12 }, /* str=000000001111 */
            { 0xa, 11 }, /* str=00000001010 */
            { 0x9, 11 }, /* str=00000001001 */
            { 0x4, 9 },  /* str=000000100 */
        },
        { /* i_total 10 */
            { 0xb, 12 }, /* str=000000001011 */
            { 0xe, 12 }, /* str=000000001110 */
            { 0xd, 12 }, /* str=000000001101 */
            { 0xc, 11 }, /* str=00000001100 */
        },
        { /* i_total 11 */
            { 0x8, 12 }, /* str=000000001000 */
            { 0xa, 12 }, /* str=000000001010 */
            { 0x9, 12 }, /* str=000000001001 */
            { 0x8, 11 }, /* str=00000001000 */
        },
        { /* i_total 12 */
            { 0xf, 13 }, /* str=0000000001111 */
            { 0xe, 13 }, /* str=0000000001110 */
            { 0xd, 13 }, /* str=0000000001101 */
            { 0xc, 12 }, /* str=000000001100 */
        },
        { /* i_total 13 */
            { 0xb, 13 }, /* str=0000000001011 */
            { 0xa, 13 }, /* str=0000000001010 */
            { 0x9, 13 }, /* str=0000000001001 */
            { 0xc, 13 }, /* str=0000000001100 */
        },
        { /* i_total 14 */
            { 0x7, 13 }, /* str=0000000000111 */
            { 0xb, 14 }, /* str=00000000001011 */
            { 0x6, 13 }, /* str=0000000000110 */
            { 0x8, 13 }, /* str=0000000001000 */
        },
        { /* i_total 15 */
            { 0x9, 14 }, /* str=00000000001001 */
            { 0x8, 14 }, /* str=00000000001000 */
            { 0xa, 14 }, /* str=00000000001010 */
            { 0x1, 13 }, /* str=0000000000001 */
        },
        { /* i_total 16 */
            { 0x7, 14 }, /* str=00000000000111 */
            { 0x6, 14 }, /* str=00000000000110 */
            { 0x5, 14 }, /* str=00000000000101 */
            { 0x4, 14 }, /* str=00000000000100 */
        },
    },
    { /* table 2 */
        { /* i_total 1 */
            { 0xf, 6 }, /* str=001111 */
            { 0xe, 4 }, /* str=1110 */
        },
        { /* i_total 2 */
            { 0xb, 6 }, /* str=001011 */
            { 0xf, 5 }, /* str=01111 */
            { 0xd, 4 }, /* str=1101 */
        },
        { /* i_total 3 */
            { 0x8, 6 }, /* str=001000 */
            { 0xc, 5 }, /* str=01100 */
            { 0xe, 5 }, /* str=01110 */
            { 0xc, 4 }, /* str=1100 */
        },
        { /* i_total 4 */
            { 0xf, 7 }, /* str=0001111 */
            { 0xa, 5 }, /* str=01010 */
            { 0xb, 5 }, /* str=01011 */
            { 0xb, 4 }, /* str=1011 */
        },
        { /* i_total 5 */
            { 0xb, 7 }, /* str=0001011 */
            { 0x8, 5 }, /* str=01000 */
            { 0x9, 5 }, /* str=01001 */
            { 0xa, 4 }, /* str=1010 */
        },
        { /* i_total 6 */
            { 0x9, 7 }, /* str=0001001 */
            { 0xe, 6 }, /* str=001110 */
            { 0xd, 6 }, /* str=001101 */
            { 0x9, 4 }, /* str=1001 */
        },
        { /* i_total 7 */
            { 0x8, 7 }, /* str=0001000 */
            { 0xa, 6 }, /* str=001010 */
            { 0x9, 6 }, /* str=001001 */
            { 0x8, 4 }, /* str=1000 */
        },
        { /* i_total 8 */
            { 0xf, 8 }, /* str=00001111 */
            { 0xe, 7 }, /* str=0001110 */
            { 0xd, 7 }, /* str=0001101 */
            { 0xd, 5 }, /* str=01101 */
        },
        { /* i_total 9 */
            { 0xb, 8 }, /* str=00001011 */
            { 0xe, 8 }, /* str=00001110 */
            { 0xa, 7 }, /* str=0001010 */
            { 0xc, 6 }, /* str=001100 */
        },
        { /* i_total 10 */
            { 0xf, 9 }, /* str=000001111 */
            { 0xa, 8 }, /* str=00001010 */
            { 0xd, 8 }, /* str=00001101 */
            { 0xc, 7 }, /* str=0001100 */
        },
        { /* i_total 11 */
            { 0xb, 9 }, /* str=000001011 */
            { 0xe, 9 }, /* str=000001110 */
            { 0x9, 8 }, /* str=00001001 */
            { 0xc, 8 }, /* str=00001100 */
        },
        { /* i_total 12 */
            { 0x8, 9 }, /* str=000001000 */
            { 0xa, 9 }, /* str=000001010 */
            { 0xd, 9 }, /* str=000001101 */
            { 0x8, 8 }, /* str=00001000 */
        },
        { /* i_total 13 */
            { 0xd, 10 }, /* str=0000001101 */
            { 0x7, 9 },  /* str=000000111 */
            { 0x9, 9 },  /* str=000001001 */
            { 0xc, 9 },  /* str=000001100 */
        },
        { /* i_total 14 */
            { 0x9, 10 }, /* str=0000001001 */
            { 0xc, 10 }, /* str=0000001100 */
            { 0xb, 10 }, /* str=0000001011 */
            { 0xa, 10 }, /* str=0000001010 */
        },
        { /* i_total 15 */
            { 0x5, 10 }, /* str=0000000101 */
            { 0x8, 10 }, /* str=0000001000 */
            { 0x7, 10 }, /* str=0000000111 */
            { 0x6, 10 }, /* str=0000000110 */
        },
        { /* i_total 16 */
            { 0x1, 10 }, /* str=0000000001 */
            { 0x4, 10 }, /* str=0000000100 */
            { 0x3, 10 }, /* str=0000000011 */
            { 0x2, 10 }, /* str=0000000010 */
        },
    },
    { /* table 3 */
        { /* i_total 1 */
            { 0x0, 6 }, /* str=000000 */
            { 0x1, 6 }, /* str=000001 */
        },
        { /* i_total 2 */
            { 0x4, 6 }, /* str=000100 */
            { 0x5, 6 }, /* str=000101 */
            { 0x6, 6 }, /* str=000110 */
        },
        { /* i_total 3 */
            { 0x8, 6 }, /* str=001000 */
            { 0x9, 6 }, /* str=001001 */
            { 0xa, 6 }, /* str=001010 */
            { 0xb, 6 }, /* str=001011 */
        },
        { /* i_total 4 */
            { 0xc, 6 }, /* str=001100 */
            { 0xd, 6 }, /* str=001101 */
            { 0xe, 6 }, /* str=001110 */
            { 0xf, 6 }, /* str=001111 */
        },
        { /* i_total 5 */
            { 0x10, 6 }, /* str=010000 */
            { 0x11, 6 }, /* str=010001 */
            { 0x12, 6 }, /* str=010010 */
            { 0x13, 6 }, /* str=010011 */
        },
        { /* i_total 6 */
            { 0x14, 6 }, /* str=010100 */
            { 0x15, 6 }, /* str=010101 */
            { 0x16, 6 }, /* str=010110 */
            { 0x17, 6 }, /* str=010111 */
        },
        { /* i_total 7 */
            { 0x18, 6 }, /* str=011000 */
            { 0x19, 6 }, /* str=011001 */
            { 0x1a, 6 }, /* str=011010 */
            { 0x1b, 6 }, /* str=011011 */
        },
        { /* i_total 8 */
            { 0x1c, 6 }, /* str=011100 */
            { 0x1d, 6 }, /* str=011101 */
            { 0x1e, 6 }, /* str=011110 */
            { 0x1f, 6 }, /* str=011111 */
        },
        { /* i_total 9 */
            { 0x20, 6 }, /* str=100000 */
            { 0x21, 6 }, /* str=100001 */
            { 0x22, 6 }, /* str=100010 */
            { 0x23, 6 }, /* str=100011 */
        },
        { /* i_total 10 */
            { 0x24, 6 }, /* str=100100 */
            { 0x25, 6 }, /* str=100101 */
            { 0x26, 6 }, /* str=100110 */
            { 0x27, 6 }, /* str=100111 */
        },
        { /* i_total 11 */
            { 0x28, 6 }, /* str=101000 */
            { 0x29, 6 }, /* str=101001 */
            { 0x2a, 6 }, /* str=101010 */
            { 0x2b, 6 }, /* str=101011 */
        },
        { /* i_total 12 */
            { 0x2c, 6 }, /* str=101100 */
            { 0x2d, 6 }, /* str=101101 */
            { 0x2e, 6 }, /* str=101110 */
            { 0x2f, 6 }, /* str=101111 */
        },
        { /* i_total 13 */
            { 0x30, 6 }, /* str=110000 */
            { 0x31, 6 }, /* str=110001 */
            { 0x32, 6 }, /* str=110010 */
            { 0x33, 6 }, /* str=110011 */
        },
        { /* i_total 14 */
            { 0x34, 6 }, /* str=110100 */
            { 0x35, 6 }, /* str=110101 */
            { 0x36, 6 }, /* str=110110 */
            { 0x37, 6 }, /* str=110111 */
        },
        { /* i_total 15 */
            { 0x38, 6 }, /* str=111000 */
            { 0x39, 6 }, /* str=111001 */
            { 0x3a, 6 }, /* str=111010 */
            { 0x3b, 6 }, /* str=111011 */
        },
        { /* i_total 16 */
            { 0x3c, 6 }, /* str=111100 */
            { 0x3d, 6 }, /* str=111101 */
            { 0x3e, 6 }, /* str=111110 */
            { 0x3f, 6 }, /* str=111111 */
        },
    },
    { /* table 4 */
        { /* i_total 1 */
            { 0x7, 6 }, /* str=000111 */
            { 0x1, 1 }, /* str=1 */
        },
        { /* i_total 2 */
            { 0x4, 6 }, /* str=000100 */
            { 0x6, 6 }, /* str=000110 */
            { 0x1, 3 }, /* str=001 */
        },
        { /* i_total 3 */
            { 0x3, 6 }, /* str=000011 */
            { 0x3, 7 }, /* str=0000011 */
            { 0x2, 7 }, /* str=0000010 */
            { 0x5, 6 }, /* str=000101 */
        },
        { /* i_total 4 */
            { 0x2, 6 }, /* str=000010 */
            { 0x3, 8 }, /* str=00000011 */
            { 0x2, 8 }, /* str=00000010 */
            { 0x0, 7 }, /* str=0000000 */
        },
    },
    { /* table 5 */
        { /* i_total 1 */
            { 0xf, 7 }, /* str=0001111 */
            { 0x1, 2 }, /* str=01 */
        },
        { /* i_total 2 */
            { 0xe, 7 }, /* str=0001110 */
            { 0xd, 7 }, /* str=0001101 */
            { 0x1, 3 }, /* str=001 */
        },
        { /* i_total 3 */
            { 0x7, 9 }, /* str=000000111 */
            { 0xc, 7 }, /* str=0001100 */
            { 0xb, 7 }, /* str=0001011 */
            { 0x1, 5 }, /* str=00001 */
        },
        { /* i_total 4 */
            { 0x6, 9 }, /* str=000000110 */
            { 0x5, 9 }, /* str=000000101 */
            { 0xa, 7 }, /* str=0001010 */
            { 0x1, 6 }, /* str=000001 */
        },
        { /* i_total 5 */
            { 0x7, 10 }, /* str=0000000111 */
            { 0x6, 10 }, /* str=0000000110 */
            { 0x4, 9 },  /* str=000000100 */
            { 0x9, 7 },  /* str=0001001 */
        },
        { /* i_total 6 */
            { 0x7, 11 }, /* str=00000000111 */
            { 0x6, 11 }, /* str=00000000110 */
            { 0x5, 10 }, /* str=0000000101 */
            { 0x8, 7 },  /* str=0001000 */
        },
        { /* i_total 7 */
            { 0x7, 12 }, /* str=000000000111 */
            { 0x6, 12 }, /* str=000000000110 */
            { 0x5, 11 }, /* str=00000000101 */
            { 0x4, 10 }, /* str=0000000100 */
        },
        { /* i_total 8 */
            { 0x7, 13 }, /* str=0000000000111 */
            { 0x5, 12 }, /* str=000000000101 */
            { 0x4, 12 }, /* str=000000000100 */
            { 0x4, 11 }, /* str=00000000100 */
        },
    },
};

/* [i_total_coeff-1][i_total_zeros] */
const vlc_t x264_total_zeros[15][16] =
{
    { /* i_total 1 */
        { 0x1, 1 }, /* str=1 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 3 }, /* str=010 */
        { 0x3, 4 }, /* str=0011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x3, 5 }, /* str=00011 */
        { 0x2, 5 }, /* str=00010 */
        { 0x3, 6 }, /* str=000011 */
        { 0x2, 6 }, /* str=000010 */
        { 0x3, 7 }, /* str=0000011 */
        { 0x2, 7 }, /* str=0000010 */
        { 0x3, 8 }, /* str=00000011 */
        { 0x2, 8 }, /* str=00000010 */
        { 0x3, 9 }, /* str=000000011 */
        { 0x2, 9 }, /* str=000000010 */
        { 0x1, 9 }, /* str=000000001 */
    },
    { /* i_total 2 */
        { 0x7, 3 }, /* str=111 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x5, 4 }, /* str=0101 */
        { 0x4, 4 }, /* str=0100 */
        { 0x3, 4 }, /* str=0011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x3, 5 }, /* str=00011 */
        { 0x2, 5 }, /* str=00010 */
        { 0x3, 6 }, /* str=000011 */
        { 0x2, 6 }, /* str=000010 */
        { 0x1, 6 }, /* str=000001 */
        { 0x0, 6 }, /* str=000000 */
    },
    { /* i_total 3 */
        { 0x5, 4 }, /* str=0101 */
        { 0x7, 3 }, /* str=111 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 4 }, /* str=0100 */
        { 0x3, 4 }, /* str=0011 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x3, 5 }, /* str=00011 */
        { 0x2, 5 }, /* str=00010 */
        { 0x1, 6 }, /* str=000001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x0, 6 }, /* str=000000 */
    },
    { /* i_total 4 */
        { 0x3, 5 }, /* str=00011 */
        { 0x7, 3 }, /* str=111 */
        { 0x5, 4 }, /* str=0101 */
        { 0x4, 4 }, /* str=0100 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 4 }, /* str=0011 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x2, 5 }, /* str=00010 */
        { 0x1, 5 }, /* str=00001 */
        { 0x0, 5 }, /* str=00000 */
    },
    { /* i_total 5 */
        { 0x5, 4 }, /* str=0101 */
        { 0x4, 4 }, /* str=0100 */
        { 0x3, 4 }, /* str=0011 */
        { 0x7, 3 }, /* str=111 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x1, 5 }, /* str=00001 */
        { 0x1, 4 }, /* str=0001 */
        { 0x0, 5 }, /* str=00000 */
    },
    { /* i_total 6 */
        { 0x1, 6 }, /* str=000001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x7, 3 }, /* str=111 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 6 }, /* str=000000 */
    },
    { /* i_total 7 */
        { 0x1, 6 }, /* str=000001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 6 }, /* str=000000 */
    },
    { /* i_total 8 */
        { 0x1, 6 }, /* str=000001 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x3, 3 }, /* str=011 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 6 }, /* str=000000 */
    },
    { /* i_total 9 */
        { 0x1, 6 }, /* str=000001 */
        { 0x0, 6 }, /* str=000000 */
        { 0x1, 4 }, /* str=0001 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x1, 3 }, /* str=001 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 5 }, /* str=00001 */
    },
    { /* i_total 10 */
        { 0x1, 5 }, /* str=00001 */
        { 0x0, 5 }, /* str=00000 */
        { 0x1, 3 }, /* str=001 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 4 }, /* str=0001 */
    },
    { /* i_total 11 */
        { 0x0, 4 }, /* str=0000 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 3 }, /* str=001 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 1 }, /* str=1 */
        { 0x3, 3 }, /* str=011 */
    },
    { /* i_total 12 */
        { 0x0, 4 }, /* str=0000 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 1 }, /* str=1 */
        { 0x1, 3 }, /* str=001 */
    },
    { /* i_total 13 */
        { 0x0, 3 }, /* str=000 */
        { 0x1, 3 }, /* str=001 */
        { 0x1, 1 }, /* str=1 */
        { 0x1, 2 }, /* str=01 */
    },
    { /* i_total 14 */
        { 0x0, 2 }, /* str=00 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 1 }, /* str=1 */
    },
    { /* i_total 15 */
        { 0x0, 1 }, /* str=0 */
        { 0x1, 1 }, /* str=1 */
    },
};

/* [i_total_coeff-1][i_total_zeros] */
const vlc_t x264_total_zeros_2x2_dc[3][4] =
{
    { /* i_total 1 */
        { 0x1, 1 }, /* str=1 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 3 }  /* str=000 */
    },
    { /* i_total 2 */
        { 0x1, 1 }, /* str=1 */
        { 0x1, 2 }, /* str=01 */
        { 0x0, 2 }, /* str=00 */
    },
    { /* i_total 3 */
        { 0x1, 1 }, /* str=1 */
        { 0x0, 1 }, /* str=0 */
    },
};

/* [i_total_coeff-1][i_total_zeros] */
const vlc_t x264_total_zeros_2x4_dc[7][8] =
{
    { /* i_total 1 */
        { 0x1, 1 }, /* str=1 */
        { 0x2, 3 }, /* str=010 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 4 }, /* str=0010 */
        { 0x3, 4 }, /* str=0011 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x0, 5 }, /* str=00000 */
    },
    { /* i_total 2 */
        { 0x0, 3 }, /* str=000 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 3 }, /* str=001 */
        { 0x4, 3 }, /* str=100 */
        { 0x5, 3 }, /* str=101 */
        { 0x6, 3 }, /* str=110 */
        { 0x7, 3 }, /* str=111 */
    },
    { /* i_total 3 */
        { 0x0, 3 }, /* str=000 */
        { 0x1, 3 }, /* str=001 */
        { 0x1, 2 }, /* str=01 */
        { 0x2, 2 }, /* str=10 */
        { 0x6, 3 }, /* str=110 */
        { 0x7, 3 }, /* str=111 */
    },
    { /* i_total 4 */
        { 0x6, 3 }, /* str=110 */
        { 0x0, 2 }, /* str=00 */
        { 0x1, 2 }, /* str=01 */
        { 0x2, 2 }, /* str=10 */
        { 0x7, 3 }, /* str=111 */
    },
    { /* i_total 5 */
        { 0x0, 2 }, /* str=00 */
        { 0x1, 2 }, /* str=01 */
        { 0x2, 2 }, /* str=10 */
        { 0x3, 2 }, /* str=11 */
    },
    { /* i_total 6 */
        { 0x0, 2 }, /* str=00 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 1 }, /* str=1 */
    },
    { /* i_total 7 */
        { 0x0, 1 }, /* str=0 */
        { 0x1, 1 }, /* str=1 */
    }
};

/* [MIN( i_zero_left-1, 6 )][run_before] */
static const vlc_t run_before[7][16] =
{
    { /* i_zero_left 1 */
        { 0x1, 1 }, /* str=1 */
        { 0x0, 1 }, /* str=0 */
    },
    { /* i_zero_left 2 */
        { 0x1, 1 }, /* str=1 */
        { 0x1, 2 }, /* str=01 */
        { 0x0, 2 }, /* str=00 */
    },
    { /* i_zero_left 3 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x1, 2 }, /* str=01 */
        { 0x0, 2 }, /* str=00 */
    },
    { /* i_zero_left 4 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x1, 2 }, /* str=01 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 3 }, /* str=000 */
    },
    { /* i_zero_left 5 */
        { 0x3, 2 }, /* str=11 */
        { 0x2, 2 }, /* str=10 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 3 }, /* str=001 */
        { 0x0, 3 }, /* str=000 */
    },
    { /* i_zero_left 6 */
        { 0x3, 2 }, /* str=11 */
        { 0x0, 3 }, /* str=000 */
        { 0x1, 3 }, /* str=001 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 3 }, /* str=010 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
    },
    { /* i_zero_left >6 */
        { 0x7, 3 }, /* str=111 */
        { 0x6, 3 }, /* str=110 */
        { 0x5, 3 }, /* str=101 */
        { 0x4, 3 }, /* str=100 */
        { 0x3, 3 }, /* str=011 */
        { 0x2, 3 }, /* str=010 */
        { 0x1, 3 }, /* str=001 */
        { 0x1, 4 }, /* str=0001 */
        { 0x1, 5 }, /* str=00001 */
        { 0x1, 6 }, /* str=000001 */
        { 0x1, 7 }, /* str=0000001 */
        { 0x1, 8 }, /* str=00000001 */
        { 0x1, 9 }, /* str=000000001 */
        { 0x1, 10 }, /* str=0000000001 */
        { 0x1, 11 }, /* str=00000000001 */
    },
};

vlc_large_t x264_level_token[7][LEVEL_TABLE_SIZE];
uint32_t x264_run_before[1<<16];

void x264_cavlc_init( x264_t *h )
{
    for( int i_suffix = 0; i_suffix < 7; i_suffix++ )
        for( int16_t level = -LEVEL_TABLE_SIZE/2; level < LEVEL_TABLE_SIZE/2; level++ )
        {
            int mask = level >> 15;
            int abs_level = (level^mask)-mask;
            int i_level_code = abs_level*2-mask-2;
            int i_next = i_suffix;
            vlc_large_t *vlc = &x264_level_token[i_suffix][level+LEVEL_TABLE_SIZE/2];

            if( ( i_level_code >> i_suffix ) < 14 )
            {
                vlc->i_size = (i_level_code >> i_suffix) + 1 + i_suffix;
                vlc->i_bits = (1<<i_suffix) + (i_level_code & ((1<<i_suffix)-1));
            }
            else if( i_suffix == 0 && i_level_code < 30 )
            {
                vlc->i_size = 19;
                vlc->i_bits = (1<<4) + (i_level_code - 14);
            }
            else if( i_suffix > 0 && ( i_level_code >> i_suffix ) == 14 )
            {
                vlc->i_size = 15 + i_suffix;
                vlc->i_bits = (1<<i_suffix) + (i_level_code & ((1<<i_suffix)-1));
            }
            else
            {
                i_level_code -= 15 << i_suffix;
                if( i_suffix == 0 )
                    i_level_code -= 15;
                vlc->i_size = 28;
                vlc->i_bits = (1<<12) + i_level_code;
            }
            if( i_next == 0 )
                i_next++;
            if( abs_level > (3 << (i_next-1)) && i_next < 6 )
                i_next++;
            vlc->i_next = i_next;
        }

    for( int i = 1; i < (1<<16); i++ )
    {
        x264_run_level_t runlevel;
        ALIGNED_ARRAY_16( dctcoef, dct, [16] );
        int size = 0;
        int bits = 0;
        for( int j = 0; j < 16; j++ )
            dct[j] = i&(1<<j);
        int total = h->quantf.coeff_level_run[DCT_LUMA_4x4]( dct, &runlevel );
        int zeros = runlevel.last + 1 - total;
        uint32_t mask = i << (x264_clz( i ) + 1);
        for( int j = 0; j < total-1 && zeros > 0; j++ )
        {
            int idx = X264_MIN(zeros, 7) - 1;
            int run = x264_clz( mask );
            int len = run_before[idx][run].i_size;
            size += len;
            bits <<= len;
            bits |= run_before[idx][run].i_bits;
            zeros -= run;
            mask <<= run + 1;
        }
        x264_run_before[i] = (bits << 5) + size;
    }
}
