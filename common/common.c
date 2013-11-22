/*****************************************************************************
 * common.c: misc common functions
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

#include "common.h"

#include <stdarg.h>
#include <ctype.h>

#if HAVE_MALLOC_H
#include <malloc.h>
#endif

const int x264_bit_depth = BIT_DEPTH;

const int x264_chroma_format = X264_CSP_I420;

static void x264_log_default( void *, int, const char *, va_list );

/****************************************************************************
 * x264_param_default:
 ****************************************************************************/
void x264_param_default( x264_param_t *param )
{
    /* */
    memset( param, 0, sizeof( x264_param_t ) );

    /* CPU autodetect */
    param->cpu = 0;

    /* Video properties */
    param->i_csp           = X264_CSP_I420;
    param->i_width         = 0;
    param->i_height        = 0;
    param->vui.i_sar_width = 0;
    param->vui.i_sar_height= 0;
    param->vui.i_overscan  = 0;  /* undef */
    param->vui.i_vidformat = 5;  /* undef */
    param->vui.b_fullrange = -1; /* default depends on input */
    param->vui.i_colorprim = 2;  /* undef */
    param->vui.i_transfer  = 2;  /* undef */
    param->vui.i_colmatrix = -1; /* default depends on input */
    param->vui.i_chroma_loc= 0;  /* left center */
    param->i_fps_num       = 25;
    param->i_fps_den       = 1;
    param->i_level_idc     = -1;
    param->i_slice_max_size = 0;
    param->i_slice_max_mbs = 0;
    param->i_slice_count = 0;

    /* Encoder parameters */
    param->i_frame_reference = 3;
    param->i_keyint_max = 250;
    param->i_keyint_min = X264_KEYINT_MIN_AUTO;
    param->i_scenecut_threshold = 40;
    param->b_constrained_intra = 0;

    param->b_deblocking_filter = 1;
    param->i_deblocking_filter_alphac0 = 0;
    param->i_deblocking_filter_beta = 0;

    param->rc.i_rc_method = X264_RC_CRF;
    param->rc.i_bitrate = 0;
    param->rc.f_rate_tolerance = 1.0;
    param->rc.i_vbv_max_bitrate = 0;
    param->rc.i_vbv_buffer_size = 0;
    param->rc.f_vbv_buffer_init = 0.9;
    param->rc.i_qp_constant = 23;
    param->rc.f_rf_constant = 23;
    param->rc.i_qp_min = 0;
    param->rc.i_qp_max = QP_MAX;
    param->rc.i_qp_step = 4;
    param->rc.f_ip_factor = 1.4;
    param->rc.f_pb_factor = 1.3;
    param->rc.i_aq_mode = X264_AQ_VARIANCE;
    param->rc.f_aq_strength = 1.0;

    param->rc.f_qcompress = 0.6;
    param->rc.b_mb_tree = 0;

    /* Log */
    param->pf_log = x264_log_default;
    param->p_log_private = NULL;
    param->i_log_level = X264_LOG_INFO;

    /* */
    param->analyse.intra = X264_ANALYSE_I4x4 ;
    param->analyse.inter = X264_ANALYSE_I4x4
                         | X264_ANALYSE_PSUB16x16 ;
    param->analyse.i_direct_mv_pred = X264_DIRECT_PRED_SPATIAL;
    param->analyse.i_me_method = X264_ME_DIA;
    param->analyse.f_psy_rd = 1.0;
    param->analyse.b_psy = 1;
    param->analyse.i_me_range = 16;
    param->analyse.i_subpel_refine = 1;
    param->analyse.b_chroma_me = 1;
    param->analyse.i_mv_range = -1; // set from level_idc
    param->analyse.i_chroma_qp_offset = 0;
    param->analyse.b_fast_pskip = 1;
    param->analyse.b_dct_decimate = 1;
    param->analyse.i_luma_deadzone[0] = 21;
    param->analyse.i_luma_deadzone[1] = 11;
    param->analyse.b_psnr = 0;
    param->analyse.b_ssim = 0;

    param->b_repeat_headers = 1;
    param->b_annexb = 1;
    param->b_aud = 0;
    param->i_nal_hrd = X264_NAL_HRD_NONE;
}

static int x264_param_apply_preset( x264_param_t *param, const char *preset )
{
    char *end;
    int i = strtol( preset, &end, 10 );
    if( *end == 0 && i >= 0 && i < sizeof(x264_preset_names)/sizeof(*x264_preset_names)-1 )
        preset = x264_preset_names[i];

    if( !strcasecmp( preset, "ultrafast" ) )
    {
        param->i_frame_reference = 1;
        param->i_scenecut_threshold = 0;
        param->b_deblocking_filter = 0;
        param->analyse.intra = 0;
        param->analyse.inter = 0;
        param->analyse.i_me_method = X264_ME_DIA;
        param->analyse.i_subpel_refine = 0;
        param->rc.i_aq_mode = 0;
    }
    else if( !strcasecmp( preset, "superfast" ) )
    {
        param->analyse.inter = X264_ANALYSE_I4x4;
        param->analyse.i_me_method = X264_ME_DIA;
        param->analyse.i_subpel_refine = 1;
        param->i_frame_reference = 1;
    }
    else
    {
        x264_log( NULL, X264_LOG_ERROR, "invalid preset '%s'\n", preset );
        return -1;
    }
    return 0;
}

static int x264_param_apply_tune( x264_param_t *param, const char *tune )
{
    char *tmp = x264_malloc( strlen( tune ) + 1 );
    if( !tmp )
        return -1;
    tmp = strcpy( tmp, tune );
    char *s = strtok( tmp, ",./-+" );
    int psy_tuning_used = 0;
    while( s )
    {
        if( !strncasecmp( s, "psnr", 4 ) )
        {
            if( psy_tuning_used++ ) goto psy_failure;
            param->rc.i_aq_mode = X264_AQ_NONE;
            param->analyse.b_psy = 0;
        }
        else if( !strncasecmp( s, "ssim", 4 ) )
        {
            if( psy_tuning_used++ ) goto psy_failure;
            param->rc.i_aq_mode = X264_AQ_AUTOVARIANCE;
            param->analyse.b_psy = 0;
        }
        else if( !strncasecmp( s, "fastdecode", 10 ) )
        {
            param->b_deblocking_filter = 0;
        }
        else if( !strncasecmp( s, "zerolatency", 11 ) )
        {
        }
        else
        {
            x264_log( NULL, X264_LOG_ERROR, "invalid tune '%s'\n", s );
            x264_free( tmp );
            return -1;
        }
        if( 0 )
        {
    psy_failure:
            x264_log( NULL, X264_LOG_WARNING, "only 1 psy tuning can be used: ignoring tune %s\n", s );
        }
        s = strtok( NULL, ",./-+" );
    }
    x264_free( tmp );
    return 0;
}

int x264_param_default_preset( x264_param_t *param, const char *preset, const char *tune )
{
    x264_param_default( param );

    if( preset && x264_param_apply_preset( param, preset ) < 0 )
        return -1;
    if( tune && x264_param_apply_tune( param, tune ) < 0 )
        return -1;
    return 0;
}

void x264_param_apply_fastfirstpass( x264_param_t *param )
{

}

static int profile_string_to_int( const char *str )
{
    if( !strcasecmp( str, "baseline" ) )
        return PROFILE_BASELINE;
    if( !strcasecmp( str, "main" ) )
        return PROFILE_MAIN;
    if( !strcasecmp( str, "high" ) )
        return PROFILE_HIGH;
    if( !strcasecmp( str, "high10" ) )
        return PROFILE_HIGH10;
    if( !strcasecmp( str, "high422" ) )
        return PROFILE_HIGH422;
    if( !strcasecmp( str, "high444" ) )
        return PROFILE_HIGH444_PREDICTIVE;
    return -1;
}

int x264_param_apply_profile( x264_param_t *param, const char *profile )
{
    if( !profile )
        return 0;

    int p = profile_string_to_int( profile );
    if( p < 0 )
    {
        x264_log( NULL, X264_LOG_ERROR, "invalid profile: %s\n", profile );
        return -1;
    }
    if( p < PROFILE_HIGH444_PREDICTIVE && ((param->rc.i_rc_method == X264_RC_CQP && param->rc.i_qp_constant <= 0) ||
        (param->rc.i_rc_method == X264_RC_CRF && (int)(param->rc.f_rf_constant) <= 0)) )
    {
        x264_log( NULL, X264_LOG_ERROR, "%s profile doesn't support lossless\n", profile );
        return -1;
    }

    return 0;
}

static int parse_enum( const char *arg, const char * const *names, int *dst )
{
    for( int i = 0; names[i]; i++ )
        if( !strcmp( arg, names[i] ) )
        {
            *dst = i;
            return 0;
        }
    return -1;
}

static int x264_atobool( const char *str, int *b_error )
{
    if( !strcmp(str, "1") ||
        !strcmp(str, "true") ||
        !strcmp(str, "yes") )
        return 1;
    if( !strcmp(str, "0") ||
        !strcmp(str, "false") ||
        !strcmp(str, "no") )
        return 0;
    *b_error = 1;
    return 0;
}

static int x264_atoi( const char *str, int *b_error )
{
    char *end;
    int v = strtol( str, &end, 0 );
    if( end == str || *end != '\0' )
        *b_error = 1;
    return v;
}

static double x264_atof( const char *str, int *b_error )
{
    char *end;
    double v = strtod( str, &end );
    if( end == str || *end != '\0' )
        *b_error = 1;
    return v;
}

#define atobool(str) ( name_was_bool = 1, x264_atobool( str, &b_error ) )
#undef atoi
#undef atof
#define atoi(str) x264_atoi( str, &b_error )
#define atof(str) x264_atof( str, &b_error )

int x264_param_parse( x264_param_t *p, const char *name, const char *value )
{
    char *name_buf = NULL;
    int b_error = 0;
    int name_was_bool;
    int value_was_null = !value;
    int i;

    if( !name )
        return X264_PARAM_BAD_NAME;
    if( !value )
        value = "true";

    if( value[0] == '=' )
        value++;

    if( strchr( name, '_' ) ) // s/_/-/g
    {
        char *c;
        name_buf = strdup(name);
        while( (c = strchr( name_buf, '_' )) )
            *c = '-';
        name = name_buf;
    }

    if( (!strncmp( name, "no-", 3 ) && (i = 3)) ||
        (!strncmp( name, "no", 2 ) && (i = 2)) )
    {
        name += i;
        value = atobool(value) ? "false" : "true";
    }
    name_was_bool = 0;

#define OPT(STR) else if( !strcmp( name, STR ) )
#define OPT2(STR0, STR1) else if( !strcmp( name, STR0 ) || !strcmp( name, STR1 ) )
    if(0);
    OPT2("level", "level-idc")
    {
        if( !strcmp(value, "1b") )
            p->i_level_idc = 9;
        else if( atof(value) < 6 )
            p->i_level_idc = (int)(10*atof(value)+.5);
        else
            p->i_level_idc = atoi(value);
    }
    OPT("sar")
    {
        b_error = ( 2 != sscanf( value, "%d:%d", &p->vui.i_sar_width, &p->vui.i_sar_height ) &&
                    2 != sscanf( value, "%d/%d", &p->vui.i_sar_width, &p->vui.i_sar_height ) );
    }
    OPT("overscan")
        b_error |= parse_enum( value, x264_overscan_names, &p->vui.i_overscan );
    OPT("videoformat")
        b_error |= parse_enum( value, x264_vidformat_names, &p->vui.i_vidformat );
    OPT("fullrange")
        b_error |= parse_enum( value, x264_fullrange_names, &p->vui.b_fullrange );
    OPT("colorprim")
        b_error |= parse_enum( value, x264_colorprim_names, &p->vui.i_colorprim );
    OPT("transfer")
        b_error |= parse_enum( value, x264_transfer_names, &p->vui.i_transfer );
    OPT("colormatrix")
        b_error |= parse_enum( value, x264_colmatrix_names, &p->vui.i_colmatrix );
    OPT("chromaloc")
    {
        p->vui.i_chroma_loc = atoi(value);
        b_error = ( p->vui.i_chroma_loc < 0 || p->vui.i_chroma_loc > 5 );
    }
    OPT("fps")
    {
        if( sscanf( value, "%u/%u", &p->i_fps_num, &p->i_fps_den ) == 2 )
            ;
        else
        {
            float fps = atof(value);
            if( fps > 0 && fps <= INT_MAX/1000 )
            {
                p->i_fps_num = (int)(fps * 1000 + .5);
                p->i_fps_den = 1000;
            }
            else
            {
                p->i_fps_num = atoi(value);
                p->i_fps_den = 1;
            }
        }
    }
    OPT2("ref", "frameref")
        p->i_frame_reference = atoi(value);
    OPT("dpb-size")
        p->i_dpb_size = atoi(value);
    OPT("keyint")
    {
        if( strstr( value, "infinite" ) )
            p->i_keyint_max = X264_KEYINT_MAX_INFINITE;
        else
            p->i_keyint_max = atoi(value);
    }
    OPT2("min-keyint", "keyint-min")
    {
        p->i_keyint_min = atoi(value);
        if( p->i_keyint_max < p->i_keyint_min )
            p->i_keyint_max = p->i_keyint_min;
    }
    OPT("scenecut")
    {
        p->i_scenecut_threshold = atobool(value);
        if( b_error || p->i_scenecut_threshold )
        {
            b_error = 0;
            p->i_scenecut_threshold = atoi(value);
        }
    }
    OPT("intra-refresh")
        p->b_intra_refresh = atobool(value);
    OPT("nf")
        p->b_deblocking_filter = !atobool(value);
    OPT2("filter", "deblock")
    {
        if( 2 == sscanf( value, "%d:%d", &p->i_deblocking_filter_alphac0, &p->i_deblocking_filter_beta ) ||
            2 == sscanf( value, "%d,%d", &p->i_deblocking_filter_alphac0, &p->i_deblocking_filter_beta ) )
        {
            p->b_deblocking_filter = 1;
        }
        else if( sscanf( value, "%d", &p->i_deblocking_filter_alphac0 ) )
        {
            p->b_deblocking_filter = 1;
            p->i_deblocking_filter_beta = p->i_deblocking_filter_alphac0;
        }
        else
            p->b_deblocking_filter = atobool(value);
    }
    OPT("slice-max-size")
        p->i_slice_max_size = atoi(value);
    OPT("slice-max-mbs")
        p->i_slice_max_mbs = atoi(value);
    OPT("slice-min-mbs")
        p->i_slice_min_mbs = atoi(value);
    OPT("slices")
        p->i_slice_count = atoi(value);
    OPT("slices-max")
        p->i_slice_count_max = atoi(value);
    OPT("constrained-intra")
        p->b_constrained_intra = atobool(value);
    OPT("log")
        p->i_log_level = atoi(value);
#if HAVE_VISUALIZE
    OPT("visualize")
        p->b_visualize = atobool(value);
#endif
    OPT("dump-yuv")
        p->psz_dump_yuv = strdup(value);
    OPT2("analyse", "partitions")
    {
        p->analyse.inter = 0;
        if( strstr( value, "none" ) )  p->analyse.inter =  0;
        if( strstr( value, "all" ) )   p->analyse.inter = ~0;

        if( strstr( value, "i4x4" ) )  p->analyse.inter |= X264_ANALYSE_I4x4;
        if( strstr( value, "p8x8" ) )  p->analyse.inter |= X264_ANALYSE_PSUB16x16;
        if( strstr( value, "p4x4" ) )  p->analyse.inter |= X264_ANALYSE_PSUB8x8;
    }
    OPT2("direct", "direct-pred")
        b_error |= parse_enum( value, x264_direct_pred_names, &p->analyse.i_direct_mv_pred );
    OPT("chroma-qp-offset")
        p->analyse.i_chroma_qp_offset = atoi(value);
    OPT("me")
        b_error |= parse_enum( value, x264_motion_est_names, &p->analyse.i_me_method );
    OPT2("merange", "me-range")
        p->analyse.i_me_range = atoi(value);
    OPT2("mvrange", "mv-range")
        p->analyse.i_mv_range = atoi(value);
    OPT2("subme", "subq")
        p->analyse.i_subpel_refine = atoi(value);
    OPT("psy-rd")
    {
        if( sscanf( value, "%f", &p->analyse.f_psy_rd ) )
        {
        }
        else
        {
            p->analyse.f_psy_rd = 0;
        }
    }
    OPT("psy")
        p->analyse.b_psy = atobool(value);
    OPT("chroma-me")
        p->analyse.b_chroma_me = atobool(value);
    OPT("fast-pskip")
        p->analyse.b_fast_pskip = atobool(value);
    OPT("dct-decimate")
        p->analyse.b_dct_decimate = atobool(value);
    OPT("deadzone-inter")
        p->analyse.i_luma_deadzone[0] = atoi(value);
    OPT("deadzone-intra")
        p->analyse.i_luma_deadzone[1] = atoi(value);
    OPT("nr")
        p->analyse.i_noise_reduction = atoi(value);
    OPT("bitrate")
    {
        p->rc.i_bitrate = atoi(value);
        p->rc.i_rc_method = X264_RC_ABR;
    }
    OPT2("qp", "qp_constant")
    {
        p->rc.i_qp_constant = atoi(value);
        p->rc.i_rc_method = X264_RC_CQP;
    }
    OPT("crf")
    {
        p->rc.f_rf_constant = atof(value);
        p->rc.i_rc_method = X264_RC_CRF;
    }
    OPT("crf-max")
        p->rc.f_rf_constant_max = atof(value);
    OPT2("qpmin", "qp-min")
        p->rc.i_qp_min = atoi(value);
    OPT2("qpmax", "qp-max")
        p->rc.i_qp_max = atoi(value);
    OPT2("qpstep", "qp-step")
        p->rc.i_qp_step = atoi(value);
    OPT("ratetol")
        p->rc.f_rate_tolerance = !strncmp("inf", value, 3) ? 1e9 : atof(value);
    OPT("vbv-maxrate")
        p->rc.i_vbv_max_bitrate = atoi(value);
    OPT("vbv-bufsize")
        p->rc.i_vbv_buffer_size = atoi(value);
    OPT("vbv-init")
        p->rc.f_vbv_buffer_init = atof(value);
    OPT2("ipratio", "ip-factor")
        p->rc.f_ip_factor = atof(value);
    OPT2("pbratio", "pb-factor")
        p->rc.f_pb_factor = atof(value);
    OPT("aq-mode")
        p->rc.i_aq_mode = atoi(value);
    OPT("aq-strength")
        p->rc.f_aq_strength = atof(value);
    OPT("qcomp")
        p->rc.f_qcompress = atof(value);
    OPT("crop-rect")
        b_error |= sscanf( value, "%u,%u,%u,%u", &p->crop_rect.i_left, &p->crop_rect.i_top,
                                                 &p->crop_rect.i_right, &p->crop_rect.i_bottom ) != 4;
    OPT("psnr")
        p->analyse.b_psnr = atobool(value);
    OPT("ssim")
        p->analyse.b_ssim = atobool(value);
    OPT("aud")
        p->b_aud = atobool(value);
    OPT("sps-id")
        p->i_sps_id = atoi(value);
    OPT("global-header")
        p->b_repeat_headers = !atobool(value);
    OPT("repeat-headers")
        p->b_repeat_headers = atobool(value);
    OPT("annexb")
        p->b_annexb = atobool(value);
    OPT("nal-hrd")
        b_error |= parse_enum( value, x264_nal_hrd_names, &p->i_nal_hrd );
    else
        return X264_PARAM_BAD_NAME;
#undef OPT
#undef OPT2
#undef atobool
#undef atoi
#undef atof

    if( name_buf )
        free( name_buf );

    b_error |= value_was_null && !name_was_bool;
    return b_error ? X264_PARAM_BAD_VALUE : 0;
}

/****************************************************************************
 * x264_log:
 ****************************************************************************/
void x264_log( x264_t *h, int i_level, const char *psz_fmt, ... )
{
    if( !h || i_level <= h->param.i_log_level )
    {
        va_list arg;
        va_start( arg, psz_fmt );
        if( !h )
            x264_log_default( NULL, i_level, psz_fmt, arg );
        else
            h->param.pf_log( h->param.p_log_private, i_level, psz_fmt, arg );
        va_end( arg );
    }
}

static void x264_log_default( void *p_unused, int i_level, const char *psz_fmt, va_list arg )
{
    char *psz_prefix;
    switch( i_level )
    {
        case X264_LOG_ERROR:
            psz_prefix = "error";
            break;
        case X264_LOG_WARNING:
            psz_prefix = "warning";
            break;
        case X264_LOG_INFO:
            psz_prefix = "info";
            break;
        case X264_LOG_DEBUG:
            psz_prefix = "debug";
            break;
        default:
            psz_prefix = "unknown";
            break;
    }
    fprintf( stderr, "x264 [%s]: ", psz_prefix );
    x264_vfprintf( stderr, psz_fmt, arg );
}

/****************************************************************************
 * x264_picture_init:
 ****************************************************************************/
void x264_picture_init( x264_picture_t *pic )
{
    memset( pic, 0, sizeof( x264_picture_t ) );
    pic->i_type = X264_TYPE_AUTO;
    pic->i_qpplus1 = X264_QP_AUTO;
}

/****************************************************************************
 * x264_picture_alloc:
 ****************************************************************************/
int x264_picture_alloc( x264_picture_t *pic, int i_csp, int i_width, int i_height )
{
    typedef struct
    {
        int planes;
        int width_fix8[3];
        int height_fix8[3];
    } x264_csp_tab_t;

    static const x264_csp_tab_t x264_csp_tab[] =
    {
        [X264_CSP_I420] = { 3, { 256*1, 256/2, 256/2 }, { 256*1, 256/2, 256/2 } },
        [X264_CSP_YV12] = { 3, { 256*1, 256/2, 256/2 }, { 256*1, 256/2, 256/2 } },
        [X264_CSP_NV12] = { 2, { 256*1, 256*1 },        { 256*1, 256/2 },       }
    };

    int csp = i_csp & X264_CSP_MASK;
    if( csp <= X264_CSP_NONE || csp >= X264_CSP_MAX )
        return -1;
    x264_picture_init( pic );
    pic->img.i_csp = i_csp;
    pic->img.i_plane = x264_csp_tab[csp].planes;
    int depth_factor = 1;
    int plane_offset[3] = {0};
    int frame_size = 0;
    for( int i = 0; i < pic->img.i_plane; i++ )
    {
        int stride = (((int64_t)i_width * x264_csp_tab[csp].width_fix8[i]) >> 8) * depth_factor;
        int plane_size = (((int64_t)i_height * x264_csp_tab[csp].height_fix8[i]) >> 8) * stride;
        pic->img.i_stride[i] = stride;
        plane_offset[i] = frame_size;
        frame_size += plane_size;
    }
    pic->img.plane[0] = x264_malloc( frame_size );
    if( !pic->img.plane[0] )
        return -1;
    for( int i = 1; i < pic->img.i_plane; i++ )
        pic->img.plane[i] = pic->img.plane[0] + plane_offset[i];
    return 0;
}

/****************************************************************************
 * x264_picture_clean:
 ****************************************************************************/
void x264_picture_clean( x264_picture_t *pic )
{
    x264_free( pic->img.plane[0] );

    /* just to be safe */
    memset( pic, 0, sizeof( x264_picture_t ) );
}

/****************************************************************************
 * x264_malloc:
 ****************************************************************************/
void *x264_malloc( int i_size )
{
    uint8_t *align_buf = NULL;
#if HAVE_MALLOC_H
        align_buf = memalign( NATIVE_ALIGN, i_size );
#else
    uint8_t *buf = malloc( i_size + (NATIVE_ALIGN-1) + sizeof(void **) );
    if( buf )
    {
        align_buf = buf + (NATIVE_ALIGN-1) + sizeof(void **);
        align_buf -= (intptr_t) align_buf & (NATIVE_ALIGN-1);
        *( (void **) ( align_buf - sizeof(void **) ) ) = buf;
    }
#endif
    if( !align_buf )
        x264_log( NULL, X264_LOG_ERROR, "malloc of size %d failed\n", i_size );
    return align_buf;
}

/****************************************************************************
 * x264_free:
 ****************************************************************************/
void x264_free( void *p )
{
    if( p )
    {
#if HAVE_MALLOC_H
        free( p );
#else
        free( *( ( ( void **) p ) - 1 ) );
#endif
    }
}

/****************************************************************************
 * x264_reduce_fraction:
 ****************************************************************************/
#define REDUCE_FRACTION( name, type )\
void name( type *n, type *d )\
{                   \
    type a = *n;    \
    type b = *d;    \
    type c;         \
    if( !a || !b )  \
        return;     \
    c = a % b;      \
    while( c )      \
    {               \
        a = b;      \
        b = c;      \
        c = a % b;  \
    }               \
    *n /= b;        \
    *d /= b;        \
}

REDUCE_FRACTION( x264_reduce_fraction  , uint32_t )
REDUCE_FRACTION( x264_reduce_fraction64, uint64_t )

/****************************************************************************
 * x264_slurp_file:
 ****************************************************************************/
char *x264_slurp_file( const char *filename )
{
    int b_error = 0;
    size_t i_size;
    char *buf;
    FILE *fh = x264_fopen( filename, "rb" );
    if( !fh )
        return NULL;
    b_error |= fseek( fh, 0, SEEK_END ) < 0;
    b_error |= ( i_size = ftell( fh ) ) <= 0;
    b_error |= fseek( fh, 0, SEEK_SET ) < 0;
    if( b_error )
        goto error;
    buf = x264_malloc( i_size+2 );
    if( !buf )
        goto error;
    b_error |= fread( buf, 1, i_size, fh ) != i_size;
    if( buf[i_size-1] != '\n' )
        buf[i_size++] = '\n';
    buf[i_size] = 0;
    fclose( fh );
    if( b_error )
    {
        x264_free( buf );
        return NULL;
    }
    return buf;
error:
    fclose( fh );
    return NULL;
}

/****************************************************************************
 * x264_param2string:
 ****************************************************************************/
char *x264_param2string( x264_param_t *p, int b_res )
{
    int len = 1000;
    char *buf, *s;
    buf = s = x264_malloc( len );
    if( !buf )
        return NULL;

    if( b_res )
    {
        s += sprintf( s, "%dx%d ", p->i_width, p->i_height );
        s += sprintf( s, "fps=%u/%u ", p->i_fps_num, p->i_fps_den );
        s += sprintf( s, "bitdepth=%d ", BIT_DEPTH );
    }

    s += sprintf( s, " ref=%d", p->i_frame_reference );
    s += sprintf( s, " deblock=%d:%d:%d", p->b_deblocking_filter,
                  p->i_deblocking_filter_alphac0, p->i_deblocking_filter_beta );
    s += sprintf( s, " analyse=%#x:%#x", p->analyse.intra, p->analyse.inter );
    s += sprintf( s, " me=%s", x264_motion_est_names[ p->analyse.i_me_method ] );
    s += sprintf( s, " subme=%d", p->analyse.i_subpel_refine );
    s += sprintf( s, " psy=%d", p->analyse.b_psy );
    if( p->analyse.b_psy )
        s += sprintf( s, " psy_rd=%.2f", p->analyse.f_psy_rd );
    s += sprintf( s, " me_range=%d", p->analyse.i_me_range );
    s += sprintf( s, " chroma_me=%d", p->analyse.b_chroma_me );
    s += sprintf( s, " deadzone=%d,%d", p->analyse.i_luma_deadzone[0], p->analyse.i_luma_deadzone[1] );
    s += sprintf( s, " fast_pskip=%d", p->analyse.b_fast_pskip );
    s += sprintf( s, " chroma_qp_offset=%d", p->analyse.i_chroma_qp_offset );
    if( p->i_slice_count )
        s += sprintf( s, " slices=%d", p->i_slice_count );
    if( p->i_slice_count_max )
        s += sprintf( s, " slices_max=%d", p->i_slice_count_max );
    if( p->i_slice_max_size )
        s += sprintf( s, " slice_max_size=%d", p->i_slice_max_size );
    if( p->i_slice_max_mbs )
        s += sprintf( s, " slice_max_mbs=%d", p->i_slice_max_mbs );
    if( p->i_slice_min_mbs )
        s += sprintf( s, " slice_min_mbs=%d", p->i_slice_min_mbs );
    s += sprintf( s, " nr=%d", p->analyse.i_noise_reduction );
    s += sprintf( s, " decimate=%d", p->analyse.b_dct_decimate );

    s += sprintf( s, " constrained_intra=%d", p->b_constrained_intra );

    if( p->i_keyint_max == X264_KEYINT_MAX_INFINITE )
        s += sprintf( s, " keyint=infinite" );
    else
        s += sprintf( s, " keyint=%d", p->i_keyint_max );
    s += sprintf( s, " keyint_min=%d scenecut=%d intra_refresh=%d",
                  p->i_keyint_min, p->i_scenecut_threshold, p->b_intra_refresh );

    s += sprintf( s, " rc_lookahead=0");

    s += sprintf( s, " rc=%s mbtree=0", p->rc.i_rc_method == X264_RC_ABR ?
                               ( p->rc.i_vbv_max_bitrate == p->rc.i_bitrate ? "cbr" : "abr" )
                               : p->rc.i_rc_method == X264_RC_CRF ? "crf" : "cqp" );
    if( p->rc.i_rc_method == X264_RC_ABR || p->rc.i_rc_method == X264_RC_CRF )
    {
        if( p->rc.i_rc_method == X264_RC_CRF )
            s += sprintf( s, " crf=%.1f", p->rc.f_rf_constant );
        else
            s += sprintf( s, " bitrate=%d ratetol=%.1f",
                          p->rc.i_bitrate, p->rc.f_rate_tolerance );
        s += sprintf( s, " qcomp=%.2f qpmin=%d qpmax=%d qpstep=%d",
                      p->rc.f_qcompress, p->rc.i_qp_min, p->rc.i_qp_max, p->rc.i_qp_step );
        if( p->rc.i_vbv_buffer_size )
        {
            s += sprintf( s, " vbv_maxrate=%d vbv_bufsize=%d",
                          p->rc.i_vbv_max_bitrate, p->rc.i_vbv_buffer_size );
            if( p->rc.i_rc_method == X264_RC_CRF )
                s += sprintf( s, " crf_max=%.1f", p->rc.f_rf_constant_max );
        }
    }
    else if( p->rc.i_rc_method == X264_RC_CQP )
        s += sprintf( s, " qp=%d", p->rc.i_qp_constant );

    if( p->rc.i_vbv_buffer_size )
        s += sprintf( s, " nal_hrd=%s", x264_nal_hrd_names[p->i_nal_hrd] );
    if( p->crop_rect.i_left | p->crop_rect.i_top | p->crop_rect.i_right | p->crop_rect.i_bottom )
        s += sprintf( s, " crop_rect=%u,%u,%u,%u", p->crop_rect.i_left, p->crop_rect.i_top,
                                                   p->crop_rect.i_right, p->crop_rect.i_bottom );

    if( !(p->rc.i_rc_method == X264_RC_CQP && p->rc.i_qp_constant == 0) )
    {
        s += sprintf( s, " ip_ratio=%.2f", p->rc.f_ip_factor );
        s += sprintf( s, " aq=%d", p->rc.i_aq_mode );
        if( p->rc.i_aq_mode )
            s += sprintf( s, ":%.2f", p->rc.f_aq_strength );
    }

    return buf;
}

