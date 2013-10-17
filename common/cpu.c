/*****************************************************************************
 * cpu.c: cpu detection
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

#define _GNU_SOURCE // for sched_getaffinity
#include "common.h"
#include "cpu.h"

#if HAVE_POSIXTHREAD && SYS_LINUX
#include <sched.h>
#endif
#if SYS_BEOS
#include <kernel/OS.h>
#endif
#if SYS_MACOSX || SYS_FREEBSD
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#if SYS_OPENBSD
#include <sys/param.h>
#include <sys/sysctl.h>
#include <machine/cpu.h>
#endif

const x264_cpu_name_t x264_cpu_names[] =
{
    {"", 0},
};

#if (ARCH_PPC && SYS_LINUX) || (ARCH_ARM && !HAVE_NEON)
#include <signal.h>
#include <setjmp.h>
static sigjmp_buf jmpbuf;
static volatile sig_atomic_t canjump = 0;

static void sigill_handler( int sig )
{
    if( !canjump )
    {
        signal( sig, SIG_DFL );
        raise( sig );
    }

    canjump = 0;
    siglongjmp( jmpbuf, 1 );
}
#endif
uint32_t x264_cpu_detect( void )
{
    return 0;
}

int x264_cpu_num_processors( void )
{
#if !HAVE_THREAD
    return 1;

#elif SYS_WINDOWS
    return x264_pthread_num_processors_np();

#elif SYS_CYGWIN || SYS_SunOS
    return sysconf( _SC_NPROCESSORS_ONLN );

#elif SYS_LINUX
    cpu_set_t p_aff;
    memset( &p_aff, 0, sizeof(p_aff) );
    if( sched_getaffinity( 0, sizeof(p_aff), &p_aff ) )
        return 1;
#if HAVE_CPU_COUNT
    return CPU_COUNT(&p_aff);
#else
    int np = 0;
    for( unsigned int bit = 0; bit < 8 * sizeof(p_aff); bit++ )
        np += (((uint8_t *)&p_aff)[bit / 8] >> (bit % 8)) & 1;
    return np;
#endif

#elif SYS_BEOS
    system_info info;
    get_system_info( &info );
    return info.cpu_count;

#elif SYS_MACOSX || SYS_FREEBSD || SYS_OPENBSD
    int ncpu;
    size_t length = sizeof( ncpu );
#if SYS_OPENBSD
    int mib[2] = { CTL_HW, HW_NCPU };
    if( sysctl(mib, 2, &ncpu, &length, NULL, 0) )
#else
    if( sysctlbyname("hw.ncpu", &ncpu, &length, NULL, 0) )
#endif
    {
        ncpu = 1;
    }
    return ncpu;

#else
    return 1;
#endif
}
