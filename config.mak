SRCPATH=.
prefix=/usr/local
exec_prefix=${prefix}
bindir=${exec_prefix}/bin
libdir=${exec_prefix}/lib
includedir=${prefix}/include
ARCH=X86
SYS=LINUX
CC=gcc
CFLAGS=-Wshadow -O0 -g -m32  -Wall -I. -I$(SRCPATH) -pg -std=gnu99 -fno-tree-vectorize
DEPMM=-MM -g0
DEPMT=-MT
LD=gcc -o 
LDFLAGS=-m32  -pg -lm
LIBX264=libx264.a
AR=ar rc 
RANLIB=ranlib
STRIP=strip
AS=
ASFLAGS= -O0 -f elf -DHAVE_ALIGNED_STACK=1 -DHIGH_BIT_DEPTH=0 -DBIT_DEPTH=8
RC=
RCFLAGS=
EXE=
HAVE_GETOPT_LONG=1
DEVNULL=/dev/null
PROF_GEN_CC=-fprofile-generate
PROF_GEN_LD=-fprofile-generate
PROF_USE_CC=-fprofile-use
PROF_USE_LD=-fprofile-use
HAVE_OPENCL=no
default: cli
install: install-cli
LDFLAGSCLI = 
CLI_LIBX264 = $(LIBX264)
