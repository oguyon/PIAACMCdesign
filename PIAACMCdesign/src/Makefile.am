bin_PROGRAMS = PIAACMCdesign
PIAACMCdesign_SOURCES = CLIcore.c CLIcore.h initmodules.c initmodules.h calc_bison.y calc_flex.l calc.h
PIAACMCdesign_CPPFLAGS = -I@abs_top_srcdir@/src -DCONFIGDIR=\"@abs_top_srcdir@/config\" -DSOURCEDIR=\"@abs_top_srcdir@/src\" -DABSSRCTOPDIR=\"@abs_top_srcdir@\" -fopenmp

ACLOCAL_AMFLAGS = -I m4

## For Bison and FLEX
AM_YFLAGS = -d
AM_LFLAGS = -i

noinst_HEADERS = calc.h
BUILT_SOURCES = calc_bison.h

SUBDIRS =
PIAACMCdesign_LDADD =



SUBDIRS += PIAACMCsimul
PIAACMCdesign_LDADD += PIAACMCsimul/libpiaacmcsimul.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/PIAACMCsimul

SUBDIRS += OptSystProp
PIAACMCdesign_LDADD += OptSystProp/liboptsystprop.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/OptSystProp

SUBDIRS += OpticsMaterials
PIAACMCdesign_LDADD += OpticsMaterials/libopticsmaterials.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/OpticsMaterials

SUBDIRS += coronagraphs
PIAACMCdesign_LDADD += coronagraphs/libcoronagraphs.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/coronagraphs

SUBDIRS += ZernikePolyn
PIAACMCdesign_LDADD += ZernikePolyn/libzernikepolyn.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/ZernikePolyn

SUBDIRS += WFpropagate
PIAACMCdesign_LDADD += WFpropagate/libwfpropagate.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/WFpropagate

SUBDIRS += image_basic
PIAACMCdesign_LDADD += image_basic/libimagebasic.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/image_basic

SUBDIRS += image_filter
PIAACMCdesign_LDADD += image_filter/libimagefilter.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/image_filter

SUBDIRS += kdtree
PIAACMCdesign_LDADD += kdtree/libkdtree.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/kdtree

SUBDIRS += image_gen
PIAACMCdesign_LDADD += image_gen/libimagegen.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/image_gen

SUBDIRS += linopt_imtools
PIAACMCdesign_LDADD += linopt_imtools/liblinoptimtools.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/linopt_imtools

SUBDIRS += statistic
PIAACMCdesign_LDADD += statistic/libstatistic.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/statistic

SUBDIRS += fft
PIAACMCdesign_LDADD += fft/libfft.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/fft

SUBDIRS += info
PIAACMCdesign_LDADD += info/libinfo.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/info

SUBDIRS += COREMOD_arith
PIAACMCdesign_LDADD += COREMOD_arith/libcoremodarith.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_iofits
PIAACMCdesign_LDADD += COREMOD_iofits/libcoremodiofits.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
PIAACMCdesign_LDADD += COREMOD_memory/libcoremodmemory.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
PIAACMCdesign_LDADD += COREMOD_tools/libcoremodtools.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

SUBDIRS += 00CORE
PIAACMCdesign_LDADD += 00CORE/lib00core.la
PIAACMCdesign_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

