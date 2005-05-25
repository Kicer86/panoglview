bin_PROGRAMS = panoglview

pkglibdir = $(PLUGIN_DIR)

AM_CXXFLAGS = $(PANOGLCXXFLAGS) $(PLUGIN_CXXFLAGS) 

panoglview_SOURCES = panoapp.cpp panoframe.cpp panointeractivecanvas.cpp panocanvas.cpp boundaries.cpp

panoglview_LDADD = $(PANOGLLDADD)

noinst_HEADERS = panoapp.h panoapp.h panoapp.h panocanvas.h panoframe.h panoglplugin.h panocanvas.h panoframe.h panoglplugin.h panopluginframe.h panocanvas.h panoframe.h panoglplugin.h panopluginframe.h panointeractivecanvas.h boundaries.h

#if WIN32_BUILD

EXTRA_DIST = panoglplugin.h panoglplugin.cpp

#else
#
#pkglib_LTLIBRARIES = libpanoviewer.la
#
#libpanoviewer_la_SOURCES = panoglplugin.cpp panocanvas.cpp panopluginframe.cpp
#
#libpanoviewer_la_LDFLAGS = @LDFLAGS@ $(PLUGIN_LDFLAGS)
#
#

EXTRA_DIST = smallicon.xpm

#endif
