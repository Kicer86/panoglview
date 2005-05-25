bin_PROGRAMS = panoglview

pkglibdir = $(PLUGIN_DIR)

panoglview_SOURCES = panoapp.cpp panoframe.cpp panointeractivecanvas.cpp panocanvas.cpp boundaries.cpp
panoglview_CXXFLAGS = $(WX_CXXFLAGS)
panoglview_LDADD = $(WX_LIBS) $(WX_GL_LIBS)

noinst_HEADERS = panoapp.h panocanvas.h panoframe.h panoglplugin.h panopluginframe.h panointeractivecanvas.h boundaries.h

#if WIN32_BUILD


#else
#
#pkglib_LTLIBRARIES = libpanoviewer.la
#
#libpanoviewer_la_SOURCES = panoglplugin.cpp panocanvas.cpp panopluginframe.cpp
#
#libpanoviewer_la_LDFLAGS = @LDFLAGS@ $(PLUGIN_LDFLAGS)
#
#

EXTRA_DIST = panoglplugin.h panoglplugin.cpp smallicon.xpm

#endif
