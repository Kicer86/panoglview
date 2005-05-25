bin_PROGRAMS = panoglview

pkglibdir = $(PLUGIN_DIR)

panoglview_SOURCES = panoapp.cpp panoframe.cpp panointeractivecanvas.cpp panocanvas.cpp boundaries.cpp
panoglview_CXXFLAGS = $(WX_CXXFLAGS)
panoglview_LDADD = $(WX_LIBS) $(WX_GL_LIBS)

noinst_HEADERS = panoapp.h panocanvas.h panoframe.h panoglplugin.h panopluginframe.h panointeractivecanvas.h boundaries.h

# Currently the plugin sources do not exist
#PLUGIN_SRC = panoglplugin.cpp panopluginframe.cpp
#
#if USE_PLUGIN
#PLUGIN_DIST = 
#pkglib_LTLIBRARIES = libpanoviewer.la
#libpanoviewer_la_CXXFLAGS = $(PLUGIN_CXXFLAGS) $(WX_CXXFLAGS)
#libpanoviewer_la_SOURCES = $(PLUGIN_SRC) panocanvas.cpp
#libpanoviewer_la_LDFLAGS = $(WX_LIBS) $(WX_GL_LIBS)
#else
#PLUGIN_DIST = $(PLUGIN_SRC) 
#endif

#EXTRA_DIST = $(PLUGIN_DIST) smallicon.xpm
EXTRA_DIST = smallicon.xpm

