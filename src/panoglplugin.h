/***************************************************************************
                          CPanoGLViewer.h  -  description
                             -------------------
    begin                : Thu Feb 17 2002
    copyright            : (C) 2002 by Fabian Wenzel
    email                : wenzel@tu-harburg.de

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 * This program is part of Fabian Wenzel's PhD - project. It may not be    *
 * used otherwise.                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef CPANOGLVIEWER_H
#define CPANOGLVIEWER_H

#include "cvistamodule.h"
#include "panopluginframe.h"
#include <tuv/CImageByteRGB.h>


class CPanoGLViewer : public CVistaModule {
VISTAMODULE(CPanoGLViewer)
public:
  CPanoGLViewer();
  virtual bool init();
  virtual void done();
  virtual bool process();
private:
  OUTPUTPIN(TUVision::CImageByteRGB           , m_image                 )
  PARAMETER(int                               , m_width                 )
  PARAMETER(int                               , m_height                )
  PARAMETER(vistaTraitsTypes::filename        , m_filename              )
  PARAMETER(double                            , m_initialtilt           )
  PARAMETER(double                            , m_initialpan            )
  PARAMETER(double                            , m_initialfov            )
  PARAMETER(double                            , m_secondsperround       )
  panoPluginFrame                              *m_frame;
};

#endif
