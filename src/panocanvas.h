/***************************************************************************
                          panocanvas.h  -  description
                             -------------------
    begin                : Mon Jun 2 2003
    copyright            : (C) 2003 by Fabian Wenzel
    email                : f.wenzel@gmx.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PANOCANVAS_H
#define PANOCANVAS_H

#include <wx/glcanvas.h>
#include <wx/image.h>
#include <wx/thread.h>
#include <wx/progdlg.h>
#include <wx/timer.h>
#include <GL/gl.h>
#include <vector>

#include <cmath>

#include "boundaries.h"

#ifndef M_PI
#define M_PI 3.14159263
#endif

#ifndef M_PI_2
#define M_PI_2 M_PI/2.0
#endif

template <class T> inline T DEG(const T rad){
  return rad * 180.0 / M_PI;
}

template <class T> inline T RAD(const T deg){
  return deg * M_PI / 180.0;
}

class CPosition{
public:
  CPosition() :
  m_pan (0.0),
  m_tilt(0.0),
  m_fov (0.0)
  {
  }

  CPosition(const double &pan, const double &tilt, const double &fov) :
  m_pan (pan),
  m_tilt(tilt),
  m_fov (fov)
  {
  }

  void shiftPan()
  {
    if(m_pan < 0.0)
      m_pan += 180.0;
    else
      m_pan -=180.0;
  }

  void setPan(const double &pan)
  {
    m_pan = pan;
  }

  void setTilt(const double &tilt)
  {
    m_tilt = tilt;
  }

  void setFov(const double &fov)
  {
    m_fov = fov;
  }

  const double &getPan() const
  {
    return m_pan;
  }

  const double &getTilt() const
  {
    return m_tilt;
  }

  const double &getFov() const
  {
    return m_fov;
  }

  void incrementPan(const double &increment)
  {
    m_pan += increment;
    if(m_pan > 360.0)
      m_pan-=360.0;
    else if (m_pan < 0.0)
      m_pan += 360.0;
  }

  void incrementTilt(const double &increment)
  {
    m_tilt += increment;
  }

  void incrementFov(const double &increment)
  {
    m_fov += increment;
  }

  const CPosition &operator+=(const CPosition &other)
  {
    incrementPan (other.getPan());
    incrementTilt(other.getTilt());
    incrementFov (other.getFov());
    return *this;
  }

  CPosition operator+(const CPosition &other) const
  {
    CPosition result = *this;
    result += other;
    return result;
  }
private:
  double m_pan;
  double m_tilt;
  double m_fov;
};

/**
  *@author Fabian Wenzel
  */

class panoCanvas : public wxGLCanvas  {
public:
  panoCanvas(wxWindow *parent, int id, const wxPoint &position=wxDefaultPosition, const wxSize &size=wxDefaultSize);
  virtual ~panoCanvas();
  void OnPaint(wxPaintEvent &event);
  void OnSize(wxSizeEvent &event);
  void OnEraseBackground(wxEraseEvent& event);

  CBoundaries calculateViewBoundaries(const CPosition &offset = CPosition());

  void incrementPosition(CPosition increment);
  void setPosition(const CPosition &position);
  const CPosition &getPosition();
  void createPanorama(const wxImage &image);
  void deletePanorama();

  const CBoundaries &getImageBoundaries() const
  {
    return m_imageboundaries;
  }

  void setImageBoundaries(const CBoundaries &newboundaries)
  {
    m_imageboundaries = newboundaries;
  }

protected:
  void position();
  void showPanorama();
  void initGL();
  CPosition getPanTilt(int x, int y, const CPosition &offset = CPosition());
    
  GLfloat m_projectionmatrix[16];
  GLdouble m_aspectratio;
  GLuint  *m_textures;
  wxSize   m_winsize;
  bool     m_hasimage;
  GLint    m_maxsize;
  wxPoint  m_numOfTexPatches;
  wxPoint  m_stepsPerTexture;
  wxPoint  m_pixelsPerTexture;
  CPosition m_position;
  bool     *m_viewableTexPatches;
  int      m_divisions;
  double   m_phiinterval;
  double   m_thetainterval;
  double   m_phistep;
  double   m_thetastep;
  double   m_phistart;
  double   m_thetastart;
  double   m_maxtexturex;
  double   m_maxtexturey;
  CBoundaries m_currentboundaries;
  CBoundaries m_imageboundaries;
  bool     m_initialized;
  DECLARE_EVENT_TABLE();
};

#endif
