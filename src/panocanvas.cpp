/***************************************************************************
                          panocanvas.cpp  -  description
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

#ifdef __GNUG__
#pragma implementation
#pragma interface
#endif

// For compilers that support precompilation, includes "wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "panocanvas.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <algorithm>

#include <wx/image.h>

    
BEGIN_EVENT_TABLE(panoCanvas, wxGLCanvas)
  EVT_PAINT(panoCanvas::OnPaint)
  EVT_SIZE (panoCanvas::OnSize)
  EVT_ERASE_BACKGROUND(panoCanvas::OnEraseBackground)
END_EVENT_TABLE()

panoCanvas::panoCanvas(wxWindow *parent, int id, const wxPoint &position, const wxSize &size) :
wxGLCanvas(parent,id,position,size),
m_position(0.0,0.0,50.0),
m_aspectratio(size.GetWidth()/(double) size.GetHeight()),
m_initialized(false),
m_hasimage(false),
m_viewableTexPatches(0),
m_currentboundaries(CPanRange  (-m_position.getFov()*m_aspectratio / 2.0, m_position.getFov()*m_aspectratio / 2.0),
                    CAngleRange( -m_position.getFov()              / 2.0, m_position.getFov()               / 2.0),
                    CAngleRange(                       0.0,                     180.0 )),
m_imageboundaries(CPanRange(-180.0,180.0),CAngleRange(-90.0,90.0),CAngleRange(0.0,180.0)),
m_divisions(128)
{
  for(int i=0;i<16;++i)
    m_projectionmatrix[i] = 0;
  for(int j=0;j<16;j+=5)
    m_projectionmatrix[j] = 1.0;
}

panoCanvas::~panoCanvas()
{
  if(m_hasimage)
    deletePanorama();
}

void panoCanvas::OnEraseBackground(wxEraseEvent& event)
{
}

void panoCanvas::OnSize(wxSizeEvent &event)
{
  wxGLCanvas::OnSize(event);

  int w,h;
  GetClientSize(&w,&h);
  SetCurrent();


  glViewport(0,0,(GLint) w, (GLint) h);
  m_aspectratio = (GLdouble) w/(GLdouble) h;
  m_winsize = wxSize(w,h);
}

void panoCanvas::OnPaint(wxPaintEvent &event)
{
  /* must always be here */
  wxPaintDC dc(this);

  SetCurrent();

  if(!m_initialized){
    initGL();
  }

  position();
  // We are looking backwards
  glRotatef(180.0,0.0,1.0,0.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if(m_hasimage){
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    showPanorama();
  }
  glFlush();
  SwapBuffers();
}

void panoCanvas::position()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(m_position.getFov(), m_aspectratio, 0.01, 10.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMultMatrixf(m_projectionmatrix);
}


void panoCanvas::deletePanorama()
{
  SetCurrent();
  glDeleteTextures(m_numOfTexPatches.x*m_numOfTexPatches.y,m_textures);
  m_hasimage = false;
  delete [] m_viewableTexPatches;
}

void panoCanvas::createPanorama(const wxImage &image)
{
  if (!m_initialized)
	  initGL();

  if(m_hasimage)
    deletePanorama();

  m_numOfTexPatches.x = (int) ceil( image.GetWidth() / (float) (m_maxsize - 2));
  m_numOfTexPatches.y = (int) ceil( image.GetHeight()/ (float) (m_maxsize - 2));

  m_pixelsPerTexture.x=image.GetWidth() /m_numOfTexPatches.x;
  m_pixelsPerTexture.y=image.GetHeight()/m_numOfTexPatches.y;

  double pixelsPerTextureX_float = (float)image.GetWidth() /(float)m_numOfTexPatches.x;
  double pixelsPerTextureY_float = (float)image.GetHeight() /(float)m_numOfTexPatches.y;

  m_textures           = new GLuint[m_numOfTexPatches.x*m_numOfTexPatches.y];
  m_viewableTexPatches = new bool  [m_numOfTexPatches.x*m_numOfTexPatches.y];

  glGenTextures(m_numOfTexPatches.x*m_numOfTexPatches.y,m_textures);

  int textureindex=0;
  bool fullpan = m_imageboundaries.getPans().getMax() - m_imageboundaries.getPans().getMin() == 360;

  unsigned char *tmp = new unsigned char [m_maxsize*m_maxsize];
  memset(tmp,m_maxsize*m_maxsize,0);

  wxProgressDialog progressDialog(wxT("Working"),wxT("Generating Panorama Image"),m_numOfTexPatches.x*m_numOfTexPatches.y);

  for(int y=0;y<m_numOfTexPatches.y;y++)
    for(int x=0;x<m_numOfTexPatches.x;x++,textureindex++){
      progressDialog.Update(textureindex);
      glBindTexture  (GL_TEXTURE_2D,  m_textures[textureindex]);
      glTexParameteri(GL_TEXTURE_2D,  GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,  GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D   (GL_TEXTURE_2D,0,GL_RGB, m_maxsize, m_maxsize, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE,tmp);

      // Main texture including right and bottom texture extension (if not at right and bottom edge of image)

      glTexSubImage2D(GL_TEXTURE_2D,0,
                  0,0,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),
                  m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1),
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  x * pixelsPerTextureX_float, y * pixelsPerTextureY_float,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),
                  m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1))).GetData());

      // Extend left and top border (rightmost and bottom texture lines because of texture 
      // wrapping) with last pixel row from previous patch if not the leftmost column or topmost
      // row patch 
      if (x) {
        // extend left texture column with last texture column from patch to the left
        // leftmost patch colmn is handled separately
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_maxsize-1, 0,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1),
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (x - 1) * pixelsPerTextureX_float + m_pixelsPerTexture.x - 1,
                    y * pixelsPerTextureY_float,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1))).GetData());
      }

      if (y) {
        // extend top texture line with last texture line from patch above
        // topmost patch row is handled separately
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    0, m_maxsize-1,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    x * pixelsPerTextureX_float,
                    (y - 1) * pixelsPerTextureY_float + m_pixelsPerTexture.y - 1,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1)).GetData());
      }

      // wrap texture on a full panorama
      if (x == m_numOfTexPatches.x - 1) {
        // extend right texture column with first pixel column in image on 360 pano
        // and right texture column with last pixel column in image partial pano
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_pixelsPerTexture.x, 0,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1),
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (fullpan ? 0 : image.GetWidth() - 1), y * pixelsPerTextureY_float,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1))).GetData());
        // extend top right corner on rightmost patch line
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_pixelsPerTexture.x, m_maxsize-1,
                    1, 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (fullpan ? 0 : image.GetWidth() - 1),
                    (y ? (y - 1) * pixelsPerTextureY_float + m_pixelsPerTexture.y - 1 : 0),
                    1, 1)).GetData());
        // extend last pixel on bottom texture line on bottom patch row
        if (y == m_numOfTexPatches.y - 1)
          glTexSubImage2D(GL_TEXTURE_2D, 0,
                      m_pixelsPerTexture.x, m_pixelsPerTexture.y,
                      1, 1,
                      GL_RGB,GL_UNSIGNED_BYTE,
                      image.GetSubImage(wxRect(
                      (fullpan ? 0 : image.GetWidth() - 1), image.GetHeight() - 1,
                      1, 1)).GetData());
      }

      if (x == 0) {
        // extend left texture column with last pixel column in image on 360 pano
        // and left texture column with first pixel column in image on partial pano
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_maxsize-1, 0,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1),
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (fullpan ? image.GetWidth() - 1 : 0), y * pixelsPerTextureY_float,
                    1, m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1))).GetData());
        // Extend top left corner on left patch column
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_maxsize-1, m_maxsize-1,
                    1, 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (fullpan ? image.GetWidth() - 1 : 0),
                    (y ? (y - 1) * pixelsPerTextureY_float + m_pixelsPerTexture.y - 1 : 0),
                    1, 1)).GetData());
      }

      // extend topmost and bottommost patch rows with topmost and bottommost texture pixel
      if (y == m_numOfTexPatches.y - 1) {
        // extend bottom texture line on bottom patch row
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    0, m_pixelsPerTexture.y,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    x * pixelsPerTextureX_float, image.GetHeight() - 1,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1)).GetData());
        // extend bottom left corner on bottom patch line
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    m_maxsize-1, m_pixelsPerTexture.y,
                    1, 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (x ? ( x - 1) * pixelsPerTextureX_float + m_pixelsPerTexture.x - 1 : (fullpan ? image.GetWidth() - 1 : 0)),
                    image.GetHeight() - 1,
                    1, 1)).GetData());
      }

      if (y == 0) {
        // extend top texture line on top patch row
        glTexSubImage2D(GL_TEXTURE_2D, 0,
                    0, m_maxsize-1,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    x * pixelsPerTextureX_float, 0,
                    m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1), 1)).GetData());
        // Extend top left corner on top patch row.
        if (x)
          glTexSubImage2D(GL_TEXTURE_2D, 0,
                      m_maxsize-1, m_maxsize-1,
                      1, 1,
                      GL_RGB,GL_UNSIGNED_BYTE,
                      image.GetSubImage(wxRect(
                      (x - 1) * pixelsPerTextureX_float + m_pixelsPerTexture.x - 1, 0,
                      1, 1)).GetData());
      }

      if (x && y) {
        // extend top left corner on all patches except top row and left column
        glTexSubImage2D(GL_TEXTURE_2D,0,
                    m_maxsize-1,m_maxsize-1,
                    1,1,
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    (x - 1) * pixelsPerTextureX_float + m_pixelsPerTexture.x - 1,
                    (y - 1) * pixelsPerTextureY_float + m_pixelsPerTexture.y - 1,
                    1,1)).GetData());
      }

    }

  delete [] tmp;

  // Number of angle steps per texture patch
  m_stepsPerTexture.x = (int) ceil( m_divisions / (float) m_numOfTexPatches.x);
  m_stepsPerTexture.y = (int) ceil( m_divisions / (float) (2 * m_numOfTexPatches.y) );

  // Angular interval for each patch
  m_phiinterval = RAD(m_imageboundaries.getPans().getMax() - m_imageboundaries.getPans().getMin()) / m_numOfTexPatches.x;
  m_thetainterval = RAD(m_imageboundaries.getTilts().getMax() - m_imageboundaries.getTilts().getMin()) / m_numOfTexPatches.y;

  // Angular increment for eatch patch step
  m_phistep   = m_phiinterval  / m_stepsPerTexture.x;
  m_thetastep = m_thetainterval/ m_stepsPerTexture.y;

  // set image boundaries
  m_thetastart = RAD(m_imageboundaries.getTilts().getMin()) + M_PI_2;
  m_phistart   = RAD(m_imageboundaries.getPans().getMin()) + M_PI;

  // Image fill level for each texture patch
  m_maxtexturex = (double) m_pixelsPerTexture.x / (double) m_maxsize;
  m_maxtexturey = (double) m_pixelsPerTexture.y / (double) m_maxsize;

  m_hasimage = true;
  Refresh();
}

CBoundaries panoCanvas::calculateViewBoundaries(const CPosition &offset)
{
  CBoundaries result;

  CPosition tmp;
  CPosition withoffset = m_position + offset;

  // Choose tilt boundaries
  if(withoffset.getTilt() + m_position.getFov() / 2.0 > 0)
    tmp = getPanTilt(m_winsize.x/2,m_winsize.y,offset);
  else
    tmp = getPanTilt(0,m_winsize.y,offset);

  result.setTiltmax(tmp.getTilt());

  if(withoffset.getTilt() - m_position.getFov() / 2.0 < 0)
    tmp = getPanTilt(m_winsize.x/2,0,offset);
  else
    tmp = getPanTilt(0,0,offset);

  result.setTiltmin(tmp.getTilt());

  if(withoffset.getTilt() + m_position.getFov()/ 2.0 > 90.0)
  {
    result.setTiltmax(179.99999);
  } else if (m_position.getTilt() - m_position.getFov ()/ 2.0 < -90.0)
  {
    result.setTiltmin(0.0);
  }

  // Choose pan boundaries
  if(withoffset.getTilt() + m_position.getFov() / 2.0 > 90.0 || withoffset.getTilt() - m_position.getFov() / 2.0 < -90.0){
    result.setPanmin(0.0);
    result.setPanmax(359.99999);
  } else {
    if (withoffset.getTilt() > 0)
      tmp = getPanTilt(0,m_winsize.y,offset);
    else
      tmp = getPanTilt(0,0,offset);

    result.setPanmin(tmp.getPan());

    double panmintmp = result.getPans().getMin();
    double pantmp = withoffset.getPan();

    if(result.getPans().getMin() < 0.0)
      result.setPanmin(result.getPans().getMin() + 360.0);

    if(pantmp < 0.0)
      pantmp += 360.0;

    while(panmintmp > pantmp)
    {
      pantmp += 360.0;
    }

    result.setPanmax(result.getPans().getMin() + 2*(pantmp - panmintmp));

    while (result.getPans().getMax() < 0.0)
      result.setPanmax(result.getPans().getMax() + 360.0);

    while (result.getPans().getMax() > 360.0)
      result.setPanmax(result.getPans().getMax() - 360.0);
  }

  return result;
}

void panoCanvas::showPanorama()
{
  glCullFace(GL_BACK);

  int textureindex=0;
  double phi;
  double theta;

  m_currentboundaries = calculateViewBoundaries();

  // Set viewable texpatches
  memset(m_viewableTexPatches,0,m_numOfTexPatches.x*m_numOfTexPatches.y*sizeof(bool));

  // Calculate indices based on current angular field of view
  double panBoundaryMax = m_currentboundaries.getPans().getMax();
  double panBoundaryMin = m_currentboundaries.getPans().getMin();
  double imgBoundaryMax = m_imageboundaries.getPans().getMax() + 180;
  double imgBoundaryMin = m_imageboundaries.getPans().getMin() + 180;

  // Correct overflow
  if (panBoundaryMax < panBoundaryMin) {
    if ((imgBoundaryMin < panBoundaryMax) && (imgBoundaryMin < panBoundaryMin))
      imgBoundaryMin += 360;
    if ((imgBoundaryMax < panBoundaryMax) && (imgBoundaryMax < panBoundaryMin))
      imgBoundaryMax += 360;
    panBoundaryMax += 360;
  }

  int panminindex = floor(RAD(panBoundaryMin - imgBoundaryMin) / m_phiinterval);
  if (panminindex < 0)
    panminindex += floor(2 * M_PI / m_phiinterval);

  int panmaxindex = floor(RAD(panBoundaryMax - imgBoundaryMin) / m_phiinterval);
  if (panmaxindex < 0)
    panmaxindex += floor(2 * M_PI / m_phiinterval);

  int thetaminindex = floor(RAD(m_currentboundaries.getTilts().getMin() - (m_imageboundaries.getTilts().getMin() + 90)) / m_thetainterval);

  int thetamaxindex = floor(RAD(m_currentboundaries.getTilts().getMax() - (m_imageboundaries.getTilts().getMin() + 90)) / m_thetainterval);

  if (thetaminindex < 0)
    thetaminindex = 0;

  if (thetamaxindex >= m_numOfTexPatches.y)
    thetamaxindex = m_numOfTexPatches.y - 1;

  for(int k=thetaminindex;  k <= thetamaxindex; k++)
    for(int l=0; l < m_numOfTexPatches.x; l++)
      if (((l >= panminindex) && (l <= panmaxindex)) ||
          ((panmaxindex < panminindex) && ((l >= panminindex) || (l <= panmaxindex))))
        m_viewableTexPatches[k * m_numOfTexPatches.x + l] = 1;

  textureindex=0;

  for(int y=0;y<m_numOfTexPatches.y;y+=1){
    for(int x=0;x<m_numOfTexPatches.x;x+=1,textureindex+=1){
      if (m_viewableTexPatches[textureindex]){
        glBindTexture  (GL_TEXTURE_2D,  m_textures[textureindex]);
        theta = y * m_thetainterval - M_PI_2 + m_thetastart;
        glBegin(GL_QUAD_STRIP);
        for(int k=0;k<m_stepsPerTexture.y;k++){

          phi = (x+1) * m_phiinterval - M_PI_2 + m_phistart;
          double nexttheta =  theta + m_thetastep;

          for(int l=0; l<=m_stepsPerTexture.x;l++){
            glColor4f(1.0 - (theta + M_PI_2) / M_PI , phi / (2 * M_PI ) , 0.0 ,1.0);
            glTexCoord2f(m_maxtexturex - (l / (double) m_stepsPerTexture.x * m_maxtexturex),
                         (k+1) / (double) m_stepsPerTexture.y * m_maxtexturey);
            glVertex3d  (cos(nexttheta) * cos(phi),-sin(nexttheta),cos(nexttheta) * sin(phi));

            glTexCoord2f(m_maxtexturex - (l / (double) m_stepsPerTexture.x * m_maxtexturex),
                         k / (double) m_stepsPerTexture.y * m_maxtexturey);

            glVertex3d  (cos(theta) * cos(phi),-sin(theta),cos(theta) * sin(phi));
            phi-=m_phistep;
          }
        theta = nexttheta;
        }
      glEnd();
      }
    }
  }
}

void panoCanvas::incrementPosition(CPosition increment)
{
  GLfloat currentmatrix[16];
  glMatrixMode(GL_MODELVIEW);
  m_position += increment;

  // We have to rotate back 180 deg as tilt is stored "forwards"
  glRotatef(180.0,0.0,1.0,0.0);

  glGetFloatv(GL_MODELVIEW_MATRIX,currentmatrix);
  glRotatef(increment.getTilt(),currentmatrix[0],currentmatrix[4],currentmatrix[8]);
  glRotatef(increment.getPan(),0.0,1.0,0.0);

  glGetFloatv(GL_MODELVIEW_MATRIX,m_projectionmatrix);
  glRotatef(180.0,0.0,1.0,0.0);
}

void panoCanvas::setPosition(const CPosition &position)
{
  m_position = position;
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glRotated(m_position.getTilt(),1.0,0.0,0.0);
  glRotated(m_position.getPan(), 0.0,1.0,0.0);
  glGetFloatv(GL_MODELVIEW_MATRIX,m_projectionmatrix);
}


void panoCanvas::initGL()
{
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &m_maxsize);
  if(m_maxsize > 256)
    m_maxsize = 256;

  glClearColor(0.0,0.0,0.0,0.0);
  glShadeModel(GL_FLAT);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  m_initialized=true;
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_CULL_FACE);
}

CPosition panoCanvas::getPanTilt(int x, int y, const CPosition &offset)
{
  CPosition result     = m_position;
  CPosition withoffset = offset + m_position;

  double orig_x =  tan(RAD(withoffset.getFov())/2.0)*m_aspectratio* 2.0 * ( x / (double) m_winsize.x - 0.5);
  double orig_y =  tan(RAD(withoffset.getFov())/2.0)              * 2.0 * ( y / (double) m_winsize.y - 0.5);
  double orig_z =  -1;

  // Add 180 because we are looking backwards
  double panrad  = RAD(withoffset.getPan() + 180.0);
  double tiltrad = RAD(withoffset.getTilt());

  double cosP = cos(panrad);
  double sinP = sin(panrad);
  double cosT = cos(tiltrad);
  double sinT = sin(tiltrad);

  double rotated_x =   cosP*orig_x - sinT*sinP*orig_y - sinP*cosT*orig_z;
  double rotated_y =   cosT*orig_y - sinT*orig_z;
  double rotated_z =   sinP*orig_x + sinT*cosP*orig_y + cosP*cosT*orig_z;

  result.setPan(DEG(atan2(rotated_x,-rotated_z)));

  result.setTilt(DEG(atan2(rotated_y,sqrt(rotated_x*rotated_x+rotated_z*rotated_z))) + 90.0);
  return result;
}

const CPosition &panoCanvas::getPosition()
{
  return m_position;
}
