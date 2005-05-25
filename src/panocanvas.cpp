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
m_divisions(128)
{
  for(int i=0;i<16;++i)
    m_projectionmatrix[i] = 0;
  for(int i=0;i<16;i+=5)
    m_projectionmatrix[i] = 1.0;
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
  if(m_hasimage)
    deletePanorama();

  m_numOfTexPatches.x = (int) ceil( image.GetWidth() / (float) m_maxsize);
  m_numOfTexPatches.y = (int) ceil( image.GetHeight()/ (float) m_maxsize);

  if(image.GetWidth() % m_maxsize == 0)
    m_numOfTexPatches.x++;

  if(image.GetHeight() % m_maxsize == 0)
    m_numOfTexPatches.y++;

  m_pixelsPerTexture.x=image.GetWidth() /m_numOfTexPatches.x;
  m_pixelsPerTexture.y=image.GetHeight()/m_numOfTexPatches.y;


  m_textures           = new GLuint[m_numOfTexPatches.x*m_numOfTexPatches.y];
  m_viewableTexPatches = new bool  [m_numOfTexPatches.x*m_numOfTexPatches.y];

  glGenTextures(m_numOfTexPatches.x*m_numOfTexPatches.y,m_textures);

  int textureindex=0;

  unsigned char *tmp = new unsigned char [m_maxsize*m_maxsize];
  memset(tmp,m_maxsize*m_maxsize,0);

  wxProgressDialog progressDialog(wxT("Working"),wxT("Generating Panoramaimage"),m_numOfTexPatches.x*m_numOfTexPatches.y);

  for(int y=0;y<m_numOfTexPatches.y;y++)
    for(int x=0;x<m_numOfTexPatches.x;x++,textureindex++){
      progressDialog.Update(textureindex);
      glBindTexture  (GL_TEXTURE_2D,  m_textures[textureindex]);
      glTexParameteri(GL_TEXTURE_2D,  GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,  GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D   (GL_TEXTURE_2D,0,GL_RGB, m_maxsize, m_maxsize, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE,tmp);

  // Main texture including right and bottom extension (if not at right and bottom edge of image)

      glTexSubImage2D(GL_TEXTURE_2D,0,
                  0,0,
                  m_pixelsPerTexture.x+(x!=m_numOfTexPatches.x - 1),
                  m_pixelsPerTexture.y+(y!=m_numOfTexPatches.y - 1),
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  x*m_pixelsPerTexture.x,y*m_pixelsPerTexture.y,
                  m_pixelsPerTexture.x + ( x != m_numOfTexPatches.x - 1) ,
                  m_pixelsPerTexture.y + ( y != m_numOfTexPatches.y - 1))).GetData());

  // Left extension: left image column at rightmost texture column
      if(x){
        glTexSubImage2D(GL_TEXTURE_2D,0,
                    m_maxsize-1,0,
                    1,m_pixelsPerTexture.y+ (y != m_numOfTexPatches.y - 1),
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    x*m_pixelsPerTexture.x-1,y*m_pixelsPerTexture.y,
                    1,m_pixelsPerTexture.y+ (y != m_numOfTexPatches.y - 1))).GetData());
      } else {
   // Left extension: If x = 0: left image column = rightmost image column at rightmost texture column

        glTexSubImage2D(GL_TEXTURE_2D,0,
                    m_maxsize-1,0,
                    1,m_pixelsPerTexture.y + (y != m_numOfTexPatches.y - 1),
                    GL_RGB,GL_UNSIGNED_BYTE,
                    image.GetSubImage(wxRect(
                    image.GetWidth()-1,y*m_pixelsPerTexture.y,
                    1,m_pixelsPerTexture.y+ (y != m_numOfTexPatches.y - 1))).GetData());
     }

   // Top extension: upper image row at bottommost texture row
    if(y){
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  0,m_maxsize-1,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),1,
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  x*m_pixelsPerTexture.x,y*m_pixelsPerTexture.y-1,
                  m_pixelsPerTexture.x +(x != m_numOfTexPatches.x - 1),1)).GetData());
    } else {
    // Top extension: If y = 0: upper image row (= top image row) at bottommost texture row
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  0,m_maxsize-1,
                  m_pixelsPerTexture.x+ (x != m_numOfTexPatches.x - 1),1,
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  x*m_pixelsPerTexture.x,0,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),1)).GetData());
    }

  // Topleft extension: upper left pixel at bottommost right texture
    if(x && y){
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  m_maxsize-1,m_maxsize-1,
                  1,1,
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  x*m_pixelsPerTexture.x-1,y*m_pixelsPerTexture.y-1,
                  1,1)).GetData());

    } else {
      //Topleft extension: Correct upper left pixel for non-top y, top one for y=0 else
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  m_maxsize-1,m_maxsize-1,
                  1,1,
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  (x ? x*m_pixelsPerTexture.x-1 : image.GetWidth()-1),(y ? y*m_pixelsPerTexture.y - 1 : 0),
                  1,1)).GetData());

    }

  // If at right of bottom edge of image, extend with wraparound texture
    if( x == m_numOfTexPatches.x - 1){
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  m_pixelsPerTexture.x,0,
                  1,m_pixelsPerTexture.y + ( y != m_numOfTexPatches.y - 1),
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  0,y*m_pixelsPerTexture.y,
                  1,m_pixelsPerTexture.y + ( y != m_numOfTexPatches.y - 1))).GetData());
    }

    if( y == m_numOfTexPatches.y - 1){
      glTexSubImage2D(GL_TEXTURE_2D,0,
                  0,m_pixelsPerTexture.y,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),1,
                  GL_RGB,GL_UNSIGNED_BYTE,
                  image.GetSubImage(wxRect(
                  0,image.GetHeight()-1,
                  m_pixelsPerTexture.x + (x != m_numOfTexPatches.x - 1),1)).GetData());
    }

    // Bottom right corner is missing
    if ( y == m_numOfTexPatches.y -1 && x == m_numOfTexPatches.x - 1){
      glTexSubImage2D(GL_TEXTURE_2D,0,
		      m_pixelsPerTexture.x,m_pixelsPerTexture.y,
		      1,1,
		      GL_RGB,GL_UNSIGNED_BYTE,
		      image.GetSubImage(wxRect(0,image.GetHeight()-1,1,1)).GetData());
    }

    // Top right corner is missing:
    if ( !y && x == m_numOfTexPatches.x -1)
      glTexSubImage2D(GL_TEXTURE_2D,0,
		      m_pixelsPerTexture.x,m_maxsize-1,
		      1,1,
		      GL_RGB,GL_UNSIGNED_BYTE,
		      image.GetSubImage(wxRect(0,0,1,1)).GetData());

    // Bottom Left Corner is missing:
    if ( !x && y == m_numOfTexPatches.y -1 )
      glTexSubImage2D(GL_TEXTURE_2D,0,
		      m_maxsize-1,m_pixelsPerTexture.y,
		      1,1,
		      GL_RGB,GL_UNSIGNED_BYTE,
		      image.GetSubImage(wxRect(image.GetWidth()-1,image.GetHeight()-1,1,1)).GetData());

  }
  delete [] tmp;

  // Number of angle steps per texture patch
  m_stepsPerTexture.x = (int) ceil( m_divisions / (float) m_numOfTexPatches.x);
  m_stepsPerTexture.y = (int) ceil( m_divisions / (float) (2 * m_numOfTexPatches.y) );

  // Angular interval for each patch
  m_phiinterval   = 2 * M_PI/m_numOfTexPatches.x;
  m_thetainterval =     M_PI/m_numOfTexPatches.y;

  // Angular increment for eatch patch step
  m_phistep   = m_phiinterval  / m_stepsPerTexture.x;
  m_thetastep = m_thetainterval/ m_stepsPerTexture.y;

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
  int panminindex = (int) (RAD(m_currentboundaries.getPans().getMin())/m_phiinterval + m_numOfTexPatches.x) % m_numOfTexPatches.x;
  int panmaxindex = (int) (RAD(m_currentboundaries.getPans().getMax())/m_phiinterval + m_numOfTexPatches.x) % m_numOfTexPatches.x;

  int thetaminindex = (int) (RAD(m_currentboundaries.getTilts().getMin()) / m_thetainterval);
  int thetamaxindex = (int) (RAD(m_currentboundaries.getTilts().getMax()) / m_thetainterval);

  if(thetamaxindex < thetaminindex)
    thetamaxindex = thetamaxindex + m_numOfTexPatches.y;

  if(panmaxindex < panminindex)
    panmaxindex = panmaxindex + m_numOfTexPatches.x;

  for(int k=thetaminindex;  k <= thetamaxindex; k++)
    for(int l=panminindex; l <= panmaxindex; l++)
      m_viewableTexPatches[( k % m_numOfTexPatches.y) * m_numOfTexPatches.x + (l % m_numOfTexPatches.x)] = 1;

  textureindex=0;

  for(int y=0;y<m_numOfTexPatches.y;y++){
    for(int x=0;x<m_numOfTexPatches.x;x++,textureindex++){
      if ( m_viewableTexPatches[textureindex]){
        glBindTexture  (GL_TEXTURE_2D,  m_textures[textureindex]);
        theta = y * m_thetainterval - M_PI_2;
        glBegin(GL_QUAD_STRIP);
        for(int k=0;k<m_stepsPerTexture.y;k++){

          phi = (x+1) * m_phiinterval - M_PI_2;
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
