//
//
// C++ Implementation: $MODULE$
//
// Description: 
//
//
// Author: Fabian Wenzel <f.wenzel@gmx.net>, (C) 2003
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "panointeractivecanvas.h"
#include "panoframe.h"

#include <wx/log.h>
#include <cmath>
#include <algorithm>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

panoDropTarget::panoDropTarget(panoFrame *frame) :
p_frame(frame)
{
}

bool panoDropTarget::OnDropFiles(wxCoord x, wxCoord y, const wxArrayString& filenames)
{
  p_frame->openArgumentFile(filenames[0]); 
  return true;
}

  
static const int ZOOMIN  = '+';
static const int ZOOMOUT = '-';

enum {
  ID_TIMER=1
};

BEGIN_EVENT_TABLE(panoInteractiveCanvas, panoCanvas)
  EVT_TIMER       (ID_TIMER,panoInteractiveCanvas::OnTimer)
  EVT_MOUSE_EVENTS(panoInteractiveCanvas::OnMouse         )
  EVT_KEY_DOWN    (panoInteractiveCanvas::OnKeyDown       )
  EVT_KEY_UP      (panoInteractiveCanvas::OnKeyUp         )
  EVT_PAINT       (panoInteractiveCanvas::OnPaint         )
  EVT_SIZE        (panoInteractiveCanvas::OnSize          )
END_EVENT_TABLE()

panoInteractiveCanvas::panoInteractiveCanvas(wxWindow* parent, int id, const wxPoint& position, const wxSize& size):
panoCanvas(parent, id, position, size),
m_timerelapse(20),
m_timer(this,ID_TIMER),
m_givenboundaries(),
p_frame((panoFrame *)parent),
m_boundarymode(-1),
m_showboundaries(true),
m_useboundaries(false),
m_zoomindown(false),
m_zoomoutdown(false),
m_leftdown(false),
m_rightdown(false),
m_updown(false),
m_downdown(false),
m_leftbuttondown(false),
m_wheelrot(0)
{
  panoDropTarget *target = new panoDropTarget(p_frame);
  SetDropTarget(target);
  SetCursor(*wxCROSS_CURSOR);
}


panoInteractiveCanvas::~panoInteractiveCanvas()
{
}

void panoInteractiveCanvas::OnSize(wxSizeEvent &event)
{
  wxGLCanvas::OnSize(event);

  int w,h;
  GetClientSize(&w,&h);
  SetCurrent();

  glViewport(0,0,(GLint) w, (GLint) h);
  
  if((double) w / (double) h > m_aspectratio)
    m_position.setFov(m_position.getFov() /w * h * m_aspectratio);
  else
    m_position.setFov(m_position.getFov() * w / h / m_aspectratio);

  if (m_givenboundaries.getFovs().validMin() && m_position.getFov() < m_givenboundaries.getFovs().getMin())
    m_position.setFov( m_givenboundaries.getFovs().getMin() );
  else if (m_position.getFov() < 1.0)
    m_position.setFov(1.0);

  m_aspectratio = (GLdouble) w/(GLdouble) h;
  m_winsize = wxSize(w,h);
}

void panoInteractiveCanvas::OnPaint(wxPaintEvent &event)
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
    glEnable(GL_TEXTURE_2D);
    showPanorama();
  }

  glDisable(GL_TEXTURE_2D);

  if(m_boundarymode != -1)
    showActiveBoundary();

  if(m_showboundaries)
    showAllBoundaries();

  glFlush();
  SwapBuffers();
}

void panoInteractiveCanvas::OnKeyDown(wxKeyEvent& event)
{
  switch(event.GetKeyCode()){

  case ZOOMIN:
  case WXK_NUMPAD_ADD:
    m_zoomindown = true;
    break;

  case ZOOMOUT:
  case WXK_NUMPAD_SUBTRACT:
    m_zoomoutdown = true;
    break;

  case WXK_ESCAPE:
    p_frame->quitFullscreen();
    break;

  case WXK_RIGHT:
  case WXK_NUMPAD_RIGHT:
    m_rightdown = true;
    break;

  case WXK_LEFT:
  case WXK_NUMPAD_LEFT:
    m_leftdown = true;
    break;

  case WXK_UP:
  case WXK_NUMPAD_UP:
    m_updown = true;
    break;

  case WXK_DOWN:
  case WXK_NUMPAD_DOWN:
    m_downdown = true;
    break;

  default:
    break;
  }
  m_timer.Start(m_timerelapse);
  event.Skip();
}

void panoInteractiveCanvas::OnKeyUp(wxKeyEvent &event)
{
  switch(event.GetKeyCode()){

  case ZOOMIN:
  case WXK_NUMPAD_ADD:    
    m_zoomindown = false;
    break;

  case ZOOMOUT:
  case WXK_NUMPAD_SUBTRACT:
    m_zoomoutdown = false;
    break;

  case WXK_LEFT:
  case WXK_NUMPAD_LEFT:
    m_leftdown = false;
    break;

  case WXK_RIGHT:
  case WXK_NUMPAD_RIGHT:
    m_rightdown = false;
    break;

  case WXK_UP:
  case WXK_NUMPAD_UP:
    m_updown = false;
    break;

  case WXK_DOWN:
  case WXK_NUMPAD_DOWN:
    m_downdown = false;
    break;

  default:
    break;
  }

  if(! (m_zoomindown || m_zoomoutdown || m_leftdown || m_rightdown || m_updown || m_downdown || m_leftbuttondown))
    m_timer.Stop();

  updateStatusText();
    event.Skip();
}

void panoInteractiveCanvas::OnMouse(wxMouseEvent& event)
{
  bool neednewtrigger = false;

  if (event.Leaving())
    m_zoomindown = m_zoomoutdown = m_leftdown = m_rightdown = m_updown = m_downdown = false;

  if (event.LeftDown())
  {
    m_leftbuttondown = true;
    neednewtrigger = true;
    m_clickposition = event.GetPosition();
    m_timer.Start(m_timerelapse);
  } else if (event.LeftUp()){
    m_leftbuttondown = false;
    if(! (m_zoomindown || m_zoomoutdown || m_leftdown || m_rightdown || m_updown || m_downdown ))
      m_timer.Stop();
    if(m_boundarymode != -1)
    {
      setBoundary(event.GetPosition());
    }
  } else {
    if(m_boundarymode != -1)
      Refresh();
  }

  if(m_leftbuttondown)
  {
    m_diff = event.GetPosition() - m_clickposition;
  }

  m_wheelrot += event.GetWheelRotation();

  int wheeldelta = event.GetWheelDelta();
  if(wheeldelta == 0)
    wheeldelta = 120;

  if(abs(m_wheelrot) >= wheeldelta){
    CPosition increment;
    if(m_wheelrot < 0){
      incrementPosition(CPosition(0.0,0.0, 0.05*m_position.getFov()));
      m_wheelrot += wheeldelta;
    } else {
      incrementPosition(CPosition(0.0,0.0,-0.05*m_position.getFov()));
      m_wheelrot -= wheeldelta;
    }

    Refresh();
  }

  m_currentpos = event.GetPosition();
  updateStatusText();
}

void panoInteractiveCanvas::OnTimer(wxTimerEvent &event)
{
  CPosition increment;

  if(m_leftbuttondown){
    if(m_diff.y){
      double angle = log(1 + 2 * abs(m_diff.y)/(float) m_winsize.x)*m_position.getFov()/4.0*( 1 + ( abs(m_diff.y) < abs(m_diff.x) ));
      increment.setTilt(m_diff.y < 0 ? -angle:angle);
    }
    if(m_diff.x){
      double angle = log(1 + 2*abs(m_diff.x)/(float) m_winsize.y)*m_position.getFov()/4.0*(1 + ( abs(m_diff.x) < abs(m_diff.y) ));
      increment.setPan(m_diff.x < 0 ? -angle:angle);
    }
  }


  if(m_zoomindown)
    increment.incrementFov(0.03*m_position.getFov());

  if(m_zoomoutdown)
    increment.incrementFov(0.03*m_position.getFov());

  if(m_leftdown)
    increment.incrementPan(-0.02*m_position.getFov()*m_aspectratio);

  if(m_rightdown)
    increment.incrementPan(0.02*m_position.getFov()*m_aspectratio);

  if(m_updown)
    increment.incrementTilt(-0.02*m_position.getFov());

  if(m_downdown)
    increment.incrementTilt(0.02*m_position.getFov());

  incrementPosition(increment);
  updateStatusText();
  Refresh();
}

void panoInteractiveCanvas::incrementPosition(CPosition increment)
{
  if(!m_useboundaries){
    CPosition withincrement = m_position + increment;
    if ( (withincrement.getTilt()) <= -90.0 || withincrement.getTilt() >= 90.0 )
      increment.setTilt(0.0);

    if(std::min(withincrement.getFov(),(withincrement.getFov())*m_aspectratio) <= 1.0  ||
       std::max(withincrement.getFov(),(withincrement.getFov())*m_aspectratio) >= 180.0 )
      increment.setFov(0.0);
  } else {

     for(int x = 0; x < m_winsize.x; x+= (m_winsize.x -1) / 2)
      for(int y = 0; y < m_winsize.y; y+= m_winsize.y -1){
        CPosition newpos = getPanTilt(x,y,CPosition(increment.getPan(),increment.getTilt(),0.0));

        newpos.shiftPan();
          
        if(!m_givenboundaries.inPanRange (newpos.getPan()) || sign(increment.getPan()) == m_stickypandirection){
          m_stickypandirection = sign(increment.getPan());
          increment.setPan(0.0);
        } else
          m_stickypandirection = 0;

        newpos = getPanTilt(x,y,CPosition(increment.getPan(),increment.getTilt(),0.0));
        newpos.shiftPan();
        
        if(!m_givenboundaries.inTiltRange (newpos.getTilt()-90.0) || !m_givenboundaries.inPanRange(newpos.getPan()) || sign(increment.getTilt()) == m_stickytiltdirection ){
          m_stickytiltdirection = sign(increment.getTilt());
          increment.setTilt(0.0);
        } else
          m_stickytiltdirection = 0;

        newpos = getPanTilt(x,y,increment);
        newpos.shiftPan();

        if(!m_givenboundaries.inPanRange(newpos.getPan()) || !m_givenboundaries.inTiltRange(newpos.getTilt()-90.0))
          increment.setFov(0.0);
      }
  }
  panoCanvas::incrementPosition(increment);
}

void panoInteractiveCanvas::updateStatusText()
{
  CPosition position;

  position = getPanTilt(m_currentpos.x,m_currentpos.y);
  position.shiftPan();

  p_frame->SetStatusText(wxString::Format(_("Mouse (P/T) = (%3.2f/%3.2f)"),position.getPan(),position.getTilt()-90.0),2);
  position = getPanTilt(m_winsize.x/2,m_winsize.y/2);
  position.shiftPan();

  p_frame->SetStatusText(wxString::Format(_("View (P/T/FOV) = (%3.2f,%3.2f,%3.2f)"),position.getPan(),position.getTilt()-90.0,m_position.getFov()),1);
}

void panoInteractiveCanvas::setBoundaryMode(int boundarymode)
{
  m_boundarymode = boundarymode;
  Refresh();
}

void panoInteractiveCanvas::showActiveBoundary()
{
  CPosition boundary = getPanTilt(m_currentpos.x,m_currentpos.y);
  showPanTiltLine(boundary.getPan(), 1.0,1.0,1.0,true);
  showPanTiltLine(boundary.getTilt(),1.0,1.0,1.0,false);
}

void panoInteractiveCanvas::showAllBoundaries()
{
  if(m_givenboundaries.getPans().validMin())
    showPanTiltLine(m_givenboundaries.getPans().getMin() - 180.0, 1.0,0.0,0.0,true);

  if(m_givenboundaries.getPans().validMin())
    showPanTiltLine(m_givenboundaries.getPans().getMax() - 180.0, 0.0,1.0,0.0,true);

  if(m_givenboundaries.getTilts().validMin())
    showPanTiltLine(m_givenboundaries.getTilts().getMin() + 90.0,1.0,0.0,0.0,false);

  if(m_givenboundaries.getTilts().validMax())
    showPanTiltLine(m_givenboundaries.getTilts().getMax() + 90.0,0.0,1.0,0.0,false);
}

void panoInteractiveCanvas::setBoundary(const wxPoint &position)
{
  CPosition boundary = getPanTilt(position.x,position.y);

  // respect backwards view and tilt offset
  boundary.shiftPan();
  boundary.incrementTilt(-90.0);

  switch(m_boundarymode){

  case 0:
    m_givenboundaries.setPanmin( boundary.getPan() );
    break;

  case 1:
    m_givenboundaries.setPanmax( boundary.getPan() );
    break;

  case 2:
    if(m_givenboundaries.getTilts().validMax() ||
       boundary.getTilt() < m_givenboundaries.getTilts().getMax() )
      m_givenboundaries.setTiltmin(boundary.getTilt());
    else
      wxMessageBox(_("Cannot set minimal tilt greater than maximal"),_("Error"),wxOK);
    break;

  case 3:
    if(m_givenboundaries.getTilts().getMin() ||
       boundary.getTilt() > m_givenboundaries.getTilts().getMin() )
      m_givenboundaries.setTiltmax(boundary.getTilt());
    else
      wxMessageBox(_("Cannot set maximal tilt less than minimal"),_("Error"),wxOK);
    break;

  default:
    break;
  }
  setBoundaryMode(-1);
}

void panoInteractiveCanvas::showPanTiltLine(double degree, float red, float green, float blue, bool pan)
{
  degree = RAD(degree);

  glColor4f(red,green,blue,0.2);

  if(pan){
    glBegin(GL_LINE_STRIP);
    for(float tilt = 0.0; tilt <= RAD(180.0); tilt += RAD(360.0/128))
      glVertex3d(sin(degree)*sin(tilt),cos(tilt),-cos(degree)*sin(tilt));
  } else {
    glBegin(GL_LINE_LOOP);
    for(float pan  = 0.0; pan  <= RAD(360.0); pan  += RAD(360.0/128))
      glVertex3d(sin(pan)*sin(degree),cos(degree),-cos(pan)*sin(degree));
  }
  glEnd();
}

void panoInteractiveCanvas::enableShowBoundaries(bool show)
{
  m_showboundaries = show;
  Refresh();
}

void panoInteractiveCanvas::enableUseBoundaries(bool use)
{
  m_useboundaries = use;
  if(use){                                                     
    CPosition position = m_position;
    m_stickypandirection = m_stickytiltdirection = 0;
    if (m_givenboundaries.getPans().validRange()){
      if( m_givenboundaries.getPans().getMin() > m_givenboundaries.getPans().getMax() ) {
        position.setPan ((m_givenboundaries.getPans().getMax()  + m_givenboundaries.getPans().getMin() - 360.0 )  / 2.0);
        position.setFov (std::min(position.getFov(), (m_givenboundaries.getPans().getMax() + 360.0 - m_givenboundaries.getPans().getMin())/m_aspectratio));
      } else {
        position.setPan ((m_givenboundaries.getPans().getMax() + m_givenboundaries.getPans().getMin() )  / 2.0);
        position.setFov (std::min(position.getFov(), (m_givenboundaries.getPans().getMax() - m_givenboundaries.getPans().getMin())/m_aspectratio));
      }
    } else
      position.setPan (0.0);
    
    if (m_givenboundaries.getTilts().validRange())
    {
      position.setTilt ((m_givenboundaries.getTilts().getMin() + m_givenboundaries.getTilts().getMax() )/ 2.0);
      position.setFov  (std::min(position.getFov(), m_givenboundaries.getTilts().getMax() - m_givenboundaries.getTilts().getMin()));
    } else
      position.setTilt (0.0);

    position.setFov(position.getFov() * 0.8);
    
    setPosition( position );
  }
  Refresh();
                
}

void panoInteractiveCanvas::resetGivenBoundaries()
{
  setBoundaryMode(-1);
  m_givenboundaries = CBoundaries();
}

int panoInteractiveCanvas::sign(const double &value) const
{
  if(value > 0.0)
    return 1;
    
  if(value < 0.0)
    return -1;

  return 0;
}