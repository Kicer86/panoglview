//
//
// C++ Interface: $MODULE$
//
// Description: 
//
//
// Author: Fabian Wenzel <f.wenzel@gmx.net>, (C) 2003
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PANOINTERACTIVECANVAS_H
#define PANOINTERACTIVECANVAS_H

#include "panocanvas.h"
#include <wx/dnd.h>

/**
@author Fabian Wenzel
*/

class panoFrame;

class panoDropTarget : public wxFileDropTarget{
public:
  panoDropTarget(panoFrame *frame);
  bool OnDropFiles(wxCoord x, wxCoord y, const wxArrayString& filenames);
private:
  panoFrame *p_frame;
};
  
class panoInteractiveCanvas : public panoCanvas
{
public:
  panoInteractiveCanvas(wxWindow* parent, int id, const wxPoint& position, const wxSize& size);
  void OnKeyDown(wxKeyEvent& event);
  void OnKeyUp(wxKeyEvent& event);
  void OnMouse(wxMouseEvent& event);
  void OnTimer(wxTimerEvent &event);
  void OnPaint(wxPaintEvent &event);
  void OnSize(wxSizeEvent &event);
  void setBoundaryMode(int boundarymode);
  void enableShowBoundaries(bool show);
  void enableUseBoundaries(bool use);
  void resetGivenBoundaries();
  bool anyKeysOrButtonsAreDown();
  void clearAllKeyAndButtonStates();

  bool getShowBoundaries()
  {
    return m_showboundaries;
  }
    
  bool getUseBoundaries()
  {
    return m_useboundaries;
  }

  void incrementPosition(CPosition increment);
  void updateStatusText();
  
  const CBoundaries &getGivenBoundaries() const
  {
    return m_givenboundaries;
  }

  void setGivenBoundaries(const CBoundaries &newboundaries)
  {
    m_givenboundaries = newboundaries;
  }
  
  ~panoInteractiveCanvas();
private:
  void showAllBoundaries();
  void showActiveBoundary();

  void showPanTiltLine(double degree, float red, float green, float blue, bool pan);
  void setBoundary(const wxPoint &position);

  int sign(const double &value) const;
  bool     m_leftbuttondown;
  bool     m_rightbuttondown;
  bool     m_zoomindown;
  bool     m_zoomoutdown;
  bool     m_leftdown;
  bool     m_rightdown;
  bool     m_updown;
  bool     m_downdown;
  bool     m_enableInertia;
  bool     m_toggleInertiaKeyDown;
  wxPoint  m_clickposition;
  wxPoint  m_diff;
  wxPoint  m_currentpos;
  wxTimer  m_timer;
  int      m_timerelapse;
  int      m_boundarymode;
  bool     m_showboundaries;
  bool     m_useboundaries;
  int      m_wheelrot;
  int      m_stickypandirection;
  int      m_stickytiltdirection;

  CPosition m_increment;

  CBoundaries m_givenboundaries;
  panoFrame *p_frame;

  DECLARE_EVENT_TABLE();
};

#endif
