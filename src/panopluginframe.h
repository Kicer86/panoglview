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
#ifndef PANOPLUGINFRAME_H
#define PANOPLUGINFRAME_H

#include <wx/frame.h>
#include <tuv/CImageByteRGB.h>
#include "panocanvas.h"

enum{
  id_snapshot=3,
};

/**
@author Fabian Wenzel
*/
class panoPluginFrame : public wxFrame
{
public:
  panoPluginFrame(const wxString &title, const wxPoint &position = wxDefaultPosition, const wxSize &size = wxDefaultSize, const wxSize &canvassize = wxDefaultSize);  
  panoCanvas *getCanvas()
  {
    return m_canvas;
  }
  ~panoPluginFrame();
  TUVision::CImageByteRGB *getSnapshot();
  void OnSnapshot(wxCommandEvent &event);
  TUVision::CImageByteRGB *m_snapshot;
private:
  panoCanvas *m_canvas;
  wxMutex     m_mutex;
  wxCondition m_condition;
  DECLARE_EVENT_TABLE()
};

#endif
