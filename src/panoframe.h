/***************************************************************************
                          panoframe.h  -  description
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

#ifndef PANOFRAME_H
#define PANOFRAME_H

#include <wx/wx.h>
#include "panointeractivecanvas.h"

/**
  *@author Fabian Wenzel
  */

enum{
  ID_QUIT = 10,
  ID_OPENIMAGE,
  ID_FULLSCREEN,
  ID_OPENPROJECT,
  ID_SAVEPROJECT,
  // Respect oder of ID's (see OnBoundary)
  ID_QUITBOUNDARY,
  ID_LEFTBOUNDARY,
  ID_RIGHTBOUNDARY,
  ID_TOPBOUNDARY,
  ID_BOTTOMBOUNDARY,
  ID_USEBOUNDARIES,
  ID_SHOWBOUNDARIES,
  ID_RESETBOUNDARIES,
  ID_RECORD,
  ID_PLAY,
  ID_STOP,
  ID_PLAYTIMER,
  ID_ABOUT
};

  
class panoFrame : public wxFrame  {
public:
  panoFrame(const wxString &title, const wxPoint &position = wxDefaultPosition, const wxSize &size = wxDefaultSize);
  panoCanvas *getCanvas()
  {
    return m_canvas;
  }
  void OnOpenImage(wxCommandEvent &event);
  void OnQuit(wxCommandEvent &event);
  void OnFullScreen(wxCommandEvent &event);
  void OnBoundary(wxCommandEvent &event);
  void OnUseBoundaries(wxCommandEvent &event);
  void OnShowBoundaries(wxCommandEvent &event);
  void OnResetBoundaries(wxCommandEvent &event);
  void OnOpenProject(wxCommandEvent &event);
  void OnSaveProject(wxCommandEvent &event);
  void OnRecord(wxCommandEvent &event);
  void OnStop(wxCommandEvent &event);
  void OnPlay(wxCommandEvent &event);
  void OnPlayRecord(wxTimerEvent &event);
  void OnAbout(wxCommandEvent &event);
  void openArgumentFile(const wxString &filename);
  void quitFullscreen();
private:
  void openProject(const wxString &filename);
  void openImage(const wxString &filename);
  
  DECLARE_EVENT_TABLE();
  panoInteractiveCanvas *m_canvas;
  wxString m_imagename;
  wxTimer m_playrecordtimer;
  std::vector< CPosition > m_movie;
  std::vector< CPosition >::iterator m_frame;
  bool m_recording;
};

#endif
