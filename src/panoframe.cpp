/***************************************************************************
                          panoframe.cpp  -  description
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

#include <wx/wxprec.h>
#include <wx/wx.h>

#include "panoframe.h"

#include <wx/config.h>
#include <wx/fileconf.h>

#include <wx/sizer.h>
#include <wx/filename.h>
#include <wx/tokenzr.h>

#include <fstream>
#include <locale.h>
#include "smallicon.xpm"

using namespace std;

BEGIN_EVENT_TABLE(panoFrame,wxFrame)
  EVT_MENU(ID_OPENIMAGE,panoFrame::OnOpenImage)    
  EVT_MENU(ID_QUIT,panoFrame::OnQuit)
  EVT_MENU(ID_FULLSCREEN,panoFrame::OnFullScreen)
  EVT_MENU(ID_OPENPROJECT,panoFrame::OnOpenProject)
  EVT_MENU(ID_SAVEPROJECT,panoFrame::OnSaveProject)
  EVT_MENU(ID_LEFTBOUNDARY,panoFrame::OnBoundary)
  EVT_MENU(ID_RIGHTBOUNDARY,panoFrame::OnBoundary)
  EVT_MENU(ID_TOPBOUNDARY,panoFrame::OnBoundary)
  EVT_MENU(ID_BOTTOMBOUNDARY,panoFrame::OnBoundary)
  EVT_MENU(ID_USEBOUNDARIES,panoFrame::OnUseBoundaries)
  EVT_MENU(ID_SHOWBOUNDARIES,panoFrame::OnShowBoundaries)
  EVT_MENU(ID_RESETBOUNDARIES,panoFrame::OnResetBoundaries)
  EVT_MENU(ID_RECORD,panoFrame::OnRecord)
  EVT_MENU(ID_PLAY,panoFrame::OnPlay)
  EVT_MENU(ID_STOP,panoFrame::OnStop)
  EVT_MENU(ID_ABOUT,panoFrame::OnAbout)
  EVT_TIMER(ID_PLAYTIMER,panoFrame::OnPlayRecord)
END_EVENT_TABLE()

panoFrame::panoFrame(const wxString &title, const wxPoint &position, const wxSize &size) :
wxFrame((wxWindow *) NULL, -1, title,position,size),
m_playrecordtimer(this,ID_PLAYTIMER),
m_imagename(wxT("")),
m_movie()
{
  wxMenu *filemenu = new wxMenu();
  filemenu->Append(ID_OPENIMAGE,_("&Open Image\tCtrl-O"),_("Opens a Panorama Image"));
  filemenu->Append(ID_OPENPROJECT,_("&Open Project"),_("Opens a Project File"));
  filemenu->Append(ID_SAVEPROJECT,_("&Save Project\tCtrl-S"),_("Saves a Project File"));
  filemenu->AppendSeparator();
  filemenu->Append(ID_QUIT,_("&Quit\tAlt-F4"),_("Quits the Application"));
  
  SetIcon(smallicon);
  
  wxMenuBar *menubar = new wxMenuBar();
  menubar->Append(filemenu,_("&File"));
  
  wxMenu *viewmenu = new wxMenu();
  viewmenu->Append(ID_FULLSCREEN,_("Fullscreen\tCtrl-F"),_("Changes to Fullscreen Mode"),wxITEM_CHECK);
  
  menubar->Append(viewmenu,_("View"));

  wxMenu *boundarymenu = new wxMenu();
  boundarymenu->Append(ID_LEFTBOUNDARY,_("Left Boundary\tCtrl-L"),_("Sets left boundary"));
  boundarymenu->Append(ID_RIGHTBOUNDARY,_("Right Boundary\tCtrl-R"),_("Sets right boundary"));
  boundarymenu->Append(ID_TOPBOUNDARY,_("Top Boundary\tCtrl-T"),_("Sets top boundary"));
  boundarymenu->Append(ID_BOTTOMBOUNDARY,_("Bottom Boundary\tCtrl-B"),_("Sets bottom boundary"));
  boundarymenu->AppendSeparator();
  boundarymenu->Append(ID_SHOWBOUNDARIES,_("Show Boundaries\tCtrl-S"),_("Show Boundaries"),wxITEM_CHECK);
  boundarymenu->Check(ID_SHOWBOUNDARIES,true);
  boundarymenu->Append(ID_USEBOUNDARIES,_("Use Boundaries\tCtrl-U"),_("Use Boundaries"),wxITEM_CHECK);
  boundarymenu->Check(ID_USEBOUNDARIES,false);
  boundarymenu->Append(ID_RESETBOUNDARIES,_("Reset Boundaries\tCtrl-N"),_("Reset Boundaries"));

  menubar->Append(boundarymenu,_("Boundaries"));  

  wxMenu *moviemenu = new wxMenu();
  moviemenu->Append(ID_RECORD,_("Record Movie"),_("Starts recording positions"));
  moviemenu->Append(ID_PLAY, _("Play Movie"),_("Starts playing positions"));
  moviemenu->Enable(ID_PLAY,false);
  
  moviemenu->AppendSeparator();
  moviemenu->Append(ID_STOP,_("Stop"),_("Stops recording or playing"));
  moviemenu->Enable(ID_STOP,false);

  menubar->Append(moviemenu,_("Movie"));
  
  wxMenu *aboutmenu = new wxMenu();
  aboutmenu->Append(ID_ABOUT,_("About"),_("Shows information about the program"));
  
  menubar->Append(aboutmenu,_("Help"));
  
  
  SetMenuBar(menubar);

  m_canvas = new panoInteractiveCanvas(this,-1,wxDefaultPosition,wxSize(50,50));
  wxSizer    *topsizer = new wxBoxSizer(wxVERTICAL);

  CreateStatusBar(3);
}

void panoFrame::openImage(const wxString &filename)
{
  wxImage image;
  if(image.LoadFile(filename)){
     m_canvas->createPanorama(image);
     m_imagename = filename;
  } else
    wxMessageBox(_("Could not load image"));
}

void panoFrame::OnOpenImage(wxCommandEvent &event)
{
  wxString imagename = wxFileSelector(_("Open PanoramaImage"),wxEmptyString,wxEmptyString,wxT("png"),_("All known image types|*.bmp;*.png;*.jpg;*.gif;*.pcx;*.pnm;*.tif;*.xpm|\
BMP files (*.bmp)|bmp|\
GIF files (*.gif)|*.gif|\
PNG files (*.png)|*png|\
JPEG files (*.jpg)|*.jpg|\
PCX files (*.pcx)|*.pcx|\
PNM files (*.pnm)|*.pnm|\
TIFF files (*.tif)|*.tif|\
XPM files (*.xpm)|*.xpm|\
All files (*.*)|*.*"),wxOPEN|wxFILE_MUST_EXIST);
  if (imagename !=wxT(""))
    openImage(imagename);
}

void panoFrame::OnQuit(wxCommandEvent &event)
{
  Close(TRUE);
}

void panoFrame::OnFullScreen(wxCommandEvent &event)
{
  ShowFullScreen(GetMenuBar()->IsChecked(ID_FULLSCREEN));
}

void panoFrame::OnOpenProject(wxCommandEvent &event)
{
  wxString filename = wxFileSelector(_("Open Project"),wxEmptyString,wxEmptyString,wxT("paf"),_("Panorama files (*.paf)|*.paf|All files (*.*)|*.*"),wxOPEN|wxFILE_MUST_EXIST);
  if (filename !=wxT("")){
    openProject(filename);
  }
}

void panoFrame::openProject(const wxString &filename)
{
  wxString imagename;
  CPosition initialpos;
  bool tmp;

  CBoundaries boundaries = m_canvas->getGivenBoundaries();
  CBoundaries img_boundaries;

  wxString configname = filename;
  
  if(wxPathOnly(filename) == wxT(""))  
    configname = wxGetCwd() + wxFileName::GetPathSeparator() + filename;

  char *locOld = setlocale(LC_NUMERIC,"C"); // or LC_ALL
  
  wxFileConfig config(wxEmptyString,wxEmptyString,
                       configname,wxEmptyString,wxCONFIG_USE_LOCAL_FILE);
  double angle;

  img_boundaries.setTiltmin(-90.0);
  img_boundaries.setTiltmax(90.0);
  img_boundaries.setPanmin(-180.0);
  img_boundaries.setPanmax(180.0);

  if(config.Read(wxT("Image Minimum Tilt"),&angle))
    if (img_boundaries.inTiltRange(angle))
      img_boundaries.setTiltmin(angle);

  if(config.Read(wxT("Image Maximum Tilt"),&angle))
    if (img_boundaries.inTiltRange(angle))
      img_boundaries.setTiltmax(angle);

  if(config.Read(wxT("Image Minimum Pan"),&angle))
    if (img_boundaries.inPanRange(angle))
      img_boundaries.setPanmin(angle);

  if(config.Read(wxT("Image Maximum Pan"),&angle))
    if (img_boundaries.inPanRange(angle))
      img_boundaries.setPanmax(angle);

  if(config.Read(wxT("Minimum Tilt"),&angle))
    boundaries.setTiltmin(angle);
    
  if(config.Read(wxT("Maximum Tilt"),&angle))
    boundaries.setTiltmax(angle);
    
  if(config.Read(wxT("Minimum Pan"),&angle))
    boundaries.setPanmin(angle);
    
  if(config.Read(wxT("Maximum Pan"),&angle))
    boundaries.setPanmax(angle);
    
  if(config.Read(wxT("Minimum FOV"),&angle))
    boundaries.setFovmin(angle);
    
  if(config.Read(wxT("Maximum FOV"),&angle))
    boundaries.setFovmax(angle);
    
  if(config.Read(wxT("Initial Tilt"),&angle))
    initialpos.setTilt(angle);
    
  if(config.Read(wxT("Initial Pan"),&angle))
    initialpos.setPan(angle);
    
  if(config.Read(wxT("Initial FOV"),&angle))
    initialpos.setFov(angle);
    
  config.Read(wxT("Panorama Image"),&imagename,wxT(""));
 
    
  tmp = GetMenuBar()->IsChecked(ID_USEBOUNDARIES);
  config.Read(wxT("Use Boundaries"),&tmp);
  GetMenuBar()->Check(ID_USEBOUNDARIES,tmp);
  m_canvas->enableUseBoundaries(tmp);
  
  tmp = GetMenuBar()->IsChecked(ID_SHOWBOUNDARIES);
  config.Read(wxT("Show Boundaries"),&tmp);
  GetMenuBar()->Check(ID_SHOWBOUNDARIES,tmp);
  m_canvas->enableShowBoundaries(tmp);
  
  if(wxPathOnly(filename) != wxT(""))
    wxSetWorkingDirectory(wxPathOnly(filename));

  m_canvas->setGivenBoundaries(boundaries);
  m_canvas->setImageBoundaries(img_boundaries);

  openImage(imagename);
  
  m_canvas->setPosition(initialpos);

  int numberofframes;
  
  // We have a movie here
  if(config.Read(wxT("Movie Frames"),&numberofframes) && numberofframes > 0 ){
    m_movie.clear();
    for(int i=0;i<numberofframes;i++){
      wxString position;
      config.Read(wxString::Format(wxT("Frame %d"),i),&position);
      wxStringTokenizer tkz(position,wxT(":"));
      if(tkz.CountTokens() != 3) {
        wxMessageBox(wxString::Format(_("Could not read movie frame %d"),i),_("Information"));
        numberofframes = i;
      } else {
        double pan;
        double tilt;
        double fov;
        if(tkz.GetNextToken().ToDouble(&pan)  &&
           tkz.GetNextToken().ToDouble(&tilt) &&
           tkz.GetNextToken().ToDouble(&fov))
          m_movie.push_back(CPosition(pan,tilt,fov));
        else {
          numberofframes = i;
          wxMessageBox(wxString::Format(_("Could not read movie frame %d"),i),_("Information"));
        }
      }
    }
    GetMenuBar()->Enable(ID_PLAY,true);
  }

  setlocale(LC_NUMERIC, locOld);
}

void panoFrame::openArgumentFile(const wxString &filename)
{
  if (filename.Right(4) == wxString(wxT(".paf")))
    openProject(filename);
  else
    openImage(filename);
}

void panoFrame::OnSaveProject(wxCommandEvent &event)
{
  wxString filename = wxFileSelector(_("Save Project"),wxEmptyString,wxEmptyString,wxT("paf"),_("Panorama files (*.paf)|*.paf|All files (*.*)|*.*"),wxOVERWRITE_PROMPT|wxSAVE);
  if (filename != wxT("")){
    if (wxFileExists(filename))
      wxRemoveFile(filename);
    wxString imagename;
    
    CPosition initialpos   = m_canvas->getPosition();
    CBoundaries boundaries = m_canvas->getGivenBoundaries();
    CBoundaries img_boundaries = m_canvas->getImageBoundaries();


    char *locOld = setlocale(LC_NUMERIC, "C"); // or LC_ALL

    wxFileConfig config(wxEmptyString,wxEmptyString,
                         filename,wxEmptyString,wxCONFIG_USE_LOCAL_FILE);
    if(img_boundaries.getTilts().validMin())
      config.Write(wxT("Image Minimum Tilt"),img_boundaries.getTilts().getMin());
    if(img_boundaries.getTilts().validMax())
      config.Write(wxT("Image Maximum Tilt"),img_boundaries.getTilts().getMax());
    if(img_boundaries.getPans().validMin())
      config.Write(wxT("Image Minimum Pan"),img_boundaries.getPans().getMin());
    if(img_boundaries.getPans().validMax())
      config.Write(wxT("Image Maximum Pan"),img_boundaries.getPans().getMax());
    if(boundaries.getTilts().validMin())
      config.Write(wxT("Minimum Tilt"),boundaries.getTilts().getMin());
    if(boundaries.getTilts().validMax())
      config.Write(wxT("Maximum Tilt"),boundaries.getTilts().getMax());
    if(boundaries.getPans().validMin())
      config.Write(wxT("Minimum Pan"),boundaries.getPans().getMin());
    if(boundaries.getPans().validMax())
      config.Write(wxT("Maximum Pan"),boundaries.getPans().getMax());
    if(boundaries.getFovs().validMin())
      config.Write(wxT("Minimum FOV"),boundaries.getFovs().getMin());
    if(boundaries.getFovs().validMax())
      config.Write(wxT("Maximum FOV"),boundaries.getFovs().getMax());

    config.Write(wxT("Initial Tilt"),initialpos.getTilt());
    config.Write(wxT("Initial Pan"),initialpos.getPan());
    config.Write(wxT("Initial FOV"),initialpos.getFov());

    wxFileName imagefilename(m_imagename);
    wxFileName projectfilename(filename);

    imagefilename.MakeRelativeTo(projectfilename.GetPath()+wxFileName::GetPathSeparator());

    config.Write(wxT("Panorama Image"),imagefilename.GetFullPath());
    config.Write(wxT("Use Boundaries"),m_canvas->getUseBoundaries());
    config.Write(wxT("Show Boundaries"),m_canvas->getShowBoundaries());

    if(m_movie.size())
    {
      config.Write(wxT("Movie Frames"),(int) m_movie.size());
      int i=0;
      for(m_frame = m_movie.begin(); m_frame != m_movie.end();++m_frame,i++)
      {
        config.Write(wxString::Format(wxT("Frame %d"),i),wxString::Format(wxT("%f:%f:%f"),m_frame->getPan(),m_frame->getTilt(),m_frame->getFov()));
      }
    }
    setlocale(LC_NUMERIC, locOld);
  }
}

void panoFrame::quitFullscreen()
{
  if(GetMenuBar()->IsChecked(ID_FULLSCREEN)){
    GetMenuBar()->Check(ID_FULLSCREEN,false);
    wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_FULLSCREEN);
    AddPendingEvent(event);
  }
}

void panoFrame::OnBoundary(wxCommandEvent &event)
{
  // Call Canvas with ID (shifted to start with 0)
  // QUIT_BOUNDARY is mapped to -1
  m_canvas->setBoundaryMode(event.GetId() - ID_LEFTBOUNDARY);
  GetMenuBar()->Check(ID_USEBOUNDARIES,false);
  m_canvas->enableUseBoundaries(false);
}

void panoFrame::OnUseBoundaries(wxCommandEvent &event)
{
  m_canvas->enableUseBoundaries(GetMenuBar()->IsChecked(ID_USEBOUNDARIES));
}

void panoFrame::OnShowBoundaries(wxCommandEvent &event)
{
  m_canvas->enableShowBoundaries(GetMenuBar()->IsChecked(ID_SHOWBOUNDARIES));
}

void panoFrame::OnResetBoundaries(wxCommandEvent &event)
{
  m_canvas->resetGivenBoundaries();
  GetMenuBar()->Check(ID_USEBOUNDARIES,false);
  m_canvas->enableUseBoundaries(false);
}

void panoFrame::OnRecord(wxCommandEvent &event)
{
  m_movie.clear();
  m_playrecordtimer.Start(20);
  m_recording = true;
  GetMenuBar()->Enable(ID_STOP  ,true );
  GetMenuBar()->Enable(ID_RECORD,false);
  GetMenuBar()->Enable(ID_PLAY  ,false);
  GetMenuBar()->Enable(ID_OPENPROJECT,false);
  GetMenuBar()->Enable(ID_SAVEPROJECT,false);
  GetMenuBar()->Enable(ID_OPENIMAGE,false);
  GetMenuBar()->Enable(ID_LEFTBOUNDARY,false);
  GetMenuBar()->Enable(ID_RIGHTBOUNDARY,false);
  GetMenuBar()->Enable(ID_TOPBOUNDARY,false);
  GetMenuBar()->Enable(ID_BOTTOMBOUNDARY,false);
}

void panoFrame::OnStop(wxCommandEvent &event)
{
  m_playrecordtimer.Stop();
  if(m_movie.size())
    GetMenuBar()->Enable(ID_PLAY,true);
  GetMenuBar()->Enable(ID_STOP,false);
  GetMenuBar()->Enable(ID_RECORD,true);
  GetMenuBar()->Enable(ID_OPENPROJECT,true);
  GetMenuBar()->Enable(ID_SAVEPROJECT,true);
  GetMenuBar()->Enable(ID_OPENIMAGE,true);
  GetMenuBar()->Enable(ID_LEFTBOUNDARY,true);
  GetMenuBar()->Enable(ID_RIGHTBOUNDARY,true);
  GetMenuBar()->Enable(ID_TOPBOUNDARY,true);
  GetMenuBar()->Enable(ID_BOTTOMBOUNDARY,true);
}

void panoFrame::OnPlay(wxCommandEvent &event)
{
  m_playrecordtimer.Start(20);
  m_recording = false;
  m_frame = m_movie.begin();
  GetMenuBar()->Enable(ID_STOP,  true );
  GetMenuBar()->Enable(ID_RECORD,false);
  GetMenuBar()->Enable(ID_PLAY,  false);
  GetMenuBar()->Enable(ID_OPENPROJECT,false);
  GetMenuBar()->Enable(ID_SAVEPROJECT,false);
  GetMenuBar()->Enable(ID_OPENIMAGE,false);
  GetMenuBar()->Enable(ID_LEFTBOUNDARY,false);
  GetMenuBar()->Enable(ID_RIGHTBOUNDARY,false);
  GetMenuBar()->Enable(ID_TOPBOUNDARY,false);
  GetMenuBar()->Enable(ID_BOTTOMBOUNDARY,false);

}

void panoFrame::OnPlayRecord(wxTimerEvent &event)
{
  if(m_recording)
  {
    m_movie.push_back(m_canvas->getPosition());
  } else {
    m_canvas->setPosition(*m_frame++);
    m_canvas->Refresh();
    if(m_frame == m_movie.end()){
      m_playrecordtimer.Stop();
      GetMenuBar()->Enable(ID_STOP,false);
      GetMenuBar()->Enable(ID_RECORD,true);
      GetMenuBar()->Enable(ID_OPENPROJECT,true);
      GetMenuBar()->Enable(ID_SAVEPROJECT,true);
      GetMenuBar()->Enable(ID_OPENIMAGE,true);
      GetMenuBar()->Enable(ID_LEFTBOUNDARY,true);
      GetMenuBar()->Enable(ID_RIGHTBOUNDARY,true);
      GetMenuBar()->Enable(ID_TOPBOUNDARY,true);
      GetMenuBar()->Enable(ID_BOTTOMBOUNDARY,true);
    }
  }
}

void panoFrame::OnAbout(wxCommandEvent &event)
{
  wxMessageBox(_("GL Panorama Viewer V 0.2.2\n\
(c) 2003 Fabian Wenzel\n\
sourceforge.net/projects/hugin\n\
\n\
Keys:\n\
Up/Down/Left/Right to pan\n\
Zoom with + and -\n\
\n\
Mouse:\n\
  Left-click pans\n\
  Right-click zooms\n\
  Comma (',') toggles inertia\n\
\n\
Control-F enters fullscreen,\n\
ESC exits fullscreen."),_("Information"));
}

