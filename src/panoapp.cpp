/***************************************************************************
                          panoapp.cpp  -  description
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

#include "panoapp.h"
#include "panoframe.h"

#include <wx/image.h>

bool panoApp::OnInit()
{
  wxInitAllImageHandlers();
  panoFrame *panoframe = new panoFrame(_("PanoramaViewer"),wxDefaultPosition,wxSize(640,480));
    
  panoframe->Show(TRUE);
  
  if (argc == 2)
    // We have been passed an argument, either an image or a text file
    panoframe->openArgumentFile(wxString(argv[1]));

  SetTopWindow(panoframe);
  return TRUE;
}

IMPLEMENT_APP(panoApp)
