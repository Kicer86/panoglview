/***************************************************************************
                          panoapp.h  -  description
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

#ifndef PANOAPP_H
#define PANOAPP_H

#include <wx/wx.h>

/**
  *@author Fabian Wenzel
  */

class panoApp : public wxApp  {
public: 
  virtual bool OnInit();
};

DECLARE_APP(panoApp)

#endif
