/***************************************************************************
                          boundaries.cpp  -  description
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

#include "boundaries.h"

CBoundaries::CBoundaries(const CPanRange &pan, const CAngleRange &tilt, const CAngleRange &fov) :
m_pans(pan),
m_tilts(tilt),
m_fovs(fov)
{
}

CBoundaries::CBoundaries()
{
}


  // Pan can wrap, so this is a function that nicyfies code
bool CAngleRange::inRange(const double &angle) const
{
  if(m_max != INVALIDDEG && angle > m_max)
    return false;
  return !(m_min != INVALIDDEG && angle < m_min);
}

bool CPanRange::inRange(const double &pan) const
{
  if(m_min < m_max)
    return (pan > m_min) && (pan < m_max);
  else
    return ((pan > -180.0 && pan < m_max) || (pan < 180.0 && pan > m_min));
}

CAngleRange CPanRange::getBorders (const CPanRange &innerRange) const
{
  
  if( (innerRange.isWrapped()) == isWrapped() )
  {
    return CAngleRange(m_min - innerRange.getMin(),m_max - innerRange.getMax());
  } else {
    if(m_min < innerRange.getMin()){
      return CAngleRange(m_min - innerRange.getMin(),m_max - innerRange.getMax() + 360.0);
    } else {
      return CAngleRange(m_min - innerRange.getMin() - 360.0,m_max - innerRange.getMax());
    }
  }
}

CAngleRange CAngleRange::getBorders(const CAngleRange &innerRange) const
{
  return CAngleRange(m_min - innerRange.getMin(),m_max - innerRange.getMax());
}

