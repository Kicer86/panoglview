/***************************************************************************
                          panocanvas.h  -  description
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

#ifndef BOUNDARIES_H
#define BOUNDARIES_H


//invalid degree
#define INVALIDDEG 6000.0

class CAngleRange{
public:
  CAngleRange() :
  m_min(INVALIDDEG),
  m_max(INVALIDDEG)
  {
  }

  CAngleRange(const double &min, const double &max) :
  m_min(min),
  m_max(max)
  {
  }

  bool inRange(const double &angle) const;
  CAngleRange getBorders(const CAngleRange &innerborders) const;

  bool validMin() const
  {
    return m_min != INVALIDDEG;
  }

  bool validMax() const
  {
    return m_max != INVALIDDEG;
  }

  bool validRange() const
  {
    return validMin() && validMax();
  }
  
  const double &getMin() const
  {
    return m_min;
  }

  const double &getMax() const
  {
    return m_max;
  }

  void setMin(const double &angle)
  {
    m_min = angle;
  }

  void setMax(const double &angle)
  {
    m_max = angle;
  }
  
protected:
  double m_min;
  double m_max;
};

class CPanRange : public CAngleRange{
public:
  CPanRange() : CAngleRange()
  {
  }

  CPanRange(const double &min, const double &max ):
  CAngleRange(min,max)
  {
  }
  
  bool inRange(const double &angle) const;
  CAngleRange getBorders(const CPanRange &innerrange) const;

  bool isWrapped() const
  {
    return m_max < m_min;
  }
};

class CBoundaries{
public:
  CBoundaries(const CPanRange &pan, const CAngleRange &tilt, const CAngleRange &fov);
  CBoundaries();

  const CPanRange &getPans() const
  {
    return m_pans;
  }

  const CAngleRange &getTilts() const
  {
    return m_tilts;
  }

  const CAngleRange &getFovs() const
  {
    return m_fovs;
  }
  
  bool inTiltRange(const double &tilt) const
  {
    return m_tilts.inRange(tilt);
  }
  
  bool inPanRange (const double &pan ) const
  {
    return m_pans.inRange(pan);
  }

  void setPanmin(const double &panmin)
  {
    m_pans.setMin(panmin);
  }

  void setPanmax(const double &panmax)
  {
    m_pans.setMax(panmax);
  }

  void setTiltmin(const double &tiltmin)
  {
    m_tilts.setMin(tiltmin);
  }

  void setTiltmax(const double &tiltmax)
  {
    m_tilts.setMax(tiltmax);
  }

  void setFovmin(const double &fovmin)
  {
    m_fovs.setMin(fovmin);
  }

  void setFovmax(const double &fovmax)
  {
    m_fovs.setMax(fovmax);
  }

  void getPanBorders (const double &panmin,  const double &panmax,  double *bordermin, double *bordermax) const;
  void getTiltBorders(const double &tiltmin, const double &tiltmax, double *bordermin, double *bordermax) const;

  
private:
  CAngleRange  m_tilts;
  CPanRange    m_pans;
  CAngleRange  m_fovs;
};


#endif
