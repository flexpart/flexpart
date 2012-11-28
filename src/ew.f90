!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

real function ew(x)

  !****************************************************************
  !SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
  !NACH DER GOFF-GRATCH-FORMEL.
  !****************************************************************

  implicit none

  real :: x, y, a, c, d

  ew=0.
  if(x.le.0.) stop 'sorry: t not in [k]'
  y=373.16/x
  a=-7.90298*(y-1.)
  a=a+(5.02808*0.43429*alog(y))
  c=(1.-(1./y))*11.344
  c=-1.+(10.**c)
  c=-1.3816*c/(10.**7)
  d=(1.-y)*3.49149
  d=-1.+(10.**d)
  d=8.1328*d/(10.**3)
  y=a+c+d
  ew=101324.6*(10.**y)       ! Saettigungsdampfdruck in Pa

end function ew
