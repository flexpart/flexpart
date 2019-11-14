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

real function psim(z,al)

  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION: CALCULATION OF THE STABILITY CORRECTION FUNCTION FOR   *
  !              MOMENTUM AS FUNCTION OF HEIGHT Z AND OBUKHOV SCALE     *
  !              HEIGHT L                                               *
  !                                                                     *
  !**********************************************************************

  use par_mod

  implicit none

  real :: z,al,zeta,x,a1,a2

  zeta=z/al
  if(zeta.le.0.) then
  ! UNSTABLE CASE
    x=(1.-15.*zeta)**0.25
    a1=((1.+x)/2.)**2
    a2=(1.+x**2)/2.
    psim=log(a1*a2)-2.*atan(x)+pi/2.
  else
  ! STABLE CASE
    psim=-4.7*zeta
  endif

end function psim
