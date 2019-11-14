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

!-----------------------------------------------------------------------
function distance(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH.
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN DEGREES
  !rlon1    - REAL LONGITUDE OF POINT 1 IN DEGREES
  !rlat2    - REAL LATITUDE OF POINT 2 IN DEGREES
  !rlon2    - REAL LONGITUDE OF POINT 2 IN DEGREES
  !
  !   OUTPUT ARGUMENT LIST:
  !distance - REAL GREAT CIRCLE DISTANCE IN KILOMETERS
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real          :: rlat1,rlon1,rlat2,rlon2,distance
  real(kind=dp) :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp, dpr=180.0_dp/pi
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ((abs(rlat1-rlat2).lt.0.03).and. &
       (abs(rlon1-rlon2).lt.0.03)) then
    distance=0.
  else
    clat1=cos(real(rlat1,kind=dp)/dpr)
    slat1=sin(real(rlat1,kind=dp)/dpr)
    clat2=cos(real(rlat2,kind=dp)/dpr)
    slat2=sin(real(rlat2,kind=dp)/dpr)
    cdlon=cos(real((rlon1-rlon2),kind=dp)/dpr)
    crd=slat1*slat2+clat1*clat2*cdlon
    distance=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance
