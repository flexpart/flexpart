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

real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,akm,bkm,plev,metdata_format)

  !********************************************************************
  !                                                                   *
  !                       Author: G. WOTAWA                           *
  !                       Date:   1994-06-27                          *
  !                                                                   *
  !     This program calculates Obukhov scale height from surface     *
  !     meteorological data and sensible heat flux.                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !  Update: A. Stohl, 2000-09-25, avoid division by zero by          *
  !  setting ustar to minimum value                                   *
  !  CHANGE: 17/11/2005 Caroline Forster NCEP GFS version             *
  !                                                                   *
  !   Unified ECMWF and GFS builds                                    *
  !   Marian Harustak, 12.5.2017                                      *
  !     - Merged obukhov and obukhov_gfs into one routine using       *
  !       if-then for meteo-type dependent code                       *
  !                                                                   *
  !   Don Morton, 13.10.2017                                          *
  !     - Adding documentation to explain the merger performed by     *
  !       Harustek.                                                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !     INPUT:                                                        *
  !                                                                   *
  !     ps             surface pressure [Pa]                          *
  !     tsurf          surface temperature [K]                        *
  !     tdsurf         surface dew point [K]                          *
  !     tlev           temperature first model level [K]              *
  !     ustar          scale velocity [m/s]                           *
  !     hf             surface sensible heat flux [W/m2]              *
  !                                                                   *
  !    ... akm and bkm are array input args used only for ECMWF ...   *
  !     akm            ECMWF vertical discretization parameter        *
  !     bkm            ECMWF vertical discretization parameter        *
  !                                                                   *
  !    ... plev is an input arg used only for GFS ...                 *
  !     plev                                                          *
  !                                                                   *
  !     metdata_format format of metdata (ecmwf/gfs)                  *
  !                                                                   *
  !********************************************************************


!!-----  Documentation added by Don Morton, 13 Oct 2017 ------- 
!
!    This function was originally represented as two separate functions,
!
!  For ECMWF:
!  real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,akm,bkm)
!
!  For GFS:
!  real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,plev)
!
!  I'm guessing on this, but it appears that the GFS version was written 
!  by Caroline Forster on 17/11/2005, by modifying the original ECMWF 
!  version.  The original akm and bkm array arguments were removed, and 
!  replaced with a scalar plev. 
!
!  In the case of ECMWF this value of plev is calculated as a function of
!  the akm and bkm arrays, but in GFS it's simply passed in as a parameter.
!
!
  use par_mod
  use class_gribfile   ! Module of GRIB-specific routines and data

  implicit none

  integer :: metdata_format
  real :: akm(nwzmax),bkm(nwzmax)
  real :: ps,tsurf,tdsurf,tlev,ustar,hf,e,ew,tv,rhoa,plev
  real :: ak1,bk1,theta,thetastar


  e=ew(tdsurf)                           ! vapor pressure
  tv=tsurf*(1.+0.378*e/ps)               ! virtual temperature
  rhoa=ps/(r_air*tv)                      ! air density

  ! For ECMWF calculate plev as function of akm and bkm
  ! For GFS, plev is provided as an input argument
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    ak1=(akm(1)+akm(2))/2.
    bk1=(bkm(1)+bkm(2))/2.
    plev=ak1+bk1*ps                        ! Pressure level 1
  end if

  theta=tlev*(100000./plev)**(r_air/cpa) ! potential temperature
  if (ustar.le.0.) ustar=1.e-8
  thetastar=hf/(rhoa*cpa*ustar)           ! scale temperature
  if(abs(thetastar).gt.1.e-10) then
     obukhov=theta*ustar**2/(karman*ga*thetastar)
  else
     obukhov=9999                        ! zero heat flux
  endif
  if (obukhov.gt. 9999.) obukhov= 9999.
  if (obukhov.lt.-9999.) obukhov=-9999.

end function obukhov
