! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
  !********************************************************************
  !                                                                   *
  !     INPUT:                                                        *
  !                                                                   *
  !     ps      surface pressure [Pa]                                 *
  !     tsurf   surface temperature [K]                               *
  !     tdsurf  surface dew point [K]                                 *
  !     tlev    temperature first model level [K]                     *
  !     ustar   scale velocity [m/s]                                  *
  !     hf      surface sensible heat flux [W/m2]                     *
  !     akm     ECMWF vertical discretization parameter               *
  !     bkm     ECMWF vertical discretization parameter               *
  !     plev                                                          *
  !     metdata_format format of metdata (ecmwf/gfs)                  *
  !                                                                   *
  !********************************************************************

  use par_mod
  use class_gribfile

  implicit none

  integer :: metdata_format
  real :: akm(nwzmax),bkm(nwzmax)
  real :: ps,tsurf,tdsurf,tlev,ustar,hf,e,ew,tv,rhoa,plev
  real :: ak1,bk1,theta,thetastar


  e=ew(tdsurf)                           ! vapor pressure
  tv=tsurf*(1.+0.378*e/ps)               ! virtual temperature
  rhoa=ps/(r_air*tv)                      ! air density
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
