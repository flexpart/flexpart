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

! Changes to the routines by A. Stohl
! xi,xi0,eta,eta0 are double precision variables to avoid problems
! at poles

module cmapf_mod

  use par_mod, only: dp

  implicit none
  private

  public :: cc2gll, cll2xy, cgszll, cxy2ll, stlmbr, stcm2p

  real,parameter :: rearth=6371.2, almst1=.9999999

  real,parameter :: pi=3.14159265358979
  real,parameter :: radpdg=pi/180., dgprad=180./pi

contains

subroutine cc2gll (strcmp, xlat,xlong, ue,vn, ug,vg)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9), xlat, xlong, ue, vn, ug, vg

  real(kind=dp) :: xpolg,ypolg,along,slong,clong,rot

  along = cspanf( xlong - strcmp(2), -180., 180.)
  if (xlat.gt.89.985) then
  !*  North polar meteorological orientation: "north" along prime meridian
    rot = - strcmp(1) * along + xlong - 180.
  elseif (xlat.lt.-89.985) then
  !*  South polar meteorological orientation: "north" along prime meridian
    rot = - strcmp(1) * along - xlong
  else
    rot = - strcmp(1) * along
  endif
  slong = sin( radpdg * rot )
  clong = cos( radpdg * rot )
  xpolg = slong * strcmp(5) + clong * strcmp(6)
  ypolg = clong * strcmp(5) - slong * strcmp(6)
  ug = ypolg * ue + xpolg * vn
  vg = ypolg * vn - xpolg * ue
  return
end subroutine cc2gll

subroutine ccrvll (strcmp, xlat,xlong, gx,gy)
  !*  Written on 9/20/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real(kind=dp) :: xpolg,ypolg,temp,along,slong,clong,ctemp, curv
  real :: strcmp(9), xlat, xlong, gx, gy

  along = cspanf( xlong - strcmp(2), -180., 180.)
  slong = sin( radpdg * strcmp(1) * along)
  clong = cos( radpdg * strcmp(1) * along)
  xpolg = - slong * strcmp(5) + clong * strcmp(6)
  ypolg = clong * strcmp(5) + slong * strcmp(6)
  temp = sin(radpdg * xlat)
  ctemp = cos(radpdg * xlat)
  curv = (strcmp(1) - temp) / ctemp / rearth
  gx = curv * xpolg
  gy = curv * ypolg
  return
end subroutine ccrvll

subroutine ccrvxy (strcmp, x,y, gx,gy)
  !*  Written on 9/20/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9), x, y, gx, gy
  real(kind=dp) :: xpolg,ypolg,temp,ymerc,efact,curv

  temp = strcmp(1) * strcmp(7) /rearth
  xpolg = strcmp(6) + temp * (strcmp(3) - x)
  ypolg = strcmp(5) + temp * (strcmp(4) - y)
  temp = sqrt ( xpolg ** 2 + ypolg ** 2 )
  if (temp.gt.0.) then
    ymerc = - log( temp) /strcmp(1)
    efact = exp(ymerc)
    curv = ( (strcmp(1) - 1.d0) * efact + &
         (strcmp(1) + 1.d0) / efact ) &
         * .5d0 / rearth
    gx = xpolg * curv / temp
    gy = ypolg * curv / temp
  else
    if (abs(strcmp(1)) .eq. 1.) then
      gx = 0.
      gy = 0.
    else
      gx = 1./rearth
      gy = 1./rearth
    endif
  endif
  return
end subroutine ccrvxy

subroutine cg2cll (strcmp, xlat,xlong, ug,vg, ue,vn)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real(kind=dp) :: xpolg,ypolg,along,slong,clong,rot
  real :: strcmp(9), xlat, xlong, ug, vg, ue, vn

  along = cspanf( xlong - strcmp(2), -180., 180.)
  if (xlat.gt.89.985) then
  !*  North polar meteorological orientation: "north" along prime meridian
    rot = - strcmp(1) * along + xlong - 180.
  elseif (xlat.lt.-89.985) then
  !*  South polar meteorological orientation: "north" along prime meridian
    rot = - strcmp(1) * along - xlong
  else
    rot = - strcmp(1) * along
  endif
  slong = sin( radpdg * rot )
  clong = cos( radpdg * rot )
  xpolg = slong * strcmp(5) + clong * strcmp(6)
  ypolg = clong * strcmp(5) - slong * strcmp(6)
  ue = ypolg * ug - xpolg * vg
  vn = ypolg * vg + xpolg * ug
  return
end subroutine cg2cll

subroutine cg2cxy (strcmp, x,y, ug,vg, ue,vn)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9) , x, y, ug, vg, ue, vn

  real :: clong, radial, rot, slong, xlat, xlong
  real(kind=dp) :: xpolg,ypolg,temp,xi0,eta0,xi,eta

  xi0 = ( x - strcmp(3) ) * strcmp(7) / rearth
  eta0 = ( y - strcmp(4) ) * strcmp(7) /rearth
  xi = xi0 * strcmp(5) - eta0 * strcmp(6)
  eta = eta0 * strcmp(5) + xi0 * strcmp(6)
  radial = 2. * eta - strcmp(1) * (xi*xi + eta*eta)
  if (radial.gt.strcmp(8)) then
  !*  Case north of 89 degrees. Meteorological wind direction definition
  !*      changes.
    call cnxyll(strcmp, xi,eta, xlat,xlong)
  !*  North polar meteorological orientation: "north" along prime meridian
    rot = strcmp(1) * (xlong - strcmp(2)) - xlong - 180.
    slong = - sin( radpdg * rot )
    clong = cos( radpdg * rot )
    xpolg = slong * strcmp(5) + clong * strcmp(6)
    ypolg = clong * strcmp(5) - slong * strcmp(6)
  else if (radial.lt.strcmp(9)) then
  !*  Case south of -89 degrees. Meteorological wind direction definition
  !*      changes.
    call cnxyll(strcmp, xi,eta, xlat,xlong)
  !*  South polar meteorological orientation: "north" along prime meridian
    rot = strcmp(1) * (xlong - strcmp(2)) + xlong
    slong = - sin( radpdg * rot )
    clong = cos( radpdg * rot )
    xpolg = slong * strcmp(5) + clong * strcmp(6)
    ypolg = clong * strcmp(5) - slong * strcmp(6)
  else
  !* Normal case. Meteorological direction of wind related to true north.
    xpolg = strcmp(6) - strcmp(1) * xi0
    ypolg = strcmp(5) - strcmp(1) * eta0
    temp = sqrt ( xpolg ** 2 + ypolg ** 2 )
    xpolg = xpolg / temp
    ypolg = ypolg / temp
  end if
  ue = ( ypolg * ug - xpolg * vg )
  vn = ( ypolg * vg + xpolg * ug )
  return
end subroutine cg2cxy

real function cgszll (strcmp, xlat,xlong)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9), xlat, xlong

  real(kind=dp) :: slat,ymerc,efact

  if (xlat .gt. 89.985) then
  !* Close to north pole
    if (strcmp(1) .gt. 0.9999) then
  !* and to gamma == 1.
      cgszll = 2. * strcmp(7)
      return
    endif
    efact = cos(radpdg * xlat)
    if (efact .le. 0.) then
      cgszll = 0.
      return
    else
      ymerc = - log( efact /(1. + sin(radpdg * xlat)))
    endif
  else if (xlat .lt. -89.985) then
  !* Close to south pole
    if (strcmp(1) .lt. -0.9999) then
  !* and to gamma == -1.0
      cgszll = 2. * strcmp(7)
      return
    endif
    efact = cos(radpdg * xlat)
    if (efact .le. 0.) then
      cgszll = 0.
      return
    else
      ymerc = log( efact /(1. - sin(radpdg * xlat)))
    endif
  else
  slat = sin(radpdg * xlat)
  ymerc = log((1. + slat) / (1. - slat))/2.
  !efact = exp(ymerc)
  !cgszll = 2. * strcmp(7) * exp (strcmp(1) * ymerc)
  !c			 / (efact + 1./efact)
  endif
  cgszll = strcmp(7) * cos(radpdg * xlat) * exp(strcmp(1) *ymerc)
  return
end function cgszll

real function cgszxy (strcmp, x,y)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9) , x, y
  real(kind=dp) :: ymerc,efact, radial, temp
  real(kind=dp) :: xi0,eta0,xi,eta


  xi0 = ( x - strcmp(3) ) * strcmp(7) / rearth
  eta0 = ( y - strcmp(4) ) * strcmp(7) /rearth
  xi = xi0 * strcmp(5) - eta0 * strcmp(6)
  eta = eta0 * strcmp(5) + xi0 * strcmp(6)
  radial = 2. * eta - strcmp(1) * (xi*xi + eta*eta)
  efact = strcmp(1) * radial
  if (efact .gt. almst1) then
    if (strcmp(1).gt.almst1) then
      cgszxy = 2. * strcmp(7)
    else
      cgszxy = 0.
    endif
    return
  endif
  if (abs(efact) .lt. 1.e-2) then
    temp = (efact / (2. - efact) )**2
    ymerc = radial / (2. - efact) * (1.    + temp * &
         (1./3. + temp * &
         (1./5. + temp * &
         (1./7. ))))
  else
    ymerc = - log( 1. - efact ) /2. /strcmp(1)
  endif
  if (ymerc .gt. 6.) then
    if (strcmp(1) .gt. almst1) then
      cgszxy = 2. * strcmp(7)
    else
      cgszxy = 0.
    endif
  else if (ymerc .lt. -6.) then
    if (strcmp(1) .lt. -almst1) then
      cgszxy = 2. * strcmp(7)
    else
      cgszxy = 0.
    endif
  else
    efact = exp(ymerc)
    cgszxy = 2. * strcmp(7) * exp (strcmp(1) * ymerc) &
         / (efact + 1./efact)
  endif
  return
end function cgszxy

subroutine cll2xy (strcmp, xlat,xlong, x,y)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  implicit none

  real :: strcmp(9) , xlat, xlong, x, y, xi, eta

  call cnllxy(strcmp, xlat,xlong, xi,eta)
  x = strcmp(3) + rearth/strcmp(7) * &
       (xi * strcmp(5) + eta * strcmp(6) )
  y = strcmp(4) + rearth/strcmp(7) * &
       (eta * strcmp(5) - xi * strcmp(6) )
  return
end subroutine cll2xy

subroutine cnllxy (strcmp, xlat,xlong, xi,eta)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL
  !  main transformation routine from latitude-longitude to
  !  canonical (equator-centered, radian unit) coordinates

  use par_mod, only: dp

  implicit none

  real :: strcmp(9), xlat, xlong, xi, eta, &
       gdlong, sndgam, csdgam, rhog1
  real(kind=dp) :: gamma
  real(kind=dp) :: dlong,dlat,slat,mercy,gmercy

  gamma = strcmp(1)
  dlat = xlat
  dlong = cspanf(xlong - strcmp(2), -180., 180.)
  dlong = dlong * radpdg
  gdlong = gamma * dlong
  if (abs(gdlong) .lt. .01) then
  !  Code for gamma small or zero.  This avoids round-off error or divide-
  !  by zero in the case of mercator or near-mercator projections.
    gdlong = gdlong * gdlong
    sndgam = dlong * (1. - 1./6. * gdlong * &
         (1. - 1./20. * gdlong * &
         (1. - 1./42. * gdlong )))
    csdgam = dlong * dlong * .5 * &
         (1. - 1./12. * gdlong * &
         (1. - 1./30. * gdlong * &
         (1. - 1./56. * gdlong )))
  else
  ! Code for moderate values of gamma
    sndgam = sin (gdlong) /gamma
    csdgam = (1. - cos(gdlong) )/gamma /gamma
  endif
  slat = sin(radpdg * dlat)
  if ((slat .ge. almst1) .or. (slat .le. -almst1)) then
    eta = 1./strcmp(1)
    xi = 0.
    return
  endif
  mercy = .5 * log( (1. + slat) / (1. - slat) )
  gmercy = gamma * mercy
  if (abs(gmercy) .lt. .001) then
  !  Code for gamma small or zero.  This avoids round-off error or divide-
  !  by zero in the case of mercator or near-mercator projections.
    rhog1 = mercy * (1. - .5 * gmercy * &
         (1. - 1./3. * gmercy * &
         (1. - 1./4. * gmercy ) ) )
  else
  ! Code for moderate values of gamma
    rhog1 = (1. - exp(-gmercy)) / gamma
  endif
  eta = rhog1 + (1. - gamma * rhog1) * gamma * csdgam
  xi = (1. - gamma * rhog1 ) * sndgam
end subroutine cnllxy

subroutine cnxyll (strcmp, xi,eta, xlat,xlong)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL
  !  main transformation routine from canonical (equator-centered,
  !  radian unit) coordinates

  use par_mod, only: dp

  implicit none

  real :: strcmp(9), xlat, xlong, odist
  real(kind=dp) :: gamma,temp,arg1,arg2,ymerc,along,gxi,cgeta
  real(kind=dp) :: xi,eta

  gamma = strcmp(1)
  !  Calculate equivalent mercator coordinate
  odist = xi*xi + eta*eta
  arg2 = 2. * eta - gamma * (xi*xi + eta*eta)
  arg1 = gamma * arg2
  ! Change by A. Stohl to avoid problems close to the poles
  ! if (arg1 .ge. almst1) then
  !  distance to north (or south) pole is zero (or imaginary ;) )
  ! xlat = sign(90.,strcmp(1))
  ! xlong = strcmp(2)
  ! return
  ! endif
  if (abs(arg1) .lt. .01) then
  !  Code for gamma small or zero.  This avoids round-off error or divide-
  !  by zero in the case of mercator or near-mercator projections.
    temp = (arg1 / (2. - arg1) )**2
    ymerc = arg2 / (2. - arg1) * (1.    + temp * &
         (1./3. + temp * &
         (1./5. + temp * &
         (1./7. ))))
  else
  ! Code for moderate values of gamma
    ymerc = - log ( 1. - arg1 ) /2. / gamma
  endif
  ! Convert ymerc to latitude
  temp = exp( - abs(ymerc) )
  xlat = sign(atan2((1. - temp) * (1. + temp), 2. * temp), ymerc)
  ! Find longitudes
  gxi = gamma*xi
  cgeta = 1. - gamma * eta
  if ( abs(gxi) .lt. .01*cgeta ) then
  !  Code for gamma small or zero.  This avoids round-off error or divide-
  !  by zero in the case of mercator or near-mercator projections.
    temp = ( gxi /cgeta )**2
    along = xi / cgeta * (1.    - temp * &
         (1./3. - temp * &
         (1./5. - temp * &
         (1./7.   ))))
  else
  ! Code for moderate values of gamma
    along = atan2( gxi, cgeta) / gamma
  endif
  xlong = sngl(strcmp(2) + dgprad * along)
  xlat = xlat * dgprad
  return
end subroutine cnxyll

subroutine cpolll (strcmp, xlat,xlong, enx,eny,enz)
  !*  Written on 11/23/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real(kind=dp) :: xpolg,ypolg,along,slong,clong,rot
  real :: strcmp(9), xlat, xlong, enx, eny, enz, clat

  along = cspanf( xlong - strcmp(2), -180., 180.)
  rot = - strcmp(1) * along
  slong = sin( radpdg * rot )
  clong = cos( radpdg * rot )
  xpolg = slong * strcmp(5) + clong * strcmp(6)
  ypolg = clong * strcmp(5) - slong * strcmp(6)
  clat = cos(radpdg * xlat)
  enx = clat * xpolg
  eny = clat * ypolg
  enz = sin(radpdg * xlat)
  return
end subroutine cpolll

subroutine cpolxy (strcmp, x,y, enx,eny,enz)
  !*  Written on 11/26/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real :: strcmp(9) , x, y, enx, eny, enz
  real(kind=dp) :: xpol,ypol,temp,xi0,eta0,xi,eta,radial
  real(kind=dp) :: temp2,ymerc,arg,oarg,clat

  xi0 = ( x - strcmp(3) ) * strcmp(7) / rearth
  eta0 = ( y - strcmp(4) ) * strcmp(7) /rearth
  xi = xi0 * strcmp(5) - eta0 * strcmp(6)
  eta = eta0 * strcmp(5) + xi0 * strcmp(6)
  radial = 2. * eta -  strcmp(1) * (xi*xi + eta*eta)
  temp = strcmp(1) * radial
  if (temp .ge. 1.) then
    enx = 0.
    eny = 0.
    enz = sign(1.,strcmp(1))
    return
  endif
  if (abs(temp).lt.1.e-2) then
    temp2 = (temp / (2. - temp))**2
    ymerc = radial / (2. - temp) * (1. + temp2 * &
         (1./3. + temp2 * &
         (1./5. + temp2 * &
         (1./7.))))
  else
    ymerc = -.5 * log(1. - temp) / strcmp(1)
  endif
  arg = exp( ymerc )
  oarg = 1./arg
  clat = 2./(arg + oarg)
  enz = (arg - oarg) * clat /2.
  temp = clat / sqrt(1. - temp)
  xpol = - xi * strcmp(1) * temp
  ypol = (1. - eta * strcmp(1) ) * temp
  enx = xpol * strcmp(5) + ypol * strcmp(6)
  eny = ypol * strcmp(5) - xpol * strcmp(6)
  return
end subroutine cpolxy

real function cspanf (value, begin, end)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL
  !* real function cspanf returns a value in the interval (begin,end]
  !* which is equivalent to value, mod (end - begin).  It is used to
  !* reduce periodic variables to a standard range.  It adjusts for the
  !* behavior of the mod function which provides positive results for
  !* positive input, and negative results for negative input
  !* input:
  !*       value - real number to be reduced to the span
  !*       begin - first value of the span
  !*       end   - last value of the span
  !* returns:
  !*       the reduced value
  !* examples:
  !*      along = cspanf(xlong, -180., +180.)
  !*      dir  = cspanf(angle, 0., 360.)

  implicit none

  real :: first,last, value, begin, end, val

  first = min(begin,end)
  last = max(begin,end)
  val = mod( value - first , last - first)
  if ( val .le. 0.) then
    cspanf = val + last
  else
    cspanf = val + first
  endif
  return
end function cspanf

subroutine cxy2ll (strcmp, x,y, xlat,xlong)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  use par_mod, only: dp

  implicit none

  real(kind=dp) :: xi0,eta0,xi,eta
  real :: strcmp(9), x, y, xlat, xlong

  xi0 = ( x - strcmp(3) ) * strcmp(7) / rearth
  eta0 = ( y - strcmp(4) ) * strcmp(7) /rearth
  xi = xi0 * strcmp(5) - eta0 * strcmp(6)
  eta = eta0 * strcmp(5) + xi0 * strcmp(6)
  call cnxyll(strcmp, xi,eta, xlat,xlong)
  xlong = cspanf(xlong, -180., 180.)
  return
end subroutine cxy2ll

real function eqvlat (xlat1,xlat2)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  implicit none

  real :: xlat1, xlat2, x, ssind, sinl1, sinl2, al1, al2, tau

  ssind(x) = sin (radpdg*x)
  sinl1 = ssind (xlat1)
  sinl2 = ssind (xlat2)
  if (abs(sinl1 - sinl2) .gt. .001) then
    al1 = log((1. - sinl1)/(1. - sinl2))
    al2 = log((1. + sinl1)/(1. + sinl2))
  else
  !  Case lat1 near or equal to lat2
    tau = - (sinl1 - sinl2)/(2. - sinl1 - sinl2)
    tau = tau*tau
    al1  = 2. / (2. - sinl1 - sinl2) * (1.    + tau * &
         (1./3. + tau * &
         (1./5. + tau * &
         (1./7.))))
    tau =   (sinl1 - sinl2)/(2. + sinl1 + sinl2)
    tau = tau*tau
    al2  = -2. / (2. + sinl1 + sinl2) * (1.    + tau * &
         (1./3. + tau * &
         (1./5. + tau * &
         (1./7.))))
  endif
  eqvlat = asin((al1 + al2) / (al1 - al2))/radpdg
  return
end function eqvlat

subroutine stcm1p(strcmp, x1,y1, xlat1,xlong1, &
       xlatg,xlongg, gridsz, orient)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  implicit none

  integer :: k
  real :: strcmp(9), x1, y1, xlat1, xlong1, turn, orient, &
       xlatg, xlongg, gridsz, x1a, y1a

  do k=3,4
    strcmp (k) = 0.
  enddo
    turn = radpdg * (orient - strcmp(1) * &
         cspanf(xlongg - strcmp(2), -180., 180.) )
  strcmp (5) = cos (turn)
  strcmp (6) = - sin (turn)
  strcmp (7) = 1.
  strcmp (7) = gridsz * strcmp(7) &
       / cgszll(strcmp, xlatg, strcmp(2))
  call cll2xy (strcmp, xlat1,xlong1, x1a,y1a)
  strcmp(3) = strcmp(3) + x1 - x1a
  strcmp(4) = strcmp(4) + y1 - y1a
  return
end subroutine stcm1p

subroutine stcm2p(strcmp, x1,y1, xlat1,xlong1, &
       x2,y2, xlat2,xlong2)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  implicit none

  real :: strcmp(9), x1, y1, xlat1, xlong1, &
       x2, y2, xlat2, xlong2

  integer :: k
  real    :: x1a, y1a, x2a, y2a, den, dena

  do k=3,6
    strcmp (k) = 0.
  enddo
  strcmp (5) = 1.
  strcmp (7) = 1.
  call cll2xy (strcmp, xlat1,xlong1, x1a,y1a)
  call cll2xy (strcmp, xlat2,xlong2, x2a,y2a)
  den = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
  dena = sqrt( (x1a - x2a)**2 + (y1a - y2a)**2 )
  strcmp(5) = ((x1a - x2a)*(x1 - x2) + (y1a - y2a) * (y1 - y2)) &
       /den /dena
  strcmp(6) = ((y1a - y2a)*(x1 - x2) - (x1a - x2a) * (y1 - y2)) &
       /den /dena
  strcmp (7) = strcmp(7) * dena / den
  call cll2xy (strcmp, xlat1,xlong1, x1a,y1a)
  strcmp(3) = strcmp(3) + x1 - x1a
  strcmp(4) = strcmp(4) + y1 - y1a
  return
end subroutine stcm2p

!*  General conformal map routines for meteorological modelers
!*  written on 3/31/94 by

!* Dr. Albion Taylor
!* NOAA / OAR / ARL                  phone: (301) 713-0295 x 132
!* rm. 3151, 1315 east-west highway  fax:   (301) 713-0119
!* silver spring, md 20910           e-mail: adtaylor@arlrisc.ssmc.noaa.gov

!*  subroutine stlmbr (strcmp, tnglat, clong)
!*    This routine initializes the map structure array strcmp to
!*    the form of a specific map projection
!*  inputs:
!*    tnglat - the latitude at which the projection will be tangent
!*             to the earth.  +90. For north polar stereographic,
!*             -90. for south polar stereographic, 0. For mercator,
!*             and other values for lambert conformal.
!*             -90 <= tnglat <= 90.
!*    clong -  a longitude in the region under consideration.  Longitudes
!*             between clong-180. and clong+180.  Will be mapped in one
!*             connected region
!*  outputs:
!*    strcmp - a 9-value map structure array for use with subsequent
!*             calls to the coordinate transform routines.
!*
!*  real function eqvlat (xlat1,xlat2)
!*    This function is provided to assist in finding the tangent latitude
!*    equivalent to the 2-reference latitude specification in the legend
!*    of most lambert conformal maps.  If the map specifies "scale
!*    1:xxxxx true at 40n and 60n", then eqvlat(40.,60.) will return the
!*    equivalent tangent latitude.
!*  inputs:
!*    xlat1,xlat2:  the two latitudes specified in the map legend
!*  returns:
!*    the equivalent tangent latitude
!*  example:  call stlmbr(strcmp, eqvlat(40.,60.), 90.)

!*  subroutine stcm2p (strcmp, x1,y1, xlat1,xlong1,
!*          x2,y2, xlat2,xlong2)
!*  subroutine stcm1p (strcmp, x1,y1, xlat1,xlong1,
!*          xlatg,xlongg, gridsz, orient)
!*    These routines complete the specification of the map structure
!*    array by conforming the map coordinates to the specifications
!*    of a particular grid.  Either stcm1p or stcm2p must be called,
!*    but not both
!*  inputs:
!*    strcmp - a 9-value map structure array, set to a particular map
!*             form by a previous call to stlmbr
!*    for stcm2p:
!*      x1,y1, x2,y2 - the map coordinates of two points on the grid
!*      xlat1,xlong1, xlat2,xlong2 - the geographic coordinates of the
!*             same two points
!*    for stcm1p:
!*      x1,y1 - the map coordinates of one point on the grid
!*      xlat1,xlong1 - the geographic coordinates of the same point
!*      xlatg,xlongg - latitude and longitude of reference point for
!*             gridsz and orientation specification.
!*      gridsz - the desired grid size in kilometers, at xlatg,xlongg
!*      orient - the angle, with respect to north, of a y-grid line, at
!*             the point xlatg,xlongg
!*  outputs:
!*    strcmp - a 9-value map structure array, fully set for use by
!*             other subroutines in this system

!*  subroutine cll2xy (strcmp, xlat,xlong, x,y)
!*  subroutine cxy2ll (strcmp, x,y, xlat,xlong)
!*     these routines convert between map coordinates x,y
!*     and geographic coordinates xlat,xlong
!*  inputs:
!*     strcmp(9) - 9-value map structure array
!*     for cll2xy:  xlat,xlong - geographic coordinates
!*     for cxy2ll:  x,y - map coordinates
!*  outputs:
!*     for cll2xy:  x,y - map coordinates
!*     for cxy2ll:  xlat,xlong - geographic coordinates

!*  subroutine cc2gxy (strcmp, x,y, ue,vn, ug,vg)
!*  subroutine cg2cxy (strcmp, x,y, ug,vg, ue,vn)
!*  subroutine cc2gll (strcmp, xlat,xlong, ue,vn, ug,vg)
!*  subroutine cg2cll (strcmp, xlat,xlong, ug,vg, ue,vn)
!*     These subroutines convert vector wind components from
!*     geographic, or compass, coordinates, to map or
!*     grid coordinates.  The site of the wind to be
!*     converted may be given either in geographic or
!*     map coordinates.  Wind components are all in kilometers
!*     per hour, whether geographic or map coordinates.
!*  inputs:
!*    strcmp(9) - 9-value map structure array
!*    for cc2gxy and cg2cxy:  x,y        -  map coordinates of site
!*    for cc2gll and cg2cll:  xlat,xlong -  geographic coordinates of site
!*    for cc2gxy and cc2gll:  ue,vn - east and north wind components
!*    for cg2cxy and cg2cll:  ug,vg - x- and y- direction wind components
!*  outputs:
!*    for cc2gxy and cc2gll:  ug,vg - x- and y- direction wind components
!*    for cg2cxy and cg2cll:  ue,vn - east and north wind components

!*  subroutine ccrvxy (strcmp, x, y,       gx,gy)
!*  subroutine ccrvll (strcmp, xlat,xlong, gx,gy)
!*    These subroutines return the curvature vector (gx,gy), as referenced
!*    to map coordinates, induced by the map transformation.  When
!*    non-linear terms in wind speed are important, a "geodesic" force
!*    should be included in the vector form [ (u,u) g - (u,g) u ] where the
!*    inner product (u,g) is defined as ux*gx + uy*gy.
!*  inputs:
!*    strcmp(9) - 9-value map structure array
!*    for ccrvxy:  x,y        -  map coordinates of site
!*    for ccrvll:  xlat,xlong -  geographic coordinates of site
!*  outputs:
!*    gx,gy       - vector coefficients of curvature, in units radians
!*                  per kilometer

!*  real function cgszll (strcmp, xlat,xlong)
!*  real function cgszxy (strcmp, x,y)
!*    These functions return the size, in kilometers, of each unit of
!*    motion in map coordinates (grid size).  The grid size at any
!*    location depends on that location; the position may be given in
!*    either map or geographic coordinates.
!*  inputs:
!*    strcmp(9) - 9-value map structure array
!*    for cgszxy:  x,y        -  map coordinates of site
!*    for cgszll:  xlat,xlong -  geographic coordinates of site
!*  returns:
!*    gridsize in kilometers at given site.

!*  subroutine cpolxy (strcmp, x,y, enx,eny,enz)
!*  subroutine cpolll (strcmp, xlat,xlong, enx,eny,enz)
!*    These subroutines provide 3-d vector components of a unit vector
!*    in the direction of the north polar axis.  When multiplied
!*    by twice the rotation rate of the earth (2 * pi/24 hr), the
!*    vertical component yields the coriolis factor.
!*  inputs:
!*    strcmp(9) - 9-value map structure array
!*    for cpolxy:  x,y        -  map coordinates of site
!*    for cpolll:  xlat,xlong -  geographic coordinates of site
!*  returns:
!*    enx,eny,enz the direction cosines of a unit vector in the
!*    direction of the rotation axis of the earth

!*  subroutine cnllxy (strcmp, xlat,xlong, xi,eta)
!*  subroutine cnxyll (strcmp, xi,eta, xlat,xlong)
!*    These subroutines perform the underlying transformations from
!*    geographic coordinates to and from canonical (equator centered)
!*    coordinates.  They are called by cxy2ll and cll2xy, but are not
!*    intended to be called directly

!*  real function cspanf (value, begin, end)
!*    This function assists other routines in providing a longitude in
!*    the proper range.  It adds to value whatever multiple of
!*    (end - begin) is needed to return a number begin < cspanf <= end

subroutine stlmbr(strcmp, tnglat, xlong)
  !*  Written on 3/31/94 by Dr. Albion Taylor  NOAA / OAR / ARL

  implicit none

  real :: strcmp(9), tnglat, xlong

  real :: eta, xi

  strcmp(1) = sin(radpdg * tnglat)
  !*  gamma = sine of the tangent latitude
  strcmp(2) = cspanf( xlong, -180., +180.)
  !* lambda_0 = reference longitude
  strcmp(3) = 0.
  !* x_0 = x- grid coordinate of origin (xi,eta) = (0.,0.)
  strcmp(4) = 0.
  !* y_0 = y-grid coordinate of origin (xi,eta) = (0.,0.)
  strcmp(5) = 1.
  !* Cosine of rotation angle from xi,eta to x,y
  strcmp(6) = 0.
  !* Sine of rotation angle from xi,eta to x,y
  strcmp(7) = rearth
  !* Gridsize in kilometers at the equator
  call cnllxy(strcmp, 89.,xlong, xi,eta)
  strcmp(8) = 2. * eta - strcmp(1) * eta * eta
  !* Radial coordinate for 1 degree from north pole
  call cnllxy(strcmp, -89.,xlong, xi,eta)
  strcmp(9) = 2. * eta - strcmp(1) * eta * eta
  !* Radial coordinate for 1 degree from south pole
  return
end subroutine stlmbr

end module cmapf_mod
