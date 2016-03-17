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

function psih (z,l)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the stability correction term                           *
  !                                                                            *
  !     AUTHOR: Matthias Langer, adapted by Andreas Stohl (6 August 1993)      *
  !             Update: G. Wotawa, 11 October 1994                             *
  !                                                                            *
  !     Literature:                                                            *
  !     [1] C.A.Paulson (1970), A Mathematical Representation of Wind Speed    *
  !           and Temperature Profiles in the Unstable Atmospheric Surface     *
  !           Layer. J.Appl.Met.,Vol.9.(1970), pp.857-861.                     *
  !                                                                            *
  !     [2] A.C.M. Beljaars, A.A.M. Holtslag (1991), Flux Parameterization over*
  !           Land Surfaces for Atmospheric Models. J.Appl.Met. Vol. 30,pp 327-*
  !           341                                                              *
  !                                                                            *
  !     Variables:                                                             *
  !     L     = Monin-Obukhov-length [m]                                       *
  !     z     = height [m]                                                     *
  !     zeta  = auxiliary variable                                             *
  !                                                                            *
  !     Constants:                                                             *
  !     eps   = 1.2E-38, SUN-underflow: to avoid division by zero errors       *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  real :: psih,x,z,zeta,l
  real,parameter :: a=1.,b=0.667,c=5.,d=0.35,eps=1.e-20

  if ((l.ge.0).and.(l.lt.eps)) then
    l=eps
  else if ((l.lt.0).and.(l.gt.(-1.*eps))) then
    l=-1.*eps
  endif

  if ((log10(z)-log10(abs(l))).lt.log10(eps)) then
    psih=0.
  else
    zeta=z/l
    if (zeta.gt.0.) then
      psih = - (1.+0.667*a*zeta)**(1.5) - b*(zeta-c/d)*exp(-d*zeta) &
           - b*c/d + 1.
    else
      x=(1.-16.*zeta)**(.25)
      psih=2.*log((1.+x*x)/2.)
    end if
  end if

end function psih
