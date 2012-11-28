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

function raerod (l,ust,z0)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the aerodynamical resistance ra from ground up to href  *
  !                                                                            *
  !     AUTHOR: Matthias Langer, modified by Andreas Stohl (6 August 1993)     *
  !                                                                            *
  !     Literature:                                                            *
  !     [1]  Hicks/Baldocchi/Meyers/Hosker/Matt (1987), A Preliminary          *
  !             Multiple Resistance Routine for Deriving Dry Deposition        *
  !             Velocities from Measured Quantities.                           *
  !             Water, Air and Soil Pollution 36 (1987), pp.311-330.           *
  !     [2]  Scire/Yamartino/Carmichael/Chang (1989),                          *
  !             CALGRID: A Mesoscale Photochemical Grid Model.                 *
  !             Vol II: User's Guide. (Report No.A049-1, June, 1989)           *
  !                                                                            *
  !     Variable list:                                                         *
  !     L     = Monin-Obukhov-length [m]                                       *
  !     ust   = friction velocity [m/sec]                                      *
  !     z0    = surface roughness length [m]                                   *
  !     href  = reference height [m], for which deposition velocity is         *
  !             calculated                                                     *
  !                                                                            *
  !     Constants:                                                             *
  !     karman    = von Karman-constant (~0.4)                                 *
  !     ramin = minimum resistence of ra (1 s/m)                               *
  !                                                                            *
  !     Subprograms and functions:                                             *
  !     function psih (z/L)                                                    *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  real :: l,psih,raerod,ust,z0

  raerod=(alog(href/z0)-psih(href,l)+psih(z0,l))/(karman*ust)

end function raerod
