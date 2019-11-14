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

module flux_mod

  ! flux eastward, westward, northward, southward, upward and downward
  ! fluxes of all species and all ageclasses
  ! areaeast,areanorth [m2] side areas of each grid cell

  implicit none

  real,allocatable, dimension (:,:,:,:,:,:,:) :: flux

  !1 fluxw west - east
  !2 fluxe east - west
  !3 fluxs south - north
  !4 fluxn north - south
  !5 fluxu upward
  !6 fluxd downward
  !real,allocatable, dimension (:,:,:) :: areanorth
  !real,allocatable, dimension (:,:,:) :: areaeast

end module flux_mod
