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

module oh_mod

  !includes OH concentration field as well as the height information
  !for this field

  implicit none

  integer :: nxOH,nyOH,nzOH
  real, allocatable, dimension(:) :: lonOH,latOH,altOH
  real, allocatable, dimension(:,:,:,:) :: OH_hourly
  real, allocatable, dimension (:,:,:,:) :: OH_field
  real, dimension(2) :: memOHtime
  real, dimension(360,180,12) :: jrate_average
  real, dimension(360) :: lonjr
  real, dimension(180) :: latjr

  integer :: nxNH3,nyNH3,nzNH3,ntNH3
  real, allocatable, dimension(:) :: lonNH3,latNH3,altNH3
  double precision, allocatable, dimension(:) :: timeNH3
  real, allocatable, dimension (:,:,:) :: NH3LOSS_field

end module oh_mod
