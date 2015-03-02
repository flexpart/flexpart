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

module unc_mod

  implicit none

  real,allocatable ,dimension (:,:,:,:,:,:,:) :: gridunc
  real,allocatable, dimension (:,:,:,:,:,:,:) :: griduncn
  real,allocatable, dimension (:,:,:,:,:,:) :: drygridunc
  real,allocatable, dimension (:,:,:,:,:,:) :: drygriduncn
  real,allocatable, dimension (:,:,:,:,:,:) :: wetgridunc
  real,allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn

! For sum of individual contributions, used for the MPI version
  real,allocatable, dimension (:,:,:,:,:,:) :: drygridunc0
  real,allocatable, dimension (:,:,:,:,:,:) :: drygriduncn0
  real,allocatable, dimension (:,:,:,:,:,:) :: wetgridunc0
  real,allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn0

  real,allocatable, dimension (:,:,:,:,:) :: init_cond

end module unc_mod
