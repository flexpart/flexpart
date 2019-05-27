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

module outg_mod

  use par_mod, only: dep_prec, sp

  implicit none

  real,allocatable, dimension (:) :: outheight
  real,allocatable, dimension (:) :: outheighthalf
  real,allocatable, dimension (:,:) :: oroout
  real,allocatable, dimension (:,:) :: orooutn
  real,allocatable, dimension (:,:) :: area
  real,allocatable, dimension (:,:) :: arean
  real,allocatable, dimension (:,:,:) :: volume
  real,allocatable, dimension (:,:,:) :: volumen
  real,allocatable, dimension (:,:,:) :: areaeast
  real,allocatable, dimension (:,:,:) :: areanorth
  real,allocatable, dimension (:,:,:) :: densityoutgrid
  real,allocatable, dimension (:,:,:) :: densitydrygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor_drygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor3d
  real,allocatable, dimension (:,:,:) :: grid
  real(sp),allocatable, dimension (:,:) :: wetgrid
  real(sp),allocatable, dimension (:,:) :: drygrid
  real,allocatable, dimension (:,:,:) :: gridsigma
  real(dep_prec),allocatable, dimension (:,:) :: drygridsigma
  real(dep_prec),allocatable, dimension (:,:) :: wetgridsigma
  real,allocatable, dimension (:) :: sparse_dump_r
  real,allocatable, dimension (:) :: sparse_dump_u
  integer,allocatable, dimension (:) :: sparse_dump_i

end module outg_mod
