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

module point_mod

  implicit none

  integer, allocatable, dimension (:) :: ireleasestart
  integer, allocatable, dimension (:) :: ireleaseend
  integer, allocatable, dimension (:) :: npart
  integer*2, allocatable, dimension (:) :: kindz

  real,allocatable, dimension (:) :: xpoint1
  real,allocatable, dimension (:) :: xpoint2
  real,allocatable, dimension (:) :: ypoint1
  real,allocatable, dimension (:) :: ypoint2
  real,allocatable, dimension (:) :: zpoint1
  real,allocatable, dimension (:) :: zpoint2

  real,allocatable, dimension (:,:) :: xmass
  real,allocatable, dimension (:) :: rho_rel

end module point_mod
