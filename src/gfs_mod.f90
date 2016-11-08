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

!*******************************************************************************
!   Include file for calculation of particle trajectories (Program FLEXPART)   *
!                                                                              *
!       This file contains GFS specific parameters used in FLEXPART            *
!       Note that module name differs from file name.                          *
!       The makefile selects either this file, or ecmwf_mod.f90, depending     *
!       on target.                                                             *
!                                                                              *
!       Author: ESO                                                            *
!                                                                              *
!       2015                                                                   *
!                                                                              *
!*******************************************************************************

module wind_mod

  implicit none

  !*********************************************
  ! Maximum dimensions of the input mother grids
  !*********************************************

!  integer,parameter :: nxmax=721,nymax=361,nuvzmax=64,nwzmax=64,nzmax=64
  integer,parameter :: nxmax=361,nymax=181,nuvzmax=64,nwzmax=64,nzmax=64
  integer,parameter :: nxshift=0     ! for GFS or FNL

  !*********************************************
  ! Maximum dimensions of the nested input grids
  !*********************************************

  integer,parameter :: maxnests=0,nxmaxn=361,nymaxn=181

  ! nxmax,nymax        maximum dimension of wind fields in x and y
  !                    direction, respectively
  ! nuvzmax,nwzmax     maximum dimension of (u,v) and (w) wind fields in z
  !                    direction (for fields on eta levels)
  ! nzmax              maximum dimension of wind fields in z direction
  !                    for the transformed Cartesian coordinates
  ! nxshift            for global grids (in x), the grid can be shifted by
  !                    nxshift grid points, in order to accomodate nested
  !                    grids, and output grids overlapping the domain "boundary"
  !                    nxshift must not be negative; "normal" setting would be 0
end module wind_mod
