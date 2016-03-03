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
!        This file contains ECMWF specific parameters used in FLEXPART         *
!        Note that module name differs from file name.                         *
!        The makefile selects either this file, or gfs_mod.f90, depending      *
!        on target.                                                            *
!                                                                              *
!        Author: ESO                                                           *
!                                                                              *
!        2015                                                                  *
!                                                                              *
!*******************************************************************************

module wind_mod

  implicit none

  !*********************************************
  ! Maximum dimensions of the input mother grids
  !*********************************************

!  integer,parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92,nzmax=92 !ECMWF new 
  integer,parameter :: nxmax=361,nymax=181,nuvzmax=138,nwzmax=138,nzmax=138 !ECMWF new 
  integer,parameter :: nxshift=359
!  integer,parameter :: nxmax=361,nymax=181,nuvzmax=138,nwzmax=138,nzmax=138 !test BUG 
!  integer,parameter :: nxshift=0

! 
  !*********************************************
  ! Maximum dimensions of the nested input grids
  !*********************************************

  integer,parameter :: maxnests=1,nxmaxn=361,nymaxn=181
!  integer,parameter :: maxnests=1,nxmaxn=86,nymaxn=31


end module wind_mod
