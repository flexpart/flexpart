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

  integer,parameter :: nxmax=361,nymax=181,nuvzmax=152,nwzmax=152,nzmax=152 !ECMWF new 
  integer,parameter :: nxshift=359

! 
  !*********************************************
  ! Maximum dimensions of the nested input grids
  !*********************************************

!  integer,parameter :: maxnests=0,nxmaxn=351,nymaxn=351 !ECMWF


end module wind_mod
