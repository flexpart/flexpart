!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann, Ignacio Pisso                                      *
!                                                                     *
! This file is part of FLEXPART-NorESM                                *
!                                                                     *
! FLEXPART-NorESM is free software: you can redistribute it           *
! and/or modify                                                       *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART-NorESM is distributed in the hope that it will be useful,  *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART-NorESM.                                         *
!  If not, see <http://www.gnu.org/licenses/>.                        * 
!**********************************************************************

module noresm_variables
  !*****************************************************************************
  !                                                                            *
  !     This module defines some support variables used to read netcdf files   *
  !     content
  !                                                                            *
  !     Author:                                                                *
  !     M. Cassiani  2016                                                      *
  !                                                                            *
  !*****************************************************************************
 
    implicit none
    save
    
    integer, allocatable, dimension(:) :: istart,icount

    !--------------0D fields---------------------------------    
    double precision dumvar
    real dumvar_real
    integer dumvar_int
    !c------------ 1D fields --------------------------------
    double precision, allocatable, dimension(:) :: dumarray1D
    real, allocatable, dimension(:) :: dumarray1D_real
    integer, allocatable, dimension(:) :: dumarray1D_int
     character, allocatable, dimension (:) :: dumarray1D_char
    !c------------ 2D fields --------------------------------
    double precision, allocatable, dimension(:,:) :: dumarray2D
    real, allocatable, dimension(:,:) :: dumarray2D_real
    integer, allocatable, dimension(:,:) :: dumarray2D_int
    character, allocatable, dimension (:,:) :: dumarray2D_char    
    !c------------ 3D fields ---------------------------------
    double precision, allocatable, dimension(:,:,:) :: dumarray3D
    real, allocatable, dimension(:,:,:) :: dumarray3D_real
    integer, allocatable, dimension(:,:,:) :: dumarray3D_int 
    !c------------ 4D fields ---------------------------------
    double precision, allocatable, dimension(:,:,:,:) :: dumarray4D
    real, allocatable, dimension(:,:,:,:) :: dumarray4D_real
    integer, allocatable, dimension(:,:,:,:) :: dumarray4D_int
    
end module noresm_variables
