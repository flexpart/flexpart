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

subroutine readOHfield

!*****************************************************************************
!                                                                            *
! Reads the OH field into memory                                             *
!                                                                            *
! AUTHOR: R.L. Thompson, Nov 2014                                            *
!                                                                            *
! UPDATES:                                                                   *
!   03/2018 SEC: Converted original netCDF files to binary format            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
!                                                                            *
! path(numpath)              contains the path names                         *
! lonOH(nxOH)                longitude of OH fields                          *
! latOH(nyOH)                latitude of OH fields                           *
! altOH(nzOH)                altitude of OH fields                           *
! etaOH(nzOH)                eta-levels of OH fields                         *
! OH_field(nxOH,nyOH,nzOH,m) OH concentration (molecules/cm3)                *
!                                                                            *
!                                                                            *
!*****************************************************************************

  use oh_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i,j,k,l,ierr
  real, dimension(:), allocatable :: etaOH

!  real, parameter :: gasct=8.314   ! gas constant
!  real, parameter :: mct=0.02894   ! kg mol-1
!  real, parameter :: g=9.80665     ! m s-2
!  real, parameter :: lrate=0.0065  ! K m-1
  real, parameter :: scalehgt=7000. ! scale height in metres


  open(unitOH,file=trim(ohfields_path) &
       //'OH_FIELDS/OH_variables.bin',status='old', &
       form='UNFORMATTED', iostat=ierr, convert='little_endian')

  if(ierr.ne.0) then
    write(*,*) 'Cannot read binary OH fields in ',trim(ohfields_path)//'OH_FIELDS/OH_variables.bin'
    stop
  endif

  read(unitOH) nxOH
  read(unitOH) nyOH
  read(unitOH) nzOH
  write(*,*) nxOH,nyOH,nzOH

! allocate variables
  allocate(lonOH(nxOH))
  allocate(latOH(nyOH))
  allocate(etaOH(nzOH))
  allocate(altOH(nzOH))
  allocate(OH_field(nxOH,nyOH,nzOH,12))
  allocate(OH_hourly(nxOH,nyOH,nzOH,2))

  read(unitOH) (lonjr(i),i=1,360)
  read(unitOH) (latjr(i),i=1,180)
  read(unitOH) (((jrate_average(i,j,k),i=1,360),j=1,180),k=1,12)
  read(unitOH) (lonOH(i),i=1,nxOH)
  read(unitOH) (latOH(i),i=1,nyOH)
  read(unitOH) (lonOH(i),i=1,nxOH)

  read(unitOH) (altOH(i),i=1,nzOH)
  read(unitOH) ((((OH_field(i,j,k,l),i=1,nxOH),j=1,nyOH),k=1,nzOH),l=1,12)
  read(unitOH) ((((OH_hourly(i,j,k,l),i=1,nxOH),j=1,nyOH),k=1,nzOH),l=1,2)

end subroutine readOHfield

