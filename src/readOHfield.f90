! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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

