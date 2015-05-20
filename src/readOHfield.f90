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

  include 'netcdf.inc'

  character(len=150) :: thefile
  character(len=2) :: mm
  integer :: nid,ierr,xid,yid,zid,vid,m
  real, dimension(:), allocatable :: etaOH

!  real, parameter :: gasct=8.314   ! gas constant
!  real, parameter :: mct=0.02894   ! kg mol-1
!  real, parameter :: g=9.80665     ! m s-2
!  real, parameter :: lrate=0.0065  ! K m-1
  real, parameter :: scalehgt=7000. ! scale height in metres

  ! Read OH fields and level heights
  !********************************

  do m=1,12
  
    ! open netcdf file
    write(mm,fmt='(i2.2)') m
!    thefile=trim(path(1))//'OH_FIELDS/'//'geos-chem.OH.2005'//mm//'01.nc'
    thefile=trim(ohfields_path)//'OH_FIELDS/'//'geos-chem.OH.2005'//mm//'01.nc'
    ierr=nf_open(trim(thefile),NF_NOWRITE,nid)
    if(ierr.ne.0) then
      write (*,*) 'The OH field at: '//thefile// ' could not be opened'
      write (*,*) 'please copy the OH fields there, or change the path in the'
      write (*,*) 'COMMAND namelist!'
      write(*,*) nf_strerror(ierr)
      stop
    endif

    ! inquire about variables
    ierr=nf_inq_dimid(nid,'Lon-000',xid)
    if(ierr.ne.0) then
      write(*,*) nf_strerror(ierr)
      stop
    endif
    ierr=nf_inq_dimid(nid,'Lat-000',yid)
    if(ierr.ne.0) then
      write(*,*) nf_strerror(ierr)
      stop
    endif
    ierr=nf_inq_dimid(nid,'Alt-000',zid)
    if(ierr.ne.0) then
      write(*,*) nf_strerror(ierr)
      stop
    endif

    if(m.eq.1) then

      ! read dimension sizes
      ierr=nf_inq_dimlen(nid,xid,nxOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif
      ierr=nf_inq_dimlen(nid,yid,nyOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif
      ierr=nf_inq_dimlen(nid,zid,nzOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif  

      ! allocate variables
      allocate(lonOH(nxOH))
      allocate(latOH(nyOH))
      allocate(etaOH(nzOH))
      allocate(altOH(nzOH))
      allocate(OH_field(nxOH,nyOH,nzOH,12))
      allocate(OH_hourly(nxOH,nyOH,nzOH,2))

      ! read dimension variables
      ierr=nf_inq_varid(nid,'LON',xid)
      ierr=nf_get_var_real(nid,xid,lonOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif
      ierr=nf_inq_varid(nid,'LAT',yid)
      ierr=nf_get_var_real(nid,yid,latOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif
      ierr=nf_inq_varid(nid,'ETAC',zid)
      ierr=nf_get_var_real(nid,zid,etaOH)
      if(ierr.ne.0) then
        write(*,*) nf_strerror(ierr)
        stop
      endif

      ! convert eta-level to altitude (assume surface pressure of 1010 hPa)
      altOH=log(1010./(etaOH*1010.))*scalehgt

    endif ! m.eq.1

    ! read OH_field
    ierr=nf_inq_varid(nid,'CHEM-L_S__OH',vid)
    ierr=nf_get_var_real(nid,vid,OH_field(:,:,:,m))
    if(ierr.ne.0) then
      write(*,*) nf_strerror(ierr)
      stop
    endif

    ierr=nf_close(nid)

  end do 
 
  deallocate(etaOH)

  ! Read J(O1D) photolysis rates
  !********************************  

  ! open netcdf file
!  thefile=trim(path(1))//'OH_FIELDS/jrate_average.nc'
  thefile=trim(ohfields_path)//'OH_FIELDS/jrate_average.nc'
  ierr=nf_open(trim(thefile),NF_NOWRITE,nid)
  if(ierr.ne.0) then
    write(*,*) nf_strerror(ierr)
    stop
  endif

  ! read dimension variables
  ierr=nf_inq_varid(nid,'longitude',xid)
  ierr=nf_get_var_real(nid,xid,lonjr)
  if(ierr.ne.0) then
    write(*,*) nf_strerror(ierr)
    stop
  endif
  ierr=nf_inq_varid(nid,'latitude',yid)
  ierr=nf_get_var_real(nid,yid,latjr)
  if(ierr.ne.0) then
    write(*,*) nf_strerror(ierr)
    stop
  endif

  ! read jrate_average
  ierr=nf_inq_varid(nid,'jrate',vid)
  ierr=nf_get_var_real(nid,vid,jrate_average)
  if(ierr.ne.0) then
    write(*,*) nf_strerror(ierr)
    stop
  endif

  ierr=nf_close(nid) 

  return

end subroutine readOHfield

