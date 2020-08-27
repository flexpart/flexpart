! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine getfields(itime,nstop,metdata_format)
!                       i     o
!*****************************************************************************
!                                                                            *
!  This subroutine manages the 3 data fields to be kept in memory.           *
!  During the first time step of petterssen it has to be fulfilled that the  *
!  first data field must have |wftime|<itime, i.e. the absolute value of     *
!  wftime must be smaller than the absolute value of the current time in [s].*
!  The other 2 fields are the next in time after the first one.              *
!  Pointers (memind) are used, because otherwise one would have to resort the*
!  wind fields, which costs a lot of computing time. Here only the pointers  *
!  are resorted.                                                             *
!                                                                            *
!     Author: A. Stohl                                                       *
!                                                                            *
!     29 April 1994                                                          *
!                                                                            *
!*****************************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:                                     *
!   Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.        *
!   Function of nstop extended.                                              *
!                                                                            *
!   Unified ECMWF and GFS builds                                             *
!   Marian Harustak, 12.5.2017                                               *
!     - Added passing of metdata_format as it was needed by called routines  *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! lwindinterval [s]    time difference between the two wind fields read in   *
! indj                 indicates the number of the wind field to be read in  *
! indmin               remembers the number of wind fields already treated   *
! memind(2)            pointer, on which place the wind fields are stored    *
! memtime(2) [s]       times of the wind fields, which are kept in memory    *
! itime [s]            current time since start date of trajectory calcu-    *
!                      lation                                                *
! nstop                > 0, if trajectory has to be terminated               *
! nx,ny,nuvz,nwz       field dimensions in x,y and z direction               *
! uu(0:nxmax,0:nymax,nuvzmax,2)   wind components in x-direction [m/s]       *
! vv(0:nxmax,0:nymax,nuvzmax,2)   wind components in y-direction [m/s]       *
! ww(0:nxmax,0:nymax,nwzmax,2)    wind components in z-direction [deltaeta/s]*
! tt(0:nxmax,0:nymax,nuvzmax,2)   temperature [K]                            *
! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                      *
! metdata_format     format of metdata (ecmwf/gfs)                           *
!                                                                            *
! Constants:                                                                 *
! idiffmax             maximum allowable time difference between 2 wind      *
!                      fields                                                *
!                                                                            *
!*****************************************************************************

  use par_mod
  use com_mod
  use class_gribfile

  implicit none

  integer :: indj,itime,nstop,memaux
  integer :: metdata_format

  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
! RLT added partial pressure water vapor 
  real :: pwater(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  integer :: kz, ix
  character(len=100) :: rowfmt

  integer :: indmin = 1


! Check, if wind fields are available for the current time step
!**************************************************************

  nstop=0

  if ((ldirect*wftime(1).gt.ldirect*itime).or. &
       (ldirect*wftime(numbwf).lt.ldirect*itime)) then
    write(*,*) 'FLEXPART WARNING: NO WIND FIELDS ARE AVAILABLE.'
    write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
    nstop=4
    return
  endif


  if ((ldirect*memtime(1).le.ldirect*itime).and. &
       (ldirect*memtime(2).gt.ldirect*itime)) then

! The right wind fields are already in memory -> don't do anything
!*****************************************************************

    continue

  else if ((ldirect*memtime(2).le.ldirect*itime).and. &
       (memtime(2).ne.999999999)) then


! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

    memaux=memind(1)
    memind(1)=memind(2)
    memind(2)=memaux
    memtime(1)=memtime(2)


! Read a new wind field and store it on place memind(2)
!******************************************************

    do indj=indmin,numbwf-1
      if (ldirect*wftime(indj+1).gt.ldirect*itime) then
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(2),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(2),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
        memtime(2)=wftime(indj+1)
        nstop = 1
        goto 40
      endif
    end do
40  indmin=indj

    if (WETBKDEP) then
      call writeprecip(itime,memind(1))
    endif

  else

! No wind fields, which can be used, are currently in memory
! -> read both wind fields
!***********************************************************

    do indj=indmin,numbwf-1
      if ((ldirect*wftime(indj).le.ldirect*itime).and. &
           (ldirect*wftime(indj+1).gt.ldirect*itime)) then
        memind(1)=1
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call readwind_ecmwf(indj,memind(1),uuh,vvh,wwh)
        else
          call readwind_gfs(indj,memind(1),uuh,vvh,wwh)
        end if
        call readwind_nests(indj,memind(1),uuhn,vvhn,wwhn)
        call calcpar(memind(1),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(1),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(1),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(1),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(1),uuhn,vvhn,wwhn,pvhn)
        memtime(1)=wftime(indj)
        memind(2)=2
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(2),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(2),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
        memtime(2)=wftime(indj+1)
        nstop = 1
        goto 60
      endif
    end do
60  indmin=indj

    if (WETBKDEP) then
      call writeprecip(itime,memind(1))
    endif

  end if

  ! RLT calculate dry air density
  pwater=qv*prs/((r_air/r_water)*(1.-qv)+qv)
  rho_dry=(prs-pwater)/(r_air*tt)

  ! test density
!  write(rowfmt,'(A,I6,A)') '(',nymax,'(E11.4,1X))'
!  if(itime.eq.0) then
!    open(500,file=path(2)(1:length(2))//'rho_dry.txt',status='replace',action='write')
!    do kz=1,nzmax
!      do ix=1,nxmax
!        write(500,fmt=rowfmt) rho_dry(ix,:,kz,1)
!      end do
!    end do
!    close(500)
!    open(500,file=path(2)(1:length(2))//'rho.txt',status='replace',action='write')
!    do kz=1,nzmax
!      do ix=1,nxmax
!        write(500,fmt=rowfmt) rho(ix,:,kz,1)
!      end do
!    end do
!    close(500)
!  endif

  lwindinterv=abs(memtime(2)-memtime(1))

  if (lwindinterv.gt.idiffmax) nstop=3

end subroutine getfields
