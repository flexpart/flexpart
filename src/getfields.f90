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

subroutine getfields(itime,nstop,id_centre)
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
!     - Added passing of id_centre as it was needed by called routines       *
!                                                                            *
!  Petra Seibert, 2018-06-26: simplified version met data format detection   *
!                                                                            *
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
! id_centre            format of metdata (ecmwf/gfs)                         *
!                                                                            *
! Constants:                                                                 *
! idiffmax             maximum allowable time difference between 2 wind      *
!                      fields                                                *
!                                                                            *
!*****************************************************************************

  use par_mod
  use com_mod
  use check_gribfile_mod

  implicit none

  integer :: indj,itime,nstop,memaux
  integer :: id_centre

  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)

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
        if (id_centre.eq.icg_id_ecmwf) then
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,id_centre)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,id_centre)
        if (id_centre.eq.icg_id_ecmwf) then
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
        if (id_centre.eq.icg_id_ecmwf) then
          call readwind_ecmwf(indj,memind(1),uuh,vvh,wwh)
        else
          call readwind_gfs(indj,memind(1),uuh,vvh,wwh)
        end if
        call readwind_nests(indj,memind(1),uuhn,vvhn,wwhn)
        call calcpar(memind(1),uuh,vvh,pvh,id_centre)
        call calcpar_nests(memind(1),uuhn,vvhn,pvhn,id_centre)
        if (id_centre.eq.icg_id_ecmwf) then
          call verttransform_ecmwf(memind(1),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(1),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(1),uuhn,vvhn,wwhn,pvhn)
        memtime(1)=wftime(indj)
        memind(2)=2
        if (id_centre.eq.icg_id_ecmwf) then
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,id_centre)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,id_centre)
        if (id_centre.eq.icg_id_ecmwf) then
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

  endif

  lwindinterv=abs(memtime(2)-memtime(1))

  if (lwindinterv.gt.idiffmax) nstop=3

end subroutine getfields
