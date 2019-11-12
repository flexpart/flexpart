subroutine getfields(itime,nstop,memstat,metdata_format)
!                       i     o       o
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
!  Changes, Bernd C. Krueger, Feb. 2001:
!   Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.
!   Function of nstop extended.
!
!  eso 2014:
!    MPI version. 
!    If running with number of processes >= mpi_mod::read_grp_min,
!    only one process (mpi_mod::lmpreader=.true.) does the actual reading, while
!    the rest call this routine only to update memind, memstat etc.
!
!    If mpi_mod::lmp_sync=.true., uses 3 fields instead of 2, to allow reading 
!    the newest in the background.
!
!    Return memstat, which is the sum of 
!    
!    memstat=16:     one new field to be read
!    memstat=32:     two new fields must be read
!    memstat=1,2,3:  position(s) in array to read next field
!    memstat=0:      no new fields to be read
!          
!   Unified ECMWF and GFS builds                                             
!   Marian Harustak, 12.5.2017                                               
!     - Added passing of metdata_format as it was needed by called routines  
!          
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! lwindinterval [s]    time difference between the two wind fields read in   *
! indj                 indicates the number of the wind field to be read in  *
! indmin               remembers the number of wind fields already treated   *
! memind(2[3])         pointer, on which place the wind fields are stored    *
! memtime(2[3]) [s]    times of the wind fields, which are kept in memory    *
! itime [s]            current time since start date of trajectory calcu-    *
!                      lation                                                *
! nstop                > 0, if trajectory has to be terminated               *
! memstat              keep track of change of status for field pointers     *
! nx,ny,nuvz,nwz       field dimensions in x,y and z direction               *
! uu(0:nxmax,0:nymax,nuvzmax,numwfmem)   wind components in x-direction [m/s]*
! vv(0:nxmax,0:nymax,nuvzmax,numwfmem)   wind components in y-direction [m/s]*
! ww(0:nxmax,0:nymax,nwzmax,numwfmem)    wind components in z-direction      *
!                                          [deltaeta/s]                      *
! tt(0:nxmax,0:nymax,nuvzmax,numwfmem)   temperature [K]                     *
! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
! metdata_format     format of metdata (ecmwf/gfs)                           *
!                                                                            *
! Constants:                                                                 *
! idiffmax             maximum allowable time difference between 2 wind      *
!                      fields                                                *
!                                                                            *
!*****************************************************************************

  use par_mod
  use com_mod
  use mpi_mod, only: lmpreader,lmp_use_reader,lmp_sync
  use class_gribfile

  implicit none

  integer :: metdata_format
  integer :: indj,itime,nstop,memaux,mindread
! mindread: index of where to read 3rd field
  integer, intent(out) :: memstat
! keep track of 3rd field index. updated when new fields are read
  integer, save :: mind3=0

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
  memstat=0

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
    memstat=0
    continue

  else if ((ldirect*memtime(2).le.ldirect*itime).and. &
       (memtime(2).ne.999999999)) then


! Current time is after 2nd wind field
! -> Resort wind field pointers, so that current time is between 1st and 2nd
!***************************************************************************

! 2 fields in memory
!*******************
    if (lmp_sync) then
      memaux=memind(1)
      memind(1)=memind(2)
      memind(2)=memaux
      memtime(1)=memtime(2)
      memstat=16+memind(2)
      mindread=memind(2)

    else 
! 3 fields in memory
! Note that the reader process never needs to keep 3 fields in memory,
! (2 is enough) so it is possible to save some memory
!*********************************************************************
      if (mind3.eq.32.or.mind3.eq.19) then
        if (lmpreader) then 
          memind(1)=2
          memind(2)=3
          memind(3)=3
        else
          memind(1)=2
          memind(2)=3
          memind(3)=1
        end if
        memstat=17
      else if (mind3.eq.17) then
        if (lmpreader) then
          memind(1)=3
          memind(2)=1
          memind(3)=1
        else
          memind(1)=3
          memind(2)=1
          memind(3)=2
        end if
        memstat=18
      else if (mind3.eq.18) then
        if (lmpreader) then
          memind(1)=1
          memind(2)=2
          memind(3)=2
        else
          memind(1)=1
          memind(2)=2
          memind(3)=3
        end if
        memstat=19
      else
        write(*,*) '#### getfields_mpi.f90> ERROR: ', &
             & 'unknown mind3, exiting.', mind3,' ####'
        stop
      end if
      mind3=memstat
      memtime(1)=memtime(2)
      mindread=memind(3)
    end if


! Read a new wind field and store it on place memind(2)
! or memind(3) for the 3-field read-ahead version
!******************************************************
    do indj=indmin,numbwf-1
      if (ldirect*wftime(indj+1).gt.ldirect*itime) then
        if (lmpreader.or..not. lmp_use_reader) then
          if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
            call readwind_ecmwf(indj+1,mindread,uuh,vvh,wwh)
          else
            call readwind_gfs(indj+1,mindread,uuh,vvh,wwh)
          end if
          call readwind_nests(indj+1,mindread,uuhn,vvhn,wwhn)
          call calcpar(mindread,uuh,vvh,pvh,metdata_format)
          call calcpar_nests(mindread,uuhn,vvhn,pvhn,metdata_format)
          if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
            call verttransform_ecmwf(mindread,uuh,vvh,wwh,pvh)
          else
            call verttransform_gfs(mindread,uuh,vvh,wwh,pvh)
          end if
          call verttransform_nests(mindread,uuhn,vvhn,wwhn,pvhn)
        end if
        memtime(2)=wftime(indj+1)
        nstop = 1
        goto 40
      endif
    end do
40  indmin=indj

   if (WETBKDEP.and.(lmpreader.or.(.not.lmp_use_reader.and.lroot))) then
        call writeprecip(itime,memind(1))
      endif

  else

! No wind fields, which can be used, are currently in memory
! -> read both wind fields
!***********************************************************
    memstat=32

    do indj=indmin,numbwf-1
      if ((ldirect*wftime(indj).le.ldirect*itime).and. &
           (ldirect*wftime(indj+1).gt.ldirect*itime)) then
        memind(1)=1
        if (lmpreader.or..not.lmp_use_reader) then
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
        end if
        memtime(1)=wftime(indj)
        memind(2)=2
        if (lmpreader.or..not.lmp_use_reader) then
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
        end if
        memtime(2)=wftime(indj+1)
! DEV: not used
        if (.not.lmp_sync) memind(3)=3 ! indicate position of next read
        nstop = 1
        goto 60
      endif
    end do
60  indmin=indj

!   if (WETBKDEP.and.lroot) then
    if (WETBKDEP.and.(lmpreader.or.(.not.lmp_use_reader.and.lroot))) then
      call writeprecip(itime,memind(1))
    endif

  endif

  lwindinterv=abs(memtime(2)-memtime(1))

  if (lwindinterv.gt.idiffmax) nstop=3

end subroutine getfields
