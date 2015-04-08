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

subroutine readavailable

  !*****************************************************************************
  !                                                                            *
  !   This routine reads the dates and times for which windfields are          *
  !   available.                                                               *
  !                                                                            *
  !     Authors: A. Stohl                                                      *
  !                                                                            *
  !     6 February 1994                                                        *
  !     8 February 1999, Use of nested fields, A. Stohl                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! bdate                beginning date as Julian date                         *
  ! beg                  beginning date for windfields                         *
  ! end                  ending date for windfields                            *
  ! fname                filename of wind field, help variable                 *
  ! ideltas [s]          duration of modelling period                          *
  ! idiff                time difference between 2 wind fields                 *
  ! idiffnorm            normal time difference between 2 wind fields          *
  ! idiffmax [s]         maximum allowable time between 2 wind fields          *
  ! jul                  julian date, help variable                            *
  ! numbwf               actual number of wind fields                          *
  ! wfname(maxwf)        file names of needed wind fields                      *
  ! wfspec(maxwf)        file specifications of wind fields (e.g., if on disc) *
  ! wftime(maxwf) [s]times of wind fields relative to beginning time           *
  ! wfname1,wfspec1,wftime1 = same as above, but only local (help variables)   *
  !                                                                            *
  ! Constants:                                                                 *
  ! maxwf                maximum number of wind fields                         *
  ! unitavailab          unit connected to file AVAILABLE                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: i,idiff,ldat,ltim,wftime1(maxwf),numbwfn(maxnests),k
  integer :: wftime1n(maxnests,maxwf),wftimen(maxnests,maxwf)
  real(kind=dp) :: juldate,jul,beg,end
  character(len=255) :: fname,spec,wfname1(maxwf),wfspec1(maxwf)
  character(len=255) :: wfname1n(maxnests,maxwf)
  character(len=40) :: wfspec1n(maxnests,maxwf)


  ! Windfields are only used, if they are within the modelling period.
  ! However, 1 additional day at the beginning and at the end is used for
  ! interpolation. -> Compute beginning and ending date for the windfields.
  !************************************************************************

  if (ideltas.gt.0) then         ! forward trajectories
    beg=bdate-1._dp
    end=bdate+real(ideltas,kind=dp)/86400._dp+real(idiffmax,kind=dp)/ &
         86400._dp
  else                           ! backward trajectories
    beg=bdate+real(ideltas,kind=dp)/86400._dp-real(idiffmax,kind=dp)/ &
         86400._dp
    end=bdate+1._dp
  endif

  ! Open the wind field availability file and read available wind fields
  ! within the modelling period.
  !*********************************************************************

  open(unitavailab,file=path(4)(1:length(4)),status='old', &
       err=999)

  do i=1,3
    read(unitavailab,*)
  end do

  numbwf=0
100   read(unitavailab,'(i8,1x,i6,2(6x,a255))',end=99) &
           ldat,ltim,fname,spec
    jul=juldate(ldat,ltim)
    if ((jul.ge.beg).and.(jul.le.end)) then
      numbwf=numbwf+1
      if (numbwf.gt.maxwf) then      ! check exceedance of dimension
       write(*,*) 'Number of wind fields needed is too great.'
       write(*,*) 'Reduce modelling period (file "COMMAND") or'
       write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
       stop
      endif

      wfname1(numbwf)=fname(1:index(fname,' '))
      wfspec1(numbwf)=spec
      wftime1(numbwf)=nint((jul-bdate)*86400._dp)
    endif
    goto 100       ! next wind field

99   continue

  close(unitavailab)

  ! Open the wind field availability file and read available wind fields
  ! within the modelling period (nested grids)
  !*********************************************************************

  do k=1,numbnests
  !print*,length(numpath+2*(k-1)+1),length(numpath+2*(k-1)+2),length(4),length(3)
  !print*,path(numpath+2*(k-1)+2)(1:length(numpath+2*(k-1)+2))
    open(unitavailab,file=path(numpath+2*(k-1)+2) &
         (1:length(numpath+2*(k-1)+2)),status='old',err=998)

    do i=1,3
      read(unitavailab,*)
    end do

    numbwfn(k)=0
700   read(unitavailab,'(i8,1x,i6,2(6x,a255))',end=699) ldat, &
           ltim,fname,spec
      jul=juldate(ldat,ltim)
      if ((jul.ge.beg).and.(jul.le.end)) then
        numbwfn(k)=numbwfn(k)+1
        if (numbwfn(k).gt.maxwf) then      ! check exceedance of dimension
       write(*,*) 'Number of nested wind fields is too great.'
       write(*,*) 'Reduce modelling period (file "COMMAND") or'
       write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
          stop
        endif

        wfname1n(k,numbwfn(k))=fname
        wfspec1n(k,numbwfn(k))=spec
        wftime1n(k,numbwfn(k))=nint((jul-bdate)*86400._dp)
      endif
      goto 700       ! next wind field

699   continue

    close(unitavailab)
  end do


  ! Check wind field times of file AVAILABLE (expected to be in temporal order)
  !****************************************************************************

  if (numbwf.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! NO WIND FIELDS    #### '
    write(*,*) ' #### AVAILABLE FOR SELECTED TIME PERIOD.     #### '
    stop
  endif

  do i=2,numbwf
    if (wftime1(i).le.wftime1(i-1)) then
      write(*,*) 'FLEXPART ERROR: FILE AVAILABLE IS CORRUPT.'
      write(*,*) 'THE WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1(i)
      stop
    endif
  end do

  ! Check wind field times of file AVAILABLE for the nested fields
  ! (expected to be in temporal order)
  !***************************************************************

  do k=1,numbnests
    if (numbwfn(k).eq.0) then
      write(*,*) '#### FLEXPART MODEL ERROR! NO WIND FIELDS  ####'
      write(*,*) '#### AVAILABLE FOR SELECTED TIME PERIOD.   ####'
      stop
    endif

    do i=2,numbwfn(k)
      if (wftime1n(k,i).le.wftime1n(k,i-1)) then
      write(*,*) 'FLEXPART ERROR: FILE AVAILABLE IS CORRUPT. '
      write(*,*) 'THE NESTED WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1n(k,i)
      write(*,*) 'AT NESTING LEVEL ',k
      stop
      endif
    end do

  end do


  ! For backward trajectories, reverse the order of the windfields
  !***************************************************************

  if (ideltas.ge.0) then
    do i=1,numbwf
      wfname(i)=wfname1(i)
      wfspec(i)=wfspec1(i)
      wftime(i)=wftime1(i)
    end do
    do k=1,numbnests
      do i=1,numbwfn(k)
        wfnamen(k,i)=wfname1n(k,i)
        wfspecn(k,i)=wfspec1n(k,i)
        wftimen(k,i)=wftime1n(k,i)
      end do
    end do
  else
    do i=1,numbwf
      wfname(numbwf-i+1)=wfname1(i)
      wfspec(numbwf-i+1)=wfspec1(i)
      wftime(numbwf-i+1)=wftime1(i)
    end do
    do k=1,numbnests
      do i=1,numbwfn(k)
        wfnamen(k,numbwfn(k)-i+1)=wfname1n(k,i)
        wfspecn(k,numbwfn(k)-i+1)=wfspec1n(k,i)
        wftimen(k,numbwfn(k)-i+1)=wftime1n(k,i)
      end do
    end do
  endif

  ! Check the time difference between the wind fields. If it is big,
  ! write a warning message. If it is too big, terminate the trajectory.
  !*********************************************************************

  do i=2,numbwf
    idiff=abs(wftime(i)-wftime(i-1))
    if (idiff.gt.idiffmax) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS TOO BIG FOR TRANSPORT CALCULATION.&
           &'
      write(*,*) 'THEREFORE, TRAJECTORIES HAVE TO BE SKIPPED.'
    else if (idiff.gt.idiffnorm) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS BIG. THIS MAY CAUSE A DEGRADATION'
      write(*,*) 'OF SIMULATION QUALITY.'
    endif
  end do

  do k=1,numbnests
    if (numbwfn(k).ne.numbwf) then
      write(*,*) 'FLEXPART ERROR: THE AVAILABLE FILES FOR THE'
      write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
      write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
      write(*,*) 'ERROR AT NEST LEVEL: ',k
      stop
    endif
    do i=1,numbwf
      if (wftimen(k,i).ne.wftime(i)) then
        write(*,*) 'FLEXPART ERROR: THE AVAILABLE FILES FOR THE'
        write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
        write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
        write(*,*) 'ERROR AT NEST LEVEL: ',k
        stop
      endif
    end do
  end do

  ! Reset the times of the wind fields that are kept in memory to no time
  !**********************************************************************

  do i=1,2
    memind(i)=i
    memtime(i)=999999999
  end do

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR! AVAILABLE FILE   #### '
  write(*,'(a)') '     '//path(numpath+2*(k-1)+2) &
       (1:length(numpath+2*(k-1)+2))
  write(*,*) ' #### CANNOT BE OPENED             #### '
  stop

999   write(*,*) ' #### FLEXPART MODEL ERROR! AVAILABLE IILE #### '
  write(*,'(a)') '     '//path(4)(1:length(4))
  write(*,*) ' #### CANNOT BE OPENED           #### '
  stop

end subroutine readavailable
