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
  ! AUTHOR: Sabine Eckhardt, June 2007                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! i                       loop indices                                       *
  ! LENGTH(numpath)         length of the path names                           *
  ! PATH(numpath)           contains the path names                            *
  ! unitoh                  unit connected with OH field                       *
  !                                                                            *
  ! -----                                                                      *
  !                                                                            *
  !*****************************************************************************

  use oh_mod
  use par_mod
  use com_mod

  implicit none

  integer :: ix,jy,lev,m


  ! Read OH field and level heights
  !********************************

! write (*,*) 'reading OH'
  open(unitOH,file=path(1)(1:length(1))//'OH_7lev_agl.dat', &
       status='old',form='UNFORMATTED', err=998)
  do m=1,12
    do lev=1,maxzOH
      do ix=0,maxxOH-1
  !      do 10 jy=0,maxyOH-1
          read(unitOH) (OH_field(m,ix,jy,lev),jy=0,maxyOH-1)
  !      if ((ix.eq.20).and.(lev.eq.1)) then
  !          write(*,*) 'reading: ', m, OH_field(m,ix,20,lev)
  !      endif
      end do
    end do
  end do
  close(unitOH)

  do lev=1,7
    OH_field_height(lev)=1000+real(lev-1)*2.*1000.
  end do

!  write (*,*) 'OH read'
  return

  ! Issue error messages
  !*********************

998   write(*,*) ' #### FLEXPART ERROR! FILE CONTAINING          ####'
  write(*,*) ' #### OH FIELD DOES NOT EXIST                  ####'
  stop

end subroutine readohfield
