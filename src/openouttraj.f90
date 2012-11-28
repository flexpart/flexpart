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

subroutine openouttraj

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the output file for the plume trajectory output       *
  !   produced by the cluster analysis.                                        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     27 January 2001                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i
  real :: xp1,yp1,xp2,yp2


  ! Open output file for trajectory output
  !***************************************

  open(unitouttraj,file=path(2)(1:length(2))//'trajectories.txt', &
       form='formatted',err=998)

  if (ldirect.eq.1) then
  write(unitouttraj,'(i8,1x,i6,1x,a)') ibdate,ibtime,'FLEXPART V8.2'
  else
  write(unitouttraj,'(i8,1x,i6,1x,a)') iedate,ietime,'FLEXPART V8.2'
  endif
  write(unitouttraj,*) method,lsubgrid,lconvection
  write(unitouttraj,*) numpoint
  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitouttraj,*) ireleasestart(i),ireleaseend(i), &
         xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i),kindz(i),npart(i)
    if (numpoint.le.1000) then
      write(unitouttraj,'(a)') compoint(i)(1:40)
    else
      write(unitouttraj,'(a)') compoint(1001)(1:40)
    endif
  end do

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### trajectories.txt                         #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine openouttraj
