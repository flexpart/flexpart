! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
  write(unitouttraj,'(i8,1x,i6,1x,a)') ibdate,ibtime, trim(flexversion)
  else
  write(unitouttraj,'(i8,1x,i6,1x,a)') iedate,ietime, trim(flexversion)
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
