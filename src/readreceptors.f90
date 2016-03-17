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

subroutine readreceptors

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the receptor points.    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 August 1996                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! receptorarea(maxreceptor)  area of dx*dy at location of receptor           *
  ! receptorname(maxreceptor)  names of receptors                              *
  ! xreceptor,yreceptor  coordinates of receptor points                        *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitreceptor         unit connected to file RECEPTORS                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: j
  real :: x,y,xm,ym
  character(len=16) :: receptor


  ! For backward runs, do not allow receptor output. Thus, set number of receptors to zero
  !*****************************************************************************

  if (ldirect.lt.0) then
    numreceptor=0
    return
  endif


  ! Open the RECEPTORS file and read output grid specifications
  !************************************************************

  open(unitreceptor,file=path(1)(1:length(1))//'RECEPTORS', &
       status='old',err=999)

  call skplin(5,unitreceptor)


  ! Read the names and coordinates of the receptors
  !************************************************

  j=0
100   j=j+1
    read(unitreceptor,*,end=99)
    read(unitreceptor,*,end=99)
    read(unitreceptor,*,end=99)
    read(unitreceptor,'(4x,a16)',end=99) receptor
    call skplin(3,unitreceptor)
    read(unitreceptor,'(4x,f11.4)',end=99) x
    call skplin(3,unitreceptor)
    read(unitreceptor,'(4x,f11.4)',end=99) y
    if ((x.eq.0.).and.(y.eq.0.).and. &
         (receptor.eq.'                ')) then
      j=j-1
      goto 100
    endif
    if (j.gt.maxreceptor) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY RECEPTOR #### '
    write(*,*) ' #### POINTS ARE GIVEN.                       #### '
    write(*,*) ' #### MAXIMUM NUMBER IS ',maxreceptor,'       #### '
    write(*,*) ' #### PLEASE MAKE CHANGES IN FILE RECEPTORS   #### '
    endif
    receptorname(j)=receptor
    xreceptor(j)=(x-xlon0)/dx       ! transform to grid coordinates
    yreceptor(j)=(y-ylat0)/dy
    xm=r_earth*cos(y*pi/180.)*dx/180.*pi
    ym=r_earth*dy/180.*pi
    receptorarea(j)=xm*ym
    goto 100

99   numreceptor=j-1
  close(unitreceptor)
  return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "RECEPTORS"  #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

end subroutine readreceptors
