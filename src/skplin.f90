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

subroutine skplin(nlines,iunit)
  !                    i      i
  !*****************************************************************************
  !                                                                            *
  !   This routine reads nlines from unit iunit and discards them              *
  !                                                                            *
  !   Authors: Petra Seibert                                                   *
  !                                                                            *
  !   31 Dec 1998                                                              *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! iunit   unit number from which lines are to be skipped                     *
  ! nlines  number of lines to be skipped                                      *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,iunit, nlines

  do i=1,nlines
    read(iunit,*)
  end do

end subroutine skplin
