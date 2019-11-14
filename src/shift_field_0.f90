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

subroutine shift_field_0(field,nxf,nyf)
  !                          i/o   i   i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine shifts global fields by nxshift grid cells, in order to   *
  !  facilitate all sorts of nested wind fields, or output grids, which,       *
  !  without shifting, would overlap with the domain "boundary".               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    3 July 2002                                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  integer :: nxf,nyf,ix,jy,ixs
  real :: field(0:nxmax-1,0:nymax-1),xshiftaux(0:nxmax-1)

  ! Loop over y and z
  !******************

  do jy=0,nyf-1

  ! Shift the data
  !***************

    if (nxshift.ne.0) then
      do ix=0,nxf-1
        if (ix.ge.nxshift) then
          ixs=ix-nxshift
        else
          ixs=nxf-nxshift+ix
        endif
        xshiftaux(ixs)=field(ix,jy)
      end do
      do ix=0,nxf-1
        field(ix,jy)=xshiftaux(ix)
      end do
    endif

  ! Repeat the westernmost grid cells at the easternmost domain "boundary"
  !***********************************************************************

    field(nxf,jy)=field(0,jy)
  end do

  return
end subroutine shift_field_0
