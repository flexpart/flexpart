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

subroutine shift_field(field,nxf,nyf,nzfmax,nzf,nmax,n)
  !                        i/o   i   i    i     i   i   i
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

  integer :: nxf,nyf,nzf,n,ix,jy,kz,ixs,nzfmax,nmax
  real :: field(0:nxmax-1,0:nymax-1,nzfmax,nmax),xshiftaux(0:nxmax-1)

  ! Loop over y and z
  !******************

  do kz=1,nzf
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
          xshiftaux(ixs)=field(ix,jy,kz,n)
        end do
        do ix=0,nxf-1
          field(ix,jy,kz,n)=xshiftaux(ix)
        end do
      endif

  ! Repeat the westernmost grid cells at the easternmost domain "boundary"
  !***********************************************************************

      field(nxf,jy,kz,n)=field(0,jy,kz,n)
    end do
  end do

end subroutine shift_field
