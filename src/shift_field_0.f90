! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
