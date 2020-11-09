! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
