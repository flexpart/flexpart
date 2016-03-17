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

subroutine mean(x,xm,xs,number)

  !*****************************************************************************
  !                                                                            *
  !  This subroutine calculates mean and standard deviation of a given element.*
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 25 January 1994                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! x(number)           field of input data                                    *
  ! xm                  mean                                                   *
  ! xs                  standard deviation                                     *
  ! number              number of elements of field x                          *
  !                                                                            *
  ! Constants:                                                                 *
  ! eps                 tiny number                                            *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: number,i
  real :: x(number),xm,xs,xl,xq,xaux
  real,parameter :: eps=1.0e-30

  xl=0.
  xq=0.
  do i=1,number
    xl=xl+x(i)
    xq=xq+x(i)*x(i)
  end do

  xm=xl/real(number)

  xaux=xq-xl*xl/real(number)

  if (xaux.lt.eps) then
    xs=0.
  else
    xs=sqrt(xaux/real(number-1))
  endif

end subroutine mean
