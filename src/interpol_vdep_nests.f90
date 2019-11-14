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

subroutine interpol_vdep_nests(level,vdepo)
  !                                 i     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of the deposition velocity on 2-d model layer.             *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     30 May 1994                                                           *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! level                number of species for which interpolation is done    *
  !                                                                           *
  !****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod

  implicit none

  integer :: level,indexh,m
  real :: y(2),vdepo

  ! a) Bilinear horizontal interpolation

  do m=1,2
    indexh=memind(m)

    y(m)=p1*vdepn(ix ,jy ,level,indexh,ngrid) &
         +p2*vdepn(ixp,jy ,level,indexh,ngrid) &
         +p3*vdepn(ix ,jyp,level,indexh,ngrid) &
         +p4*vdepn(ixp,jyp,level,indexh,ngrid)
  end do


  ! b) Temporal interpolation

  vdepo=(y(1)*dt2+y(2)*dt1)*dtt

  depoindicator(level)=.false.


end subroutine interpol_vdep_nests
