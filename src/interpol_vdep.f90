! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine interpol_vdep(level,vdepo)
  !                           i     o
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
! write(*,*) 'interpol: ',dt1,dt2,dtt,lsynctime,ix,jy
  do m=1,2
    indexh=memind(m)

    y(m)=p1*vdep(ix ,jy ,level,indexh) &
         +p2*vdep(ixp,jy ,level,indexh) &
         +p3*vdep(ix ,jyp,level,indexh) &
         +p4*vdep(ixp,jyp,level,indexh)
  end do



  ! b) Temporal interpolation

  vdepo=(y(1)*dt2+y(2)*dt1)*dtt

  depoindicator(level)=.false.


end subroutine interpol_vdep
