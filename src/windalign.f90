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

subroutine windalign(u,v,ffap,ffcp,ux,vy)
  !                     i i  i    i   o  o
  !*****************************************************************************
  !                                                                            *
  !  Transformation from along- and cross-wind components to u and v           *
  !  components.                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     3 June 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ffap  turbulent wind in along wind direction                               *
  ! ffcp  turbulent wind in cross wind direction                               *
  ! u     main wind component in x direction                                   *
  ! ux    turbulent wind in x direction                                        *
  ! v     main wind component in y direction                                   *
  ! vy    turbulent wind in y direction                                        *
  !                                                                            *
  !*****************************************************************************

  implicit none

  real :: u,v,ffap,ffcp,ux,vy,ffinv,ux1,ux2,vy1,vy2,sinphi,cosphi
  real,parameter :: eps=1.e-30


  ! Transform along wind components
  !********************************

  ffinv=1./max(sqrt(u*u+v*v),eps)
  sinphi=v*ffinv
  vy1=sinphi*ffap
  cosphi=u*ffinv
  ux1=cosphi*ffap


  ! Transform cross wind components
  !********************************

  ux2=-sinphi*ffcp
  vy2=cosphi*ffcp


  ! Add contributions from along and cross wind components
  !*******************************************************

  ux=ux1+ux2
  vy=vy1+vy2

end subroutine windalign
