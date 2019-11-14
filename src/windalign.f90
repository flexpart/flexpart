! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
