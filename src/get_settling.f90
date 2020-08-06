! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine get_settling(itime,xt,yt,zt,nsp,settling)
  !                          i   i  i  i   i     o
  !*****************************************************************************
  !                                                                            *
  !  This subroutine calculates particle settling velocity.                    *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    May 2010                                                                *
  !                                                                            *
  !  Improvement over traditional settling calculation in FLEXPART:            *
  !  generalize to higher Reynolds numbers and also take into account the      *
  !  temperature dependence of dynamic viscosity.                              *
  !                                                                            *
  !  Based on:                                                                 *
  !  Naeslund E., and Thaning, L. (1991): On the settling velocity in a        *
  !  nonstationary atmosphere, Aerosol Science and Technology 14, 247-256.     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zt           coordinates position for which wind data shall be cal-  *
  !                    culated                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: itime,indz
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: rho1(2),tt1(2),temperature,airdens,vis_dyn,vis_kin,viscosity
  real :: settling,settling_old,reynolds,c_d
  integer :: i,n,nix,njy,indzh,nsp


  !*****************************************************************************
  ! 1. Interpolate temperature and density: nearest neighbor interpolation sufficient
  !*****************************************************************************

  nix=int(xt)
  njy=int(yt)

  ! Determine the level below the current position for u,v
  !*******************************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      goto 6
    endif
  end do
6   continue


  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz


  ! Bilinear horizontal interpolation
  !**********************************

  ! Loop over 2 levels
  !*******************

  do n=1,2
    indzh=indz+n-1
    rho1(n)=rho(nix,njy,indzh,1)
    tt1(n)=tt(nix,njy,indzh,1)
  end do


  ! Linear vertical interpolation
  !******************************

  temperature=dz2*tt1(1)+dz1*tt1(2)
  airdens=dz2*rho1(1)+dz1*rho1(2)


  vis_dyn=viscosity(temperature)
  vis_kin=vis_dyn/airdens

  reynolds=dquer(nsp)/1.e6*abs(vsetaver(nsp))/vis_kin



  ! Iteration to determine both Reynolds number and settling velocity
  !******************************************************************

  settling_old=vsetaver(nsp)    ! initialize iteration with Stokes' law, constant viscosity estimate

  do i=1,20    ! do a few iterations

    if (reynolds.lt.1.917) then
      c_d=24./reynolds
    else if (reynolds.lt.500.) then
      c_d=18.5/(reynolds**0.6)
    else
      c_d=0.44
    endif

    settling=-1.* &
         sqrt(4*ga*dquer(nsp)/1.e6*density(nsp)*cunningham(nsp)/ &
         (3.*c_d*airdens))

    if (abs((settling-settling_old)/settling).lt.0.01) goto 11  ! stop iteration

    reynolds=dquer(nsp)/1.e6*abs(settling)/vis_kin
    settling_old=settling
  end do

11   continue

end subroutine get_settling
