! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine getrb(nc,ustar,nyl,diffh2o,reldiff,rb)
  !                 i    i    i     i       i    o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of the quasilaminar sublayer resistance to dry deposition.    *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 20 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! rb(ncmax)       sublayer resistance                                        *
  ! schmidt         Schmidt number                                             *
  ! ustar [m/s]     friction velocity                                          *
  ! diffh20 [m2/s]  diffusivity of water vapor in air                          *
  ! reldiff         diffusivity relative to H2O                                *
  !                                                                            *
  ! Constants:                                                                 *
  ! karman          von Karman constant                                        *
  ! pr              Prandtl number                                             *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  real :: ustar,diffh2o,rb(maxspec),schmidt,nyl
  real :: reldiff(maxspec)
  integer :: ic,nc
  real,parameter :: pr=0.72

  do ic=1,nc
    if (reldiff(ic).gt.0.) then
      schmidt=nyl/diffh2o*reldiff(ic)
      rb(ic)=2.0*(schmidt/pr)**0.67/(karman*ustar)
    endif
  end do

end subroutine getrb
