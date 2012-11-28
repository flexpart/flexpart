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
