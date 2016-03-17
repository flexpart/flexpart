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

subroutine getrc(nc,i,j,t,gr,rh,rr,rc)
  !                 i  i i i i  i  i  o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of the surface resistance according to the procedure given    *
  !  in:                                                                       *
  !  Wesely (1989): Parameterization of surface resistances to gaseous         *
  !  dry deposition in regional-scale numerical models.                        *
  !  Atmos. Environ. 23, 1293-1304.                                            *
  !                                                                            *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 19 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! reldiff(maxspec)  diffusivity of H2O/diffusivity of component i            *
  ! gr [W/m2]       global radiation                                           *
  ! i               index of seasonal category                                 *
  ! j               index of landuse class                                     *
  ! ldep(maxspec)          1, if deposition shall be calculated for species i  *
  ! nc                   actual number of chemical components                  *
  ! rcl(maxspec,5,8) [s/m] Lower canopy resistance                             *
  ! rgs(maxspec,5,8) [s/m] Ground resistance                                   *
  ! rlu(maxspec,5,8) [s/m] Leaf cuticular resistance                           *
  ! rm(maxspec) [s/m]      Mesophyll resistance                                *
  ! t [C]           temperature                                                *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: i,j,ic,nc
  real :: gr,rh,rr,t,rs,rsm,corr,rluc,rclc,rgsc,rdc,rluo
  real :: rc(maxspec)


  ! Compute stomatal resistance
  !****************************
  ! Sabine Eckhardt, Dec 06: use 1E25 instead of 99999. for infinite res.

  if ((t.gt.0.).and.(t.lt.40.)) then
    rs=ri(i,j)*(1.+(200./(gr+0.1))**2)*(400./(t*(40.-t)))
  else
    rs=1.E25
  !  rs=99999.
  endif


  ! Correct stomatal resistance for effect of dew and rain
  !*******************************************************

  if ((rh.gt.0.9).or.(rr.gt.0.)) rs=rs*3.

  ! Compute the lower canopy resistance
  !************************************

  rdc=100.*(1.+1000./(gr+10.))


  corr=1000.*exp(-1.*t-4.)
  do ic=1,nc
    if (reldiff(ic).gt.0.) then

  ! Compute combined stomatal and mesophyll resistance
  !***************************************************

      rsm=rs*reldiff(ic)+rm(ic)

  ! Correct leaf cuticular, lower canopy and ground resistance
  !***********************************************************

      rluc=rlu(ic,i,j)+corr
      rclc=rcl(ic,i,j)+corr
      rgsc=rgs(ic,i,j)+corr

  ! Correct leaf cuticular resistance for effect of dew and rain
  !*************************************************************

      if (rr.gt.0.) then
        rluo=1./(1./1000.+1./(3.*rluc))
        rluc=1./(1./(3.*rluc)+1.e-7*henry(ic)+f0(ic)/rluo)
      else if (rh.gt.0.9) then
        rluo=1./(1./3000.+1./(3.*rluc))
        rluc=1./(1./(3.*rluc)+1.e-7*henry(ic)+f0(ic)/rluo)
      endif

  ! Combine resistances to give total resistance
  !*********************************************

      rc(ic)=1./(1./rsm+1./rluc+1./(rdc+rclc)+1./(rac(i,j)+rgsc))
  ! Sabine Eckhardt, Dec 06: avoid possible excessively high vdep
      if (rc(ic).lt.10.) rc(ic)=10.
    endif
  end do

end subroutine getrc
