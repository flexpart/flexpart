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

subroutine advance_rec(itime,xt,yt,zt,prob)
  !                    i     i  i  i  o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of turbulent particle trajectories utilizing a                *
  !  zero-acceleration scheme, which is corrected by a numerically more        *
  !  accurate Petterssen scheme whenever possible.                             *
  !                                                                            *
  !  Particle positions are read in, incremented, and returned to the calling  *
  !  program.                                                                  *
  !                                                                            *
  !  In different regions of the atmosphere (PBL vs. free troposphere),        *
  !  different parameters are needed for advection, parameterizing turbulent   *
  !  velocities, etc. For efficiency, different interpolation routines have    *
  !  been written for these different cases, with the disadvantage that there  *
  !  exist several routines doing almost the same. They all share the          *
  !  included file 'interpol_mod'. The following                               *
  !  interpolation routines are used:                                          *
  !                                                                            *
  !  interpol_all(_nests)     interpolates everything (called inside the PBL)  *
  !  interpol_misslev(_nests) if a particle moves vertically in the PBL,       *
  !                           additional parameters are interpolated if it     *
  !                           crosses a model level                            *
  !  interpol_wind(_nests)    interpolates the wind and determines the         *
  !                           standard deviation of the wind (called outside   *
  !                           PBL) also interpolates potential vorticity       *
  !  interpol_wind_short(_nests) only interpolates the wind (needed for the    *
  !                           Petterssen scheme)                               *
  !  interpol_vdep(_nests)    interpolates deposition velocities               *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Changes:                                                                  *
  !                                                                            *
  !  8 April 2000: Deep convection parameterization                            *
  !                                                                            *
  !  May 2002: Petterssen scheme introduced                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          time at which this subroutine is entered                *
  ! itimec [s]         actual time, which is incremented in this subroutine    *
  ! href [m]           height for which dry deposition velocity is calculated  *
  ! ldirect            1 forward, -1 backward                                  *
  ! ldt [s]            Time step for the next integration                      *
  ! lsynctime [s]      Synchronisation interval of FLEXPART                    *
  ! ngrid              index which grid is to be used                          *
  ! prob               probability of absorption due to dry deposition         *
  ! vdepo              Deposition velocities for all species                   *
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use interpol_mod

  implicit none

  real(kind=dp) :: xt,yt
  real :: zt,xtn,ytn
  integer :: itime,i,j,k,memindnext
  integer :: nix,njy,ks
  real :: prob(maxspec),vdepo(maxspec)
  real,parameter :: eps=nxmax/3.e5

  if (DRYDEP) then    ! reset probability for deposition
    do ks=1,nspec
      depoindicator(ks)=.true.
      prob(ks)=0.
    end do
  endif


  ! Determine whether lat/long grid or polarstereographic projection
  ! is to be used
  ! Furthermore, determine which nesting level to be used
  !*****************************************************************

  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngrid=j
        goto 23
      endif
    end do
23   continue
  endif


  !***************************
  ! Interpolate necessary data
  !***************************

  if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
    memindnext=1
  else
    memindnext=2
  endif

  ! Determine nested grid coordinates
  !**********************************

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)
  else
    ix=int(xt)
    jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
  endif
  ixp=ix+1
  jyp=jy+1


  ! Determine probability of deposition
  !************************************

      if ((DRYDEP).and.(zt.lt.2.*href)) then
        do ks=1,nspec
          if (DRYDEPSPEC(ks)) then
            if (depoindicator(ks)) then
              if (ngrid.le.0) then
                call interpol_vdep(ks,vdepo(ks))
              else
                call interpol_vdep_nests(ks,vdepo(ks))
              endif
            endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
               prob(ks)=vdepo(ks)
!               prob(ks)=vdepo(ks)/2./href
! instead of prob - return vdepo -> result kg/m2/s
          endif
        end do
      endif


end subroutine advance_rec

