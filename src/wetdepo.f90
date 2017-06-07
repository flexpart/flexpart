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

subroutine wetdepo(itime,ltsample,loutnext)
!                  i      i        i
!*****************************************************************************
!                                                                            *
! Calculation of wet deposition using the concept of scavenging coefficients.*
! For lack of detailed information, washout and rainout are jointly treated. *
! It is assumed that precipitation does not occur uniformly within the whole *
! grid cell, but that only a fraction of the grid cell experiences rainfall. *
! This fraction is parameterized from total cloud cover and rates of large   *
! scale and convective precipitation.                                        *
!                                                                            *
!    Author: A. Stohl                                                        *
!                                                                            *
!    1 December 1996                                                         *
!                                                                            *
! Correction by Petra Seibert, Sept 2002:                                    *
! use centred precipitation data for integration                             *
! Code may not be correct for decay of deposition!                           *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! ix,jy              indices of output grid cell for each particle           *
! itime [s]          actual simulation time [s]                              *
! jpart              particle index                                          *
! ldeltat [s]        interval since radioactive decay was computed           *
! loutnext [s]       time for which gridded deposition is next output        *
! loutstep [s]       interval at which gridded deposition is output          *
! ltsample [s]       interval over which mass is deposited                   *
! wetdeposit         mass that is wet deposited                              *
! wetgrid            accumulated deposited mass on output grid               *
! wetscav            scavenging coefficient                                  *
!                                                                            *
! Constants:                                                                 *
!                                                                            *
!*****************************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat
  integer :: itage,nage
  integer :: ks, kp
  integer(selected_int_kind(16)), dimension(nspec) :: blc_count, inc_count
  real :: grfraction(3),wetscav
  real :: wetdeposit(maxspec),restmass
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled

! Compute interval since radioactive decay of deposited mass was computed
!************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

! Loop over all particles
!************************

  blc_count(:)=0
  inc_count(:)=0

  do jpart=1,numpart

    if (itra1(jpart).eq.-999999999) goto 20
    if(ldirect.eq.1)then
      if (itra1(jpart).gt.itime) goto 20
    else
      if (itra1(jpart).lt.itime) goto 20
    endif

! Determine age class of the particle - nage is used for the kernel
!******************************************************************
     itage=abs(itra1(jpart)-itramem(jpart))
     do nage=1,nageclass
       if (itage.lt.lage(nage)) goto 33
     end do
 33  continue

    do ks=1,nspec      ! loop over species

      if (WETDEPSPEC(ks).eqv..false.) cycle 

!**************************************************
! CALCULATE DEPOSITION 
!**************************************************
!       wetscav=0.
       
!        write(*,*) ks,dquer(ks), crain_aero(ks),csnow_aero(ks)
!       if (((dquer(ks).le.0.).and.(weta_gas(ks).gt.0..or.wetb_gas(ks).gt.0.)) &
!          .or. &
!          ((dquer(ks).gt.0.).and.(crain_aero(ks).gt.0..or.csnow_aero(ks).gt.0.).or. &
!            (ccn_aero(ks).gt0) .or. (in_aero(ks).gt.0) .or. (henry(ks).gt.0)))  then

      call get_wetscav(itime,ltsample,loutnext,jpart,ks,grfraction,inc_count,blc_count,wetscav)
      

      if (wetscav.gt.0.) then
        wetdeposit(ks)=xmass1(jpart,ks)* &
             (1.-exp(-wetscav*abs(ltsample)))*grfraction(1)  ! wet deposition
      else ! if no scavenging
        wetdeposit(ks)=0.
      endif
 
      restmass = xmass1(jpart,ks)-wetdeposit(ks)
      if (ioutputforeachrelease.eq.1) then
        kp=npoint(jpart)
      else
        kp=1
      endif
      if (restmass .gt. smallnum) then
        xmass1(jpart,ks)=restmass
!   depostatistic
!   wetdepo_sum(ks,kp)=wetdepo_sum(ks,kp)+wetdeposit(ks)
!   depostatistic
      else
        xmass1(jpart,ks)=0.
      endif
!   Correct deposited mass to the last time step when radioactive decay of
!   gridded deposited mass was calculated
      if (decay(ks).gt.0.) then
        wetdeposit(ks)=wetdeposit(ks)*exp(abs(ldeltat)*decay(ks))
      endif

!    endif ! no deposition
    end do ! loop over species

! Sabine Eckhardt, June 2008 create deposition runs only for forward runs
! Add the wet deposition to accumulated amount on output grid and nested output grid
!*****************************************************************************

    if (ldirect.eq.1) then
      call wetdepokernel(nclass(jpart),wetdeposit,real(xtra1(jpart)), &
           real(ytra1(jpart)),nage,kp)
      if (nested_output.eq.1) call wetdepokernel_nest(nclass(jpart), &
           wetdeposit,real(xtra1(jpart)),real(ytra1(jpart)),nage,kp)
    endif

20  continue
  end do ! all particles

! count the total number of below-cloud and in-cloud occurences:
  tot_blc_count(1:nspec)=tot_blc_count(1:nspec)+blc_count(1:nspec)
  tot_inc_count(1:nspec)=tot_inc_count(1:nspec)+inc_count(1:nspec)

end subroutine wetdepo
