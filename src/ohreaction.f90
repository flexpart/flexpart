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

subroutine ohreaction(itime,ltsample,loutnext)
  !                     i      i        i
  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  !    Author: R.L. Thompson                                                   *
  !                                                                            *
  !    Nov 2014                                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Variables:                                                                 *
  ! ix,jy                indices of output grid cell for each particle         *
  ! itime [s]            actual simulation time [s]                            *
  ! jpart                particle index                                        *
  ! ldeltat [s]          interval since radioactive decay was computed         *
  ! loutnext [s]         time for which gridded deposition is next output      *
  ! loutstep [s]         interval at which gridded deposition is output        *
  ! oh_average [molecule/cm^3] OH Concentration                                *
  ! ltsample [s]         interval over which mass is deposited                 *
  !                                                                            *
  !*****************************************************************************

  use oh_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat,j,k,ix,jy!,ijx,jjy
  integer :: ngrid,interp_time,n,m,h,indz,i!,ia,il
  integer :: jjjjmmdd,hhmmss,OHx,OHy,OHz,NH3x,Nh3y,NH3z
  real, dimension(nzOH) :: alttopOH
  real, dimension(nzNH3) :: alttopNH3
  real :: xlon,ylat 
  real :: xtn,ytn
  real :: restmass,ohreacted,oh_average,nh3reachted,nh3_average
  real :: ohrate,nh3lifetime,temp,nh3reacted
  real, parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real(kind=dp) :: jul

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

  jul=bdate+real(itime,kind=dp)/86400.
  call caldate(jul,jjjjmmdd,hhmmss)
  m=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
  h=hhmmss/10000

  ! Loop over particles
  !*****************************************

  do jpart=1,numpart

    ! Determine which nesting level to be used
    ngrid=0
    do j=numbnests,1,-1
      if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and. &
           (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
        ngrid=j
        goto 23
      endif
    end do
23  continue

    ! Determine nested grid coordinates
    if (ngrid.gt.0) then
      xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
      ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
    else
      ix=int(xtra1(jpart))
      jy=int(ytra1(jpart))
    endif

    interp_time=nint(itime-0.5*ltsample)
    n=2
    if(abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) n=1

    do i=2,nz
      if (height(i).gt.ztra1(jpart)) then
        indz=i-1
        goto 6
      endif
    end do
6   continue

    ! Get OH from nearest grid-cell and specific month 
    !*************************************************

    ! world coordinates
    xlon=xtra1(jpart)*dx+xlon0
    if (xlon.gt.180) then
       xlon=xlon-360
    endif
    ylat=ytra1(jpart)*dy+ylat0

    ! do first oh reaction then the NH3 loss
    if (OHREA) then
      ! get position in the OH field
      OHx=minloc(abs(lonOH-xlon),dim=1,mask=abs(lonOH-xlon).eq.minval(abs(lonOH-xlon)))
      OHy=minloc(abs(latOH-ylat),dim=1,mask=abs(latOH-ylat).eq.minval(abs(latOH-ylat)))

      ! get the level of the OH field for the particle
      ! ztra1 is the z-coord of the trajectory above model orography in metres
      ! altOH is the height of the centre of the level in the OH field above orography
      do i=2,nzOH
        alttopOH(i-1)=altOH(i)+0.5*(altOH(i)-altOH(i-1))
      end do
      alttopOH(nzOH)=altOH(nzOH)+0.5*(altOH(nzOH)-altOH(nzOH-1))

      OHz=minloc(abs(alttopOH-ztra1(jpart)),dim=1,mask=abs(alttopOH-ztra1(jpart))&
            &.eq.minval(abs(alttopOH-ztra1(jpart))))

      ! Interpolate between hourly OH fields to current time
      !*****************************************************

      oh_average=OH_hourly(OHx,OHy,OHz,1)+&
               &(OH_hourly(OHx,OHy,OHz,2)-OH_hourly(OHx,OHy,OHz,1))*&
               &(itime-memOHtime(1))/(memOHtime(2)-memOHtime(1))

      if (oh_average.gt.smallnum) then

      ! Computation of the OH reaction
      !**********************************************************

        temp=tt(ix,jy,indz,n)

        do k=1,nspec                                
          if (ohcconst(k).gt.0.) then
            ohrate=ohcconst(k)*temp**ohnconst(k)*exp(-ohdconst(k)/temp)*oh_average
            ! new particle mass
            restmass = xmass1(jpart,k)*exp(-1*ohrate*abs(ltsample))
            if (restmass .gt. smallnum) then
              xmass1(jpart,k)=restmass
            else
              xmass1(jpart,k)=0.
            endif
            ohreacted=xmass1(jpart,k)*(1-exp(-1*ohrate*abs(ltsample)))
          else
            ohreacted=0.
        endif
      end do

      endif  ! oh_average.gt.smallnum 
    endif ! OHREA

    if (NH3LOSS) then
      ! get position in the NH3 field
      NH3x=minloc(abs(lonNH3-xlon),dim=1,mask=abs(lonNH3-xlon).eq.minval(abs(lonNH3-xlon)))
      NH3y=minloc(abs(latNH3-ylat),dim=1,mask=abs(latNH3-ylat).eq.minval(abs(latNH3-ylat)))

      ! get the level of the NH3 field for the particle
      ! ztra1 is the z-coord of the trajectory above model orography in metres
      ! altNH3 is the height of the centre of the level in the NH3 field above orography
      do i=2,nzNH3
        alttopNH3(i-1)=altNH3(i-1)+0.5*(altNH3(i)-altNH3(i-1))
      end do

      alttopNH3(nzNH3)=altNH3(nzNH3)+0.5*(altNH3(nzNH3)-altNH3(nzNH3-1))

      NH3z=minloc(abs(alttopNH3-ztra1(jpart)),dim=1,mask=abs(alttopNH3-ztra1(jpart))&
            &.eq.minval(abs(alttopNH3-ztra1(jpart))))


      ! There is always the OH field closest to the actual model time step read in
      !***************************************************************************

 !      write(*,*) "NH3LOSS_field(NH3x,NH3y,NH3z)", NH3LOSS_field(92,63,10)
      
      temp = 1.0-1800.0*NH3LOSS_field(NH3x,NH3y,NH3z)
      ! Avoid taking log of value close to 1.0
      if (abs(temp-1.0) < 10000.0*epsilon(0.0)) then
        nh3lifetime=sign(999999.0, NH3LOSS_field(NH3x,NH3y,NH3z))
      else
        nh3lifetime=-1800./log(1.0-1800.0*NH3LOSS_field(NH3x,NH3y,NH3z))
      end if

      if (nh3lifetime > 999999.0) nh3lifetime = 999999.0
      if (nh3lifetime < -999999.0) nh3lifetime = -999999.0


      ! if (1.0-1800.0*NH3LOSS_field(NH3x,NH3y,NH3z) <= 0.0) then
      !   write(*,*) "WARNING: 1-1800.0*NH3LOSS_field(NH3x,NH3y,NH3z) =< 0.0 "
      ! end if
      
      

!     if ((nh3lifetime).gt.smallnum) then  !only for loss

      ! Computation of the NH3 loss reaction
      !*************************************

        do k=1,nspec                                
         if (use_NH3_loss(k).gt.0.) then
            ! new particle mass
            restmass = xmass1(jpart,k)*exp(-1/nh3lifetime*abs(ltsample))
            if (jpart.lt.-10) then
              write(*,*) 'partpos: ',NH3x,NH3y,NH3z,jpart,ztra1(jpart),xtra1(jpart),ytra1(jpart),itime,ltsample
              write(*,*) 'NH3field: ',NH3LOSS_field(NH3x,NH3y,NH3z),'xmass: ',xmass1(jpart,1),xmass1(jpart,2), &
                   &'restmass: ',restmass,'lifetime: ',nh3lifetime
            endif
            if (restmass .gt. smallnum) then
               xmass1(jpart,k)=restmass
            else
               xmass1(jpart,k)=0.
            endif
         endif
      end do

!     endif  ! nh3loss_average.gt.smallnum 
   endif ! NH3LOSS

  end do  !continue loop over all particles


end subroutine ohreaction

