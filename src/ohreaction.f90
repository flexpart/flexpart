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
  !    Author: S. Eckhardt                                                     *
  !                                                                            *
  !    June 2007                                                               *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Variables:                                                                 *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! ldeltat [s]        interval since radioactive decay was computed           *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! oh_average [mol/m^3]   OH Concentration                                    *
  ! ltsample [s]       interval over which mass is deposited                   *
  !                                                                            *
  !*****************************************************************************

  use oh_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat,j,k,ix,jy
  integer :: ngrid,il,interp_time,n,mm,indz,i
  integer :: jjjjmmdd,ihmmss,OHx,OHy,dOHx,dOHy,OHz
  real :: xtn,ytn,oh_average
  !real oh_diurn_var,sum_ang
  !real zenithangle, ang
  real :: restmass,ohreacted,OHinc
  real :: xlon, ylat, gas_const, act_energy
  real :: ohreact_temp_corr, act_temp
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real(kind=dp) :: jul

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  gas_const=8.314 ! define gas constant
  act_energy=10000 ! activation energy

  !write(*,*) 'OH reaction n:',n,ohreact(1)
  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif


    dOHx=360/(maxxOH-1)
    dOHy=180/(maxyOH-1)

    jul=bdate+real(itime,kind=dp)/86400._dp
    call caldate(jul,jjjjmmdd,ihmmss)
    mm=int((jjjjmmdd-(jjjjmmdd/10000)*10000)/100)

    do jpart=1,numpart

  ! Determine which nesting level to be used
  !*****************************************

    ngrid=0
    do j=numbnests,1,-1
      if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and. &
           (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
        ngrid=j
        goto 23
      endif
    end do
23   continue


  ! Determine nested grid coordinates
  !**********************************

    if (ngrid.gt.0) then
      xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
      ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
    else
      ix=int(xtra1(jpart))
      jy=int(ytra1(jpart))
    endif

  n=2
  if (abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) &
       n=1

  do i=2,nz
    if (height(i).gt.ztra1(jpart)) then
      indz=i-1
      goto 6
    endif
  end do
6   continue

  ! The concentration from the nearest available gridcell is taken
  ! get OH concentration for the specific month and solar angle

  !  write(*,*) OH_field(1,1,1,1),OH_field(10,1,1,10)
  !  write(*,*) OH_field(1,maxxOH-1,maxyOH-1,1)
  !  write(*,*) OH_field(10,maxxOH-1,maxyOH-1,10)
  !  write(*,*) OH_field_height(1,10,4,1),OH_field_height(10,4,10,10)
  !  write(*,*) OH_field_height(1,maxxOH-1,maxyOH-1,1)
  !  write(*,*) OH_field_height(10,maxxOH-1,maxyOH-1,10)
    interp_time=nint(itime-0.5*ltsample)

  ! World coordinates
    xlon=xtra1(jpart)*dx+xlon0
    if (xlon.gt.180) then
       xlon=xlon-360
    endif
    ylat=ytra1(jpart)*dy+ylat0
  ! get position in the OH field - assume that the OH field is global
    OHx=(180+xlon-1)/dOHx
    OHy=(90+ylat-1)/dOHy
  !  sum_ang=0
  ! get the level of the OH height field were the actual particle is in
  ! ztra1 is the z-coordinate of the trajectory above model orography in m
  ! OH_field_height is the heigth of the OH field above orography
      OHz=maxzOH
  ! assume equally distrib. OH field, OH_field_height gives the middle of
  ! the z coordinate
      OHinc=(OH_field_height(3)-OH_field_height(2))/2
      do il=2,maxzOH+1
        if ((OH_field_height(il-1)+OHinc).gt.ztra1(jpart)) goto 26
      end do
26     continue

     OHz=il-1
  !   loop was not interrupted il would be 8 (9-1)
     if (OHz.gt.maxzOH) OHz=7
  !   write (*,*) 'OH height: '
  !    +        ,ztra1(jpart),jpart,OHz,OH_field_height(OHz),OHinc,
  !    +        OH_field_height

    oh_average=OH_field(mm,OHx,OHy,OHz)
    if (oh_average.gt.smallnum) then
  !**********************************************************
  ! if there is noOH concentration no reaction
  ! for performance reason take average concentration and
  ! ignore diurnal variation
  ! do 28 il=1,24
  !      ang=70-zenithangle(ylat,xlon,jul+(24-il)/24.)
  !      if (ang.lt.0) then
  !          ang=0
  !      endif
  !      sum_ang=sum_ang+ang
  !28         enddo
  !    oh_diurn_var=(ang/sum_ang)*(oh_average*24)
  !    oh_average=oh_diurn_var
  !**********************************************************


  !    Computation of the OH reaction
  !**********************************************************
    act_temp=tt(ix,jy,indz,n)

    do k=1,nspec                                  ! loop over species
      if (ohreact(k).gt.0.) then
        ohreact_temp_corr=ohreact(k)*oh_average* &
             exp((act_energy/gas_const)*(1/298.15-1/act_temp))
        ohreacted=xmass1(jpart,k)* &
             (1.-exp(-1*ohreact_temp_corr*abs(ltsample)))
  !      new particle mass:
        restmass = xmass1(jpart,k)-ohreacted
        if (restmass .gt. smallnum) then
          xmass1(jpart,k)=restmass
  !   write (104) xlon,ylat,ztra1(jpart),k,oh_diurn_var,jjjjmmdd,
  !    +               ihmmss,restmass,ohreacted
        else
          xmass1(jpart,k)=0.
        endif
  !      write (*,*) 'restmass: ',restmass
      else
        ohreacted=0.
      endif
    end do

  endif
  !endif OH concentration gt 0
    end do
  !continue loop over all particles

end subroutine ohreaction
