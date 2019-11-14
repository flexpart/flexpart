! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
  integer :: jjjjmmdd,hhmmss,OHx,OHy,OHz
  real, dimension(nzOH) :: altOHtop
  real :: xlon,ylat 
  real :: xtn,ytn
  real :: restmass,ohreacted,oh_average
  real :: ohrate,temp 
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

    ! get position in the OH field
    OHx=minloc(abs(lonOH-xlon),dim=1,mask=abs(lonOH-xlon).eq.minval(abs(lonOH-xlon)))
    OHy=minloc(abs(latOH-ylat),dim=1,mask=abs(latOH-ylat).eq.minval(abs(latOH-ylat)))

    ! get the level of the OH field for the particle
    ! ztra1 is the z-coord of the trajectory above model orography in metres
    ! altOH is the height of the centre of the level in the OH field above orography
    do i=2,nzOH
      altOHtop(i-1)=altOH(i)+0.5*(altOH(i)-altOH(i-1))
    end do
    altOHtop(nzOH)=altOH(nzOH)+0.5*(altOH(nzOH)-altOH(nzOH-1))
    OHz=minloc(abs(altOHtop-ztra1(jpart)),dim=1,mask=abs(altOHtop-ztra1(jpart))&
            &.eq.minval(abs(altOHtop-ztra1(jpart))))

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

  end do  !continue loop over all particles


end subroutine ohreaction

