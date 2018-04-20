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

subroutine gethourlyOH(itime)
  !                     i     
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
  !                                                                            *
  !*****************************************************************************

  use oh_mod
  use par_mod
  use com_mod

  implicit none
  
  integer :: itime
  integer :: ix,jy,kz,m1,m2,d1,d2,m,dd,ndays
  integer :: ijx,jjy
  integer :: jjjjmmdd,hhmmss
  real :: sza,jrate,photo_O1D,zenithangle,OH_tmp,jrate_tmp
  real(kind=dp) :: jul1,jul2
  integer, dimension(12) :: monthdays

  monthdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)

!  print*, 'itime: ',itime
!  print*, 'memOHtime(1):',memOHtime(1)
!  print*, 'memOHtime(2):',memOHtime(2)

!! testing OH
!  open(1000,file=path(2)(1:length(2))//'OH_50N_0E.txt',action='write',status='old',access='append')
!  open(999,file=path(2)(1:length(2))//'jrate_50N_0E.txt',action='write',status='old',access='append')
!  open(998,file=path(2)(1:length(2))//'OH_50S_0E.txt',action='write',status='old',access='append')
!  open(997,file=path(2)(1:length(2))//'jrate_50S_0E.txt',action='write',status='old',access='append')


  ! Check hourly OH field is available for the current time step
  !**************************************************************

  if ((ldirect*memOHtime(1).le.ldirect*itime).and. &
       (ldirect*memOHtime(2).gt.ldirect*itime)) then

  ! The right OH fields are already in memory -> don't do anything
  !****************************************************************

    continue

  else if ((ldirect*memOHtime(2).le.ldirect*itime).and. &
       (memOHtime(2).ne.0.)) then

    ! Current time is after 2nd OH field
    !************************************

    memOHtime(1)=memOHtime(2)
    memOHtime(2)=memOHtime(1)+ldirect*3600.
    OH_hourly(:,:,:,1)=OH_hourly(:,:,:,2)

    ! Compute new hourly value of OH
    !**********************************************************

    jul2=bdate+memOHtime(2)/86400._dp  ! date for next hour
    call caldate(jul2,jjjjmmdd,hhmmss)
    m2=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
    dd=(jjjjmmdd-(jjjjmmdd/100)*100)
    ndays=monthdays(m2)

!    print*, 'jul2:',jul2
!    print*, 'm2:',m2
!    print*, 'dd:',dd
!    print*, 'ndays:',ndays

    do kz=1,nzOH
      do jy=1,nyOH
        do ix=1,nxOH
          ijx=minloc(abs(lonjr-lonOH(ix)),dim=1,mask=abs(lonjr-lonOH(ix)).eq.minval(abs(lonjr-lonOH(ix))))
          jjy=minloc(abs(latjr-latOH(jy)),dim=1,mask=abs(latjr-latOH(jy)).eq.minval(abs(latjr-latOH(jy))))
          ! calculate solar zenith angle in degrees (sza) 
          sza=zenithangle(latOH(jy),lonOH(ix),jul2)
          ! calculate J(O1D) (jrate)
          jrate=photo_O1D(sza)
          ! interpolate monthly fields to current day
          if(dd.ge.(ndays/2)) then
            m=min((m2+1),12)
            OH_tmp=OH_field(ix,jy,kz,m2)+(OH_field(ix,jy,kz,m)-OH_field(ix,jy,kz,m2))*real(dd-ndays/2)/real(ndays)
            jrate_tmp=jrate_average(ijx,jjy,m2)+(jrate_average(ijx,jjy,m)-jrate_average(ijx,jjy,m2))*real(dd-ndays/2)/real(ndays)
          else
            m=max((m2-1),1)
            OH_tmp=OH_field(ix,jy,kz,m)+(OH_field(ix,jy,kz,m2)-OH_field(ix,jy,kz,m))*real(dd+ndays/2)/real(ndays)
            jrate_tmp=jrate_average(ijx,jjy,m)+(jrate_average(ijx,jjy,m2)-jrate_average(ijx,jjy,m))*real(dd+ndays/2)/real(ndays)
          endif
          ! apply hourly correction to OH
          if(jrate_tmp.gt.0.) then
            OH_hourly(ix,jy,kz,2)=OH_tmp*jrate/jrate_tmp
          else
            OH_hourly(ix,jy,kz,2)=0.
          endif
          !! testing OH
!           if(jy.eq.36.and.ix.eq.37.and.kz.eq.1) then
!             write(999,fmt='(E12.6)') OH_hourly(ix,jy,kz,2) 
!             write(1000,fmt='(E12.6)') OH_tmp
!             write(999,fmt='(E12.6)') jrate/jrate_tmp
!           endif
!           if(jy.eq.11.and.ix.eq.37.and.kz.eq.1) then
!             write(998,fmt='(E12.6)') OH_hourly(ix,jy,kz,2)
!             write(998,fmt='(E12.6)') OH_tmp
!             write(997,fmt='(E12.6)') jrate/jrate_tmp
!           endif
        end do
      end do
    end do

  else

    ! No OH fields in memory -> compute both hourly OH fields
    !**********************************************************

    jul1=bdate  ! begin date of simulation (julian)
    call caldate(jul1,jjjjmmdd,hhmmss)
    m1=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
    memOHtime(1)=0.
    d1=(jjjjmmdd-(jjjjmmdd/100)*100)

    jul2=bdate+ldirect*real(1./24.,kind=dp)  ! date for next hour
    call caldate(jul2,jjjjmmdd,hhmmss)
    m2=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
    memOHtime(2)=ldirect*3600.
    d2=(jjjjmmdd-(jjjjmmdd/100)*100)

!    print*, 'jul1:',jul1
!    print*, 'jul2:',jul2
!    print*, 'm1,m2:',m1,m2
!    print*, 'd1,d2:',d1,d2

    do kz=1,nzOH
      do jy=1,nyOH
        do ix=1,nxOH
          ijx=minloc(abs(lonjr-lonOH(ix)),dim=1,mask=abs(lonjr-lonOH(ix)).eq.minval(abs(lonjr-lonOH(ix))))
          jjy=minloc(abs(latjr-latOH(jy)),dim=1,mask=abs(latjr-latOH(jy)).eq.minval(abs(latjr-latOH(jy))))
          ! calculate solar zenith angle in degrees (sza), beginning 
          sza=zenithangle(latOH(jy),lonOH(ix),jul1)
          ! calculate J(O1D) (jrate), beginning
          jrate=photo_O1D(sza)
          ! interpolate monthly fields to current day
          ndays=monthdays(m1)
          if(d1.ge.(ndays/2)) then
            m=min((m1+1),12)
            OH_tmp=OH_field(ix,jy,kz,m1)+(OH_field(ix,jy,kz,m)-OH_field(ix,jy,kz,m1))*real(d1-ndays/2)/real(ndays)
          else
            m=max((m1-1),1)
            OH_tmp=OH_field(ix,jy,kz,m)+(OH_field(ix,jy,kz,m1)-OH_field(ix,jy,kz,m))*real(d1+ndays/2)/real(ndays)
          endif
          ! apply hourly correction to OH
          if(jrate_average(ijx,jjy,m1).gt.0.) then
            OH_hourly(ix,jy,kz,1)=OH_tmp*jrate/jrate_average(ijx,jjy,m1)
          else
            OH_hourly(ix,jy,kz,1)=0.
          endif
          ! calculate solar zenith angle in degrees (sza), after 1-hour 
          sza=zenithangle(latOH(jy),lonOH(ix),jul2)
          ! calculate J(O1D) (jrate), after 1-hour
          jrate=photo_O1D(sza)
          ! interpolate monthly fields to current day
          ndays=monthdays(m2)
          if(d2.ge.(ndays/2)) then
            m=min((m2+1),12)
            OH_tmp=OH_field(ix,jy,kz,m2)+(OH_field(ix,jy,kz,m)-OH_field(ix,jy,kz,m2))*real(d2-ndays/2)/real(ndays)
          else
            m=max((m2-1),1)
            OH_tmp=OH_field(ix,jy,kz,m)+(OH_field(ix,jy,kz,m2)-OH_field(ix,jy,kz,m))*real(d2+ndays/2)/real(ndays)
          endif
          ! apply hourly correction to OH
          if(jrate_average(ijx,jjy,m2).gt.0.) then
            OH_hourly(ix,jy,kz,2)=OH_tmp*jrate/jrate_average(ijx,jjy,m2)
          else
            OH_hourly(ix,jy,kz,2)=0.
          endif
        end do
      end do
    end do

  endif

  !! testing OH
!  close(999)
!  close(998)


end subroutine gethourlyOH

