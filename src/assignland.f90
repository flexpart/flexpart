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

subroutine assignland

  !*****************************************************************************
  !                                                                            *
  !     This routine assigns fractions of the 13 landuse classes to each ECMWF *
  !     grid point.                                                            *
  !     The landuse inventory of                                               *
  !                                                                            *
  ! Belward, A.S., Estes, J.E., and Kline, K.D., 1999,                         *
  ! The IGBP-DIS 1-Km Land-Cover Data Set DISCover:                            *
  ! A Project Overview: Photogrammetric Engineering and Remote Sensing ,       *
  ! v. 65, no. 9, p. 1013-1020                                                 *
  !                                                                            *
  !     if there are no data in the inventory                                  *
  !     the ECMWF land/sea mask is used to distinguish                         *
  !     between sea (-> ocean) and land (-> grasslands).                       *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     5 December 1996                                                        *
  !     8 February 1999 Additional use of nests, A. Stohl                      *
  !    29 December 2006 new landuse inventory, S. Eckhardt                     *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! xlanduse          fractions of numclass landuses for each model grid point *
  ! landinvent       landuse inventory (0.3 deg resolution)                    *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: ix,jy,k,l,li,nrefine,iix,jjy
  integer,parameter :: lumaxx=1200,lumaxy=600
  integer,parameter :: xlon0lu=-180,ylat0lu=-90
  real,parameter :: dxlu=0.3
  real :: xlon,ylat,sumperc,p,xi,yj
  real :: xlandusep(lumaxx,lumaxy,numclass)
   character*2 ck

  do ix=1,lumaxx
    do jy=1,lumaxy
          do k=1,numclass
            xlandusep(ix,jy,k)=0.
          end do
          sumperc=0.
          do li=1,3
            sumperc=sumperc+landinvent(ix,jy,li+3)
          end do
          do li=1,3
            k=landinvent(ix,jy,li)
          if (sumperc.gt.0) then
            p=landinvent(ix,jy,li+3)/sumperc
          else
            p=0
          endif
  ! p has values between 0 and 1
            xlandusep(ix,jy,k)=p
          end do
    end do
  end do

  ! do 13 k=1,11
  ! write (ck,'(i2.2)') k
  ! open(4,file='xlandusetest'//ck,form='formatted')
  ! do 11 ix=1,lumaxx
  !11       write (4,*) (xlandusep(ix,jy,k),jy=1,lumaxy)
  !11       write (4,*) (landinvent(ix,jy,k),jy=1,lumaxy)
  !13     close(4)

  ! write (*,*) xlon0,ylat0,xlon0n(1),ylat0n(1),nxmin1,nymin1
  ! write (*,*) dx, dy, dxout, dyout, ylat0, xlon0
  nrefine=10
  do ix=0,nxmin1
    do jy=0,nymin1
      do k=1,numclass
        sumperc=0.
        xlanduse(ix,jy,k)=0.
      end do
        do iix=1, nrefine
          xlon=(ix+(iix-1)/real(nrefine))*dx+xlon0        ! longitude, should be between -180 and 179
          if (xlon.ge.(xlon0lu+lumaxx*dxlu))  then
               xlon=xlon-lumaxx*dxlu
          endif
          do jjy=1, nrefine
           ylat=(jy+(jjy-1)/real(nrefine))*dy+ylat0       ! and lat. of each gridpoint
           xi=int((xlon-xlon0lu)/dxlu)+1
           yj=int((ylat-ylat0lu)/dxlu)+1
           if (xi.gt.lumaxx) xi=xi-lumaxx
           if (yj.gt.lumaxy) yj=yj-lumaxy
           if (xi.lt.0) then
              write (*,*) 'problem with landuseinv sampling: ', &
                   xlon,xlon0lu,ix,iix,xlon0,dx,nxmax
              stop
           endif
           do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)+xlandusep(int(xi),int(yj),k)
             sumperc=sumperc+xlanduse(ix,jy,k)  ! just for the check if landuseinv. is available
           end do
          end do
        end do
          if (sumperc.gt.0) then                       ! detailed landuse available
          sumperc=0.
          do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)/real(nrefine*nrefine)
            sumperc=sumperc+xlanduse(ix,jy,k)
          end do
  !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
          if (sumperc.lt.1-1E-5) then
            do k=1,numclass
              xlanduse(ix,jy,k)= &
                   xlanduse(ix,jy,k)/sumperc
            end do
          endif
          else
            if (lsm(ix,jy).lt.0.1) then           ! over sea  -> ocean
              xlanduse(ix,jy,3)=1.
            else                                  ! over land -> rangeland
              xlanduse(ix,jy,7)=1.
            endif
          endif


    end do
  end do

  !***********************************
   !for test: write out xlanduse

  ! open(4,file='../landusetest',form='formatted')
  ! do 56 k=1,13
  ! do 55 ix=0,nxmin1
  !55    write (4,*) (xlanduse(ix,jy,k),jy=0,nymin1)
  !56    continue
  ! close(4)
  ! write (*,*) 'landuse written'
  ! pause
  !stop
  !  open(4,file='../landseatest'//ck,form='formatted')
  ! do 57 ix=0,nxmin1
  !57       write (4,*) (lsm(ix,jy),jy=0,nymin1)
  !  write (*,*) 'landseamask written'

  !****************************************
  ! Same as above, but for the nested grids
  !****************************************

  !************** TEST ********************
  ! dyn(1)=dyn(1)/40
  ! dxn(1)=dxn(1)/40
  ! xlon0n(1)=1
  ! ylat0n(1)=50
  !************** TEST ********************

  do l=1,numbnests
    do ix=0,nxn(l)-1
      do jy=0,nyn(l)-1
        do k=1,numclass
          sumperc=0.
          xlandusen(ix,jy,k,l)=0.
        end do
          do iix=1, nrefine
           xlon=(ix+(iix-1)/real(nrefine))*dxn(l)+xlon0n(l)
           do jjy=1, nrefine
             ylat=(jy+(jjy-1)/real(nrefine))*dyn(l)+ylat0n(l)
             xi=int((xlon-xlon0lu)/dxlu)+1
             yj=int((ylat-ylat0lu)/dxlu)+1
             if (xi.gt.lumaxx) xi=xi-lumaxx
             if (yj.gt.lumaxy) yj=yj-lumaxy
             do k=1,numclass
                xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)+ &
                     xlandusep(int(xi),int(yj),k)
               sumperc=sumperc+xlandusen(ix,jy,k,l)
             end do
           end do
          end do
          if (sumperc.gt.0) then                     ! detailed landuse available
          sumperc=0.
            do k=1,numclass
               xlandusen(ix,jy,k,l)= &
                    xlandusen(ix,jy,k,l)/real(nrefine*nrefine)
              sumperc=sumperc+xlandusen(ix,jy,k,l)
            end do
  !cc the sum of all categories should be 1 ... 100 percent ... in order to get vdep right!
            if (sumperc.lt.1-1E-5) then
              do k=1,numclass
                xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)/sumperc
              end do
            endif
          else                                    ! check land/sea mask
            if (lsmn(ix,jy,l).lt.0.1) then   ! over sea  -> ocean
              xlandusen(ix,jy,3,l)=1.
            else                        ! over land -> grasslands
              xlandusen(ix,jy,7,l)=1.
            endif
          endif
      end do
    end do
  end do

  !***********************************
  ! for test: write out xlanduse

  ! do 66 k=1,11
  ! write (ck,'(i2.2)') k
  ! open(4,file='nlandusetest'//ck,form='formatted')
  ! do 65 ix=0,nxn(1)-1
  !65       write (4,*) (xlandusen(ix,jy,k,1),jy=0,nyn(1)-1)
  !66      close(4)

  ! write (*,*) 'landuse nested written'


end subroutine assignland
