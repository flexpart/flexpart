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

subroutine calcpar(n,uuh,vvh,pvh,metdata_format)
  !                   i  i   i   o
  !*****************************************************************************
  !                                                                            *
  !     Computation of several boundary layer parameters needed for the        *
  !     dispersion calculation and calculation of dry deposition velocities.   *
  !     All parameters are calculated over the entire grid.                    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     21 May 1995                                                            *
  !                                                                            *
  ! ------------------------------------------------------------------         *
  !   Petra Seibert, Feb 2000:                                                 *
  !     convection scheme:                                                     *
  !     new variables in call to richardson                                    *
  !                                                                            *
  !   CHANGE 17/11/2005 Caroline Forster NCEP GFS version                      *
  !                                                                            *
  !   Changes, Bernd C. Krueger, Feb. 2001:                                    *
  !    Variables tth and qvh (on eta coordinates) in common block              *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Merged calcpar and calcpar_gfs into one routine using if-then        *
  !       for meteo-type dependent code                                        *
  !*****************************************************************************

  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 3)       *
  ! uuh                                                                        *
  ! vvh                                                                        *
  ! pvh                                                                        *
  ! metdata_format     format of metdata (ecmwf/gfs)                           * 
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !                                                                            *
  ! Functions:                                                                 *
  ! scalev             computation of ustar                                    *
  ! obukhov            computatio of Obukhov length                            *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use class_gribfile

  implicit none

  integer :: metdata_format
  integer :: n,ix,jy,i,kz,lz,kzmin,llev,loop_start
  real :: ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
  real :: ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
  real :: altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax),hmixdummy,akzdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real,parameter :: const=r_air/ga

  !write(*,*) 'in calcpar writting snowheight'
  !***********************************
  ! for test: write out snow depths

  ! open(4,file='slandusetest',form='formatted')
  ! do 5 ix=0,nxmin1
  !5       write (4,*) (sd(ix,jy,1,n),jy=0,nymin1)
  !  close(4)


  ! Loop over entire grid
  !**********************

  do jy=0,nymin1

  ! Set minimum height for tropopause
  !**********************************

    ylat=ylat0+real(jy)*dy
    if ((ylat.ge.-20.).and.(ylat.le.20.)) then
      altmin = 5000.
    else
      if ((ylat.gt.20.).and.(ylat.lt.40.)) then
        altmin=2500.+(40.-ylat)*125.
      else if ((ylat.gt.-40.).and.(ylat.lt.-20.)) then
        altmin=2500.+(40.+ylat)*125.
      else
        altmin=2500.
      endif
    endif

    do ix=0,nxmin1

  ! 1) Calculation of friction velocity
  !************************************

      ustar(ix,jy,1,n)=scalev(ps(ix,jy,1,n),tt2(ix,jy,1,n), &
           td2(ix,jy,1,n),surfstr(ix,jy,1,n))
      if (ustar(ix,jy,1,n).le.1.e-8) ustar(ix,jy,1,n)=1.e-8

  ! 2) Calculation of inverse Obukhov length scale
  !***********************************************

      if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
        ! NCEP version: find first level above ground
        llev = 0
        do i=1,nuvz
          if (ps(ix,jy,1,n).lt.akz(i)) llev=i
        end do
        llev = llev+1
        if (llev.gt.nuvz) llev = nuvz-1
        ! NCEP version

        ! calculate inverse Obukhov length scale with tth(llev)
        ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n), &
             tth(ix,jy,llev,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm,akz(llev),metdata_format)
      else
        llev=0
        ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n), &
            tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm,akzdummy,metdata_format)
      end if

      if (ol.ne.0.) then
        oli(ix,jy,1,n)=1./ol
      else
        oli(ix,jy,1,n)=99999.
      endif


  ! 3) Calculation of convective velocity scale and mixing height
  !**************************************************************

      do i=1,nuvz
        ulev(i)=uuh(ix,jy,i)
        vlev(i)=vvh(ix,jy,i)
        ttlev(i)=tth(ix,jy,i,n)
        qvlev(i)=qvh(ix,jy,i,n)
      end do

      if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
        ! NCEP version hmix has been read in in readwind.f, is therefore not calculated here
      call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev, &
           ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n), &
             td2(ix,jy,1,n),hmixdummy,wstar(ix,jy,1,n),hmixplus,metdata_format)
      else
        call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev, &
             ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n), &
             td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus,metdata_format)
      end if

      if(lsubgrid.eq.1) then
        subsceff=min(excessoro(ix,jy),hmixplus)
      else
        subsceff=0.0
      endif
  !
  ! CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
  !
      hmix(ix,jy,1,n)=hmix(ix,jy,1,n)+subsceff
      hmix(ix,jy,1,n)=max(hmixmin,hmix(ix,jy,1,n)) ! set minimum PBL height
      hmix(ix,jy,1,n)=min(hmixmax,hmix(ix,jy,1,n)) ! set maximum PBL height

  ! 4) Calculation of dry deposition velocities
  !********************************************

      if (DRYDEP) then
  ! Sabine Eckhardt, Dec 06: use new index for z0 for water depending on
  ! windspeed
        z0(7)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga

  ! Calculate relative humidity at surface
  !***************************************
        rh=ew(td2(ix,jy,1,n))/ew(tt2(ix,jy,1,n))

        call getvdep(n,ix,jy,ustar(ix,jy,1,n), &
             tt2(ix,jy,1,n),ps(ix,jy,1,n),1./oli(ix,jy,1,n), &
             ssr(ix,jy,1,n),rh,lsprec(ix,jy,1,n)+convprec(ix,jy,1,n), &
             sd(ix,jy,1,n),vd)

        do i=1,nspec
          vdep(ix,jy,i,n)=vd(i)
        end do

      endif

  !******************************************************
  ! Calculate height of thermal tropopause (Hoinka, 1997)
  !******************************************************

  ! 1) Calculate altitudes of model levels
  !***************************************

      tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
      pold=ps(ix,jy,1,n)
      zold=0.
      if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
        loop_start=2
      else
        loop_start=llev
      end if
      do kz=loop_start,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)  ! pressure on model layers
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))

        if (abs(tv-tvold).gt.0.2) then
         zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          zlev(kz)=zold+const*log(pold/pint)*tv
        endif
        tvold=tv
        pold=pint
        zold=zlev(kz)
      end do

  ! 2) Define a minimum level kzmin, from which upward the tropopause is
  !    searched for. This is to avoid inversions in the lower troposphere
  !    to be identified as the tropopause
  !************************************************************************

      if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
        loop_start=1
      else
        loop_start=llev
      end if

      do kz=loop_start,nuvz
        if (zlev(kz).ge.altmin) then
          kzmin=kz
          goto 45
        endif
      end do
45    continue

  ! 3) Search for first stable layer above minimum height that fulfills the
  !    thermal tropopause criterion
  !************************************************************************

      do kz=kzmin,nuvz
        do lz=kz+1,nuvz
          if ((zlev(lz)-zlev(kz)).gt.2000.) then
            if (((tth(ix,jy,kz,n)-tth(ix,jy,lz,n))/ &
                 (zlev(lz)-zlev(kz))).lt.0.002) then
              tropopause(ix,jy,1,n)=zlev(kz)
              goto 51
            endif
            goto 50
          endif
        end do
50      continue
      end do
51    continue


    end do
  end do

  ! Calculation of potential vorticity on 3-d grid
  !***********************************************

  call calcpv(n,uuh,vvh,pvh)


end subroutine calcpar
