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

subroutine verttransform(n,uuh,vvh,wwh,pvh)
!                         i  i   i   i   i
!*****************************************************************************
!                                                                            *
!     This subroutine transforms temperature, dew point temperature and      *
!     wind components from eta to meter coordinates.                         *
!     The vertical wind component is transformed from Pa/s to m/s using      *
!     the conversion factor pinmconv.                                        *
!     In addition, this routine calculates vertical density gradients        *
!     needed for the parameterization of the turbulent velocities.           *
!                                                                            *
!     Author: A. Stohl, G. Wotawa                                            *
!                                                                            *
!     12 August 1996                                                         *
!     Update: 16 January 1998                                                *
!                                                                            *
!     Major update: 17 February 1999                                         *
!     by G. Wotawa                                                           *
!                                                                            *
!     - Vertical levels for u, v and w are put together                      *
!     - Slope correction for vertical velocity: Modification of calculation  *
!       procedure                                                            *
!                                                                            *
!*****************************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:
!   Variables tth and qvh (on eta coordinates) from common block
!*****************************************************************************
! Sabine Eckhardt, March 2007
! added the variable cloud for use with scavenging - descr. in com_mod
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! nx,ny,nz                        field dimensions in x,y and z direction    *
! clouds(0:nxmax,0:nymax,0:nzmax,numwfmem) cloud field for wet deposition    *
! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
!                                          [deltaeta/s]                      *
! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
!                                                                            *
!*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll
! eso TODO:
! only used for timing of CPU measurement. Remove this (and calls to mpif_mtime below)
! as this routine should not be dependent on MPI
!  use mpi_mod
! :TODO

  implicit none

  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
  integer :: rain_cloud_above(0:nxmax-1,0:nymax-1),kz_inv,idx(0:nxmax-1,0:nymax-1)
  real :: f_qvsat,pressure,cosf(0:nymax-1)
  real :: rh,lsp,convp,tim,tim2,rhmin,precmin,prec
  real :: uvzlev(0:nxmax-1,0:nymax-1,nuvzmax),rhoh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pinmconv(0:nxmax-1,0:nymax-1,nzmax)
  real :: ew,pint(0:nxmax-1,0:nymax-1),tv(0:nxmax-1,0:nymax-1)
  real :: tvold(0:nxmax-1,0:nymax-1),pold(0:nxmax-1,0:nymax-1),dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: wzlev(0:nxmax-1,0:nymax-1,nwzmax)
  real,parameter :: const=r_air/ga
!  integer:: ihr, imin, isec, i100th,ihr2, imin2, isec2, i100th2
  parameter (precmin = 0.002) ! minimum prec in mm/h for cloud diagnostics

  logical :: init = .true.

  !ZHG SEP 2014 tests  
  integer :: cloud_ver,cloud_min, cloud_max 
  integer ::teller(5), convpteller=0, lspteller=0
  real :: cloud_col_wat, cloud_water
  !ZHG 2015 temporary variables for testing
  real :: rcw(0:nxmax-1,0:nymax-1)
  real :: rpc(0:nxmax-1,0:nymax-1)
  character(len=60) :: zhgpath='/xnilu_wrk/flex_wrk/zhg/'
  character(len=60) :: fnameA,fnameB,fnameC,fnameD,fnameE,fnameF,fnameG,fnameH
  CHARACTER(LEN=3)  :: aspec
  integer :: virr=0
  real :: tot_cloud_h
!ZHG

!*************************************************************************
! If verttransform is called the first time, initialize heights of the   *
! z levels in meter. The heights are the heights of model levels, where  *
! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
! the vertical resolution in the z system is doubled. As reference point,*
! the lower left corner of the grid is used.                             *
! Unlike in the eta system, no difference between heights for u,v and    *
! heights for w exists.                                                  *
!*************************************************************************


!eso measure CPU time
!  call mpif_mtime('verttransform',0)

  if (init) then

! Search for a point with high surface pressure (i.e. not above significant topography)
! Then, use this point to construct a reference z profile, to be used at all times
!*****************************************************************************

    do jy=0,nymin1
      do ix=0,nxmin1
        if (ps(ix,jy,1,n).gt.100000.) then
          ixm=ix
          jym=jy
          goto 3
        endif
      end do
    end do
3   continue


    tvold(ixm,jym)=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n))/ &
         ps(ixm,jym,1,n))
    pold(ixm,jym)=ps(ixm,jym,1,n)
    height(1)=0.

    do kz=2,nuvz
      pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
      tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

      if (abs(tv(ixm,jym)-tvold(ixm,jym)).gt.0.2) then
        height(kz)= &
             height(kz-1)+const*log(pold(ixm,jym)/pint(ixm,jym))* &
             (tv(ixm,jym)-tvold(ixm,jym))/log(tv(ixm,jym)/tvold(ixm,jym))
      else
        height(kz)=height(kz-1)+ &
             const*log(pold(ixm,jym)/pint(ixm,jym))*tv(ixm,jym)
      endif

      tvold(ixm,jym)=tv(ixm,jym)
      pold(ixm,jym)=pint(ixm,jym)
    end do


! Determine highest levels that can be within PBL
!************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        goto 9
      endif
    end do
9   continue

! Do not repeat initialization of the Cartesian z grid
!*****************************************************

    init=.false.

  endif


! Loop over the whole grid
!*************************


  do jy=0,nymin1
    do ix=0,nxmin1
      tvold(ix,jy)=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
    enddo
  enddo
  pold=ps(:,:,1,n)
  uvzlev(:,:,1)=0.
  wzlev(:,:,1)=0.
  rhoh(:,:,1)=pold/(r_air*tvold)


! Compute heights of eta levels
!******************************

  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(:,:,1,n)
    tv=tth(:,:,kz,n)*(1.+0.608*qvh(:,:,kz,n))
    rhoh(:,:,kz)=pint(:,:)/(r_air*tv)

    where (abs(tv-tvold).gt.0.2) 
      uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*log(pold/pint)* &
           (tv-tvold)/log(tv/tvold)
    elsewhere
      uvzlev(:,:,kz)=uvzlev(:,:,kz-1)+const*log(pold/pint)*tv
    endwhere

    tvold=tv
    pold=pint
  end do


  do kz=2,nwz-1
    wzlev(:,:,kz)=(uvzlev(:,:,kz+1)+uvzlev(:,:,kz))/2.
  end do
  wzlev(:,:,nwz)=wzlev(:,:,nwz-1)+ &
       uvzlev(:,:,nuvz)-uvzlev(:,:,nuvz-1)

! pinmconv=(h2-h1)/(p2-p1)

  pinmconv(:,:,1)=(uvzlev(:,:,2))/ &
       ((aknew(2)+bknew(2)*ps(:,:,1,n))- &
       (aknew(1)+bknew(1)*ps(:,:,1,n)))
  do kz=2,nz-1
    pinmconv(:,:,kz)=(uvzlev(:,:,kz+1)-uvzlev(:,:,kz-1))/ &
         ((aknew(kz+1)+bknew(kz+1)*ps(:,:,1,n))- &
         (aknew(kz-1)+bknew(kz-1)*ps(:,:,1,n)))
  end do
  pinmconv(:,:,nz)=(uvzlev(:,:,nz)-uvzlev(:,:,nz-1))/ &
       ((aknew(nz)+bknew(nz)*ps(:,:,1,n))- &
       (aknew(nz-1)+bknew(nz-1)*ps(:,:,1,n)))

! Levels, where u,v,t and q are given
!************************************


  uu(:,:,1,n)=uuh(:,:,1)
  vv(:,:,1,n)=vvh(:,:,1)
  tt(:,:,1,n)=tth(:,:,1,n)
  qv(:,:,1,n)=qvh(:,:,1,n)
!hg adding the cloud water 
  clwc(:,:,1,n)=clwch(:,:,1,n)
  if (.not.sumclouds) ciwc(:,:,1,n)=ciwch(:,:,1,n)   
!hg 
  pv(:,:,1,n)=pvh(:,:,1)
  rho(:,:,1,n)=rhoh(:,:,1)
  uu(:,:,nz,n)=uuh(:,:,nuvz)
  vv(:,:,nz,n)=vvh(:,:,nuvz)
  tt(:,:,nz,n)=tth(:,:,nuvz,n)
  qv(:,:,nz,n)=qvh(:,:,nuvz,n)

!hg adding the cloud water
  clwc(:,:,nz,n)=clwch(:,:,nuvz,n)
  if (.not.sumclouds) ciwc(:,:,nz,n)=ciwch(:,:,nuvz,n)
!hg
  pv(:,:,nz,n)=pvh(:,:,nuvz)
  rho(:,:,nz,n)=rhoh(:,:,nuvz)


  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).gt.uvzlev(ix,jy,nuvz)) then
          uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
          vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
          tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
          qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
!hg adding the cloud water
          clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
          if (.not.sumclouds) ciwc(ix,jy,iz,n)=ciwc(ix,jy,nz,n)
!hg
          pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
          rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
        else
          innuvz: do kz=idx(ix,jy),nuvz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                 (height(iz).le.uvzlev(ix,jy,kz))) then
              idx(ix,jy)=kz
              exit innuvz
            endif
          enddo innuvz
        endif
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).le.uvzlev(ix,jy,nuvz)) then
          kz=idx(ix,jy)
          dz1=height(iz)-uvzlev(ix,jy,kz-1)
          dz2=uvzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
          vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
          tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
               +tth(ix,jy,kz,n)*dz1)/dz
          qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
               +qvh(ix,jy,kz,n)*dz1)/dz
!hg adding the cloud water
          clwc(ix,jy,iz,n)=(clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
          if (.not.sumclouds) &
               &ciwc(ix,jy,iz,n)=(ciwch(ix,jy,kz-1,n)*dz2+ciwch(ix,jy,kz,n)*dz1)/dz
!hg
          pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
          rho(ix,jy,iz,n)=(rhoh(ix,jy,kz-1)*dz2+rhoh(ix,jy,kz)*dz1)/dz
        endif
      enddo
    enddo
  enddo


! Levels, where w is given
!*************************

  ww(:,:,1,n)=wwh(:,:,1)*pinmconv(:,:,1)
  ww(:,:,nz,n)=wwh(:,:,nwz)*pinmconv(:,:,nz)
  kmin=2
  idx=kmin
  do iz=2,nz
    do jy=0,nymin1
      do ix=0,nxmin1
        inn:         do kz=idx(ix,jy),nwz
          if(idx(ix,jy) .le. kz .and. height(iz).gt.wzlev(ix,jy,kz-1).and. &
               height(iz).le.wzlev(ix,jy,kz)) then
            idx(ix,jy)=kz
            exit inn
          endif
        enddo inn
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        kz=idx(ix,jy)
        dz1=height(iz)-wzlev(ix,jy,kz-1)
        dz2=wzlev(ix,jy,kz)-height(iz)
        dz=dz1+dz2
        ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(ix,jy,kz-1)*dz2 &
             +wwh(ix,jy,kz)*pinmconv(ix,jy,kz)*dz1)/dz
      enddo
    enddo
  enddo

! Compute density gradients at intermediate levels
!*************************************************

  drhodz(:,:,1,n)=(rho(:,:,2,n)-rho(:,:,1,n))/ &
       (height(2)-height(1))
  do kz=2,nz-1
    drhodz(:,:,kz,n)=(rho(:,:,kz+1,n)-rho(:,:,kz-1,n))/ &
         (height(kz+1)-height(kz-1))
  end do
  drhodz(:,:,nz,n)=drhodz(:,:,nz-1,n)

!    end do
!  end do


!****************************************************************
! Compute slope of eta levels in windward direction and resulting
! vertical wind correction
!****************************************************************

  do jy=1,ny-2
    cosf(jy)=1./cos((real(jy)*dy+ylat0)*pi180)
  enddo

  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=1,ny-2
      do ix=1,nx-2

        inneta: do kz=idx(ix,jy),nz
          if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
               (height(iz).le.uvzlev(ix,jy,kz))) then
            idx(ix,jy)=kz
            exit inneta
          endif
        enddo inneta
      enddo
    enddo

    do jy=1,ny-2
      do ix=1,nx-2
        kz=idx(ix,jy)
        dz1=height(iz)-uvzlev(ix,jy,kz-1)
        dz2=uvzlev(ix,jy,kz)-height(iz)
        dz=dz1+dz2
        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvzlev(ixp,jy,kz-1)-uvzlev(ix1,jy,kz-1))/2.
        dzdx2=(uvzlev(ixp,jy,kz)-uvzlev(ix1,jy,kz))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvzlev(ix,jyp,kz-1)-uvzlev(ix,jy1,kz-1))/2.
        dzdy2=(uvzlev(ix,jyp,kz)-uvzlev(ix,jy1,kz))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*uu(ix,jy,iz,n)*dxconst*cosf(jy)+dzdy*vv(ix,jy,iz,n)*dyconst)

      end do

    end do
  end do

! If north pole is in the domain, calculate wind velocities in polar
! stereographic coordinates
!*******************************************************************

  if (nglobal) then
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do


    do iz=1,nz

! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
!
!   AMSnauffer Nov 18 2004 Added check for case vv=0
!
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+ &
           vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else if (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


! Fix: Set W at pole to the zonally averaged W of the next equator-
! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

  endif


! If south pole is in the domain, calculate wind velocities in polar
! stereographic coordinates
!*******************************************************************

  if (sglobal) then
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do

    do iz=1,nz

! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
!
!   AMSnauffer Nov 18 2004 Added check for case vv=0
!
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+ &
           vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else if (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


! Fix: Set W at pole to the zonally averaged W of the next equator-
! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do
  endif


!***********************************************************************************  
  if (readclouds) then !HG METHOD
! The method is loops all grids vertically and constructs the 3D matrix for clouds
! Cloud top and cloud bottom gid cells are assigned as well as the total column 
! cloud water. For precipitating grids, the type and whether it is in or below 
! cloud scavenging are assigned with numbers 2-5 (following the old metod).
! Distinction is done for lsp and convp though they are treated the same in regards
! to scavenging. Also clouds that are not precipitating are defined which may be 
! to include future cloud processing by non-precipitating-clouds. 
!***********************************************************************************
    write(*,*) 'using cloud water from ECMWF'
    clw(:,:,:,n)=0.0
    icloud_stats(:,:,:,n)=0.0
    clouds(:,:,:,n)=0
! If water/ice are read separately into clwc and ciwc, store sum in clwc
    if (.not.sumclouds) then 
      clwc(:,:,:,n) = clwc(:,:,:,n) + ciwc(:,:,:,n)
    end if
    do jy=0,nymin1
      do ix=0,nxmin1
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        prec=lsp+convp
        tot_cloud_h=0
! Find clouds in the vertical
        do kz=1, nz-1 !go from top to bottom
          if (clwc(ix,jy,kz,n).gt.0) then      
! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3 
            clw(ix,jy,kz,n)=(clwc(ix,jy,kz,n)*rho(ix,jy,kz,n))*(height(kz+1)-height(kz))
            tot_cloud_h=tot_cloud_h+(height(kz+1)-height(kz)) 
            icloud_stats(ix,jy,4,n)= icloud_stats(ix,jy,4,n)+clw(ix,jy,kz,n)          ! Column cloud water [m3/m3]
            icloud_stats(ix,jy,3,n)= min(height(kz+1),height(kz))                     ! Cloud BOT height stats      [m]
!ZHG 2015 extra for testing
            clh(ix,jy,kz,n)=height(kz+1)-height(kz)
            icloud_stats(ix,jy,1,n)=icloud_stats(ix,jy,1,n)+(height(kz+1)-height(kz)) ! Cloud total vertical extent [m]
            icloud_stats(ix,jy,2,n)= max(icloud_stats(ix,jy,2,n),height(kz))          ! Cloud TOP height            [m]
!ZHG
          endif
        end do

! If Precipitation. Define removal type in the vertical
        if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

          do kz=nz,1,-1 !go Bottom up!
            if (clw(ix,jy,kz,n).gt. 0) then ! is in cloud
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1) 
              clouds(ix,jy,kz,n)=1                               ! is a cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3                            ! lsp in-cloud
              else
                clouds(ix,jy,kz,n)=2                             ! convp in-cloud
              endif                                              ! convective or large scale
            elseif((clw(ix,jy,kz,n).le.0) .and. (icloud_stats(ix,jy,3,n).ge.height(kz)) ) then ! is below cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5                             ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4                             ! convp dominated washout
              endif                                              ! convective or large scale 
            endif

            if (height(kz).ge. 19000) then                        ! set a max height for removal
              clouds(ix,jy,kz,n)=0
            endif !clw>0
          end do !nz
        endif ! precipitation
      end do
    end do
!**************************************************************************
  else       ! use old definitions
!**************************************************************************
!   create a cloud and rainout/washout field, clouds occur where rh>80%
!   total cloudheight is stored at level 0
    write(*,*) 'using cloud water from Parameterization'
    do jy=0,nymin1
      do ix=0,nxmin1
! OLD METHOD
        rain_cloud_above(ix,jy)=0
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        cloudsh(ix,jy,n)=0
        do kz_inv=1,nz-1
          kz=nz-kz_inv+1
          pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
          rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
          clouds(ix,jy,kz,n)=0
          if (rh.gt.0.8) then ! in cloud
            if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
              rain_cloud_above(ix,jy)=1
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+ &
                   height(kz)-height(kz-1)
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
              else
                clouds(ix,jy,kz,n)=2 ! convp dominated rainout
              endif
            else ! no precipitation
              clouds(ix,jy,kz,n)=1 ! cloud
            endif
          else ! no cloud
            if (rain_cloud_above(ix,jy).eq.1) then ! scavenging
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5 ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4 ! convp dominated washout
              endif
            endif
          endif
        end do
!END OLD METHOD
      end do
    end do
  endif !readclouds

! eso: copy the relevant data to clw4 to reduce amount of communicated data for MPI
  clw4(:,:,n) = icloud_stats(:,:,4,n)

     !********* TEST ***************
     ! WRITE OUT SOME TEST VARIABLES
     !********* TEST ************'**
!teller(:)=0
!virr=virr+1
!WRITE(aspec, '(i3.3)'), virr

!if (readclouds) then
!fnameH=trim(zhgpath)//trim(aspec)//'Vertical_placement.txt'
!else
!fnameH=trim(zhgpath)//trim(aspec)//'Vertical_placement_old.txt'
!endif
!
!OPEN(UNIT=118, FILE=fnameH,FORM='FORMATTED',STATUS = 'UNKNOWN')
!do kz_inv=1,nz-1
!  kz=nz-kz_inv+1
!  !kz=91
!  do jy=0,nymin1
!     do ix=0,nxmin1
!          if (clouds(ix,jy,kz,n).eq.1) teller(1)=teller(1)+1 ! no precipitation cloud
!          if (clouds(ix,jy,kz,n).eq.2) teller(2)=teller(2)+1 ! convp dominated rainout
!          if (clouds(ix,jy,kz,n).eq.3) teller(3)=teller(3)+1 ! lsp dominated rainout
!          if (clouds(ix,jy,kz,n).eq.4) teller(4)=teller(4)+1 ! convp dominated washout
!          if (clouds(ix,jy,kz,n).eq.5) teller(5)=teller(5)+1 ! lsp dominated washout
!          
!        !  write(*,*) height(kz),teller
!     end do
!  end do
!  write(118,*) height(kz),teller
!  teller(:)=0
!end do
!teller(:)=0
!write(*,*) teller 
!write(*,*) aspec
!
!fnameA=trim(zhgpath)//trim(aspec)//'cloudV.txt'
!fnameB=trim(zhgpath)//trim(aspec)//'cloudT.txt'
!fnameC=trim(zhgpath)//trim(aspec)//'cloudB.txt'
!fnameD=trim(zhgpath)//trim(aspec)//'cloudW.txt'
!fnameE=trim(zhgpath)//trim(aspec)//'old_cloudV.txt'
!fnameF=trim(zhgpath)//trim(aspec)//'lsp.txt'
!fnameG=trim(zhgpath)//trim(aspec)//'convp.txt'
!if (readclouds) then
!OPEN(UNIT=111, FILE=fnameA,FORM='FORMATTED',STATUS = 'UNKNOWN')
!OPEN(UNIT=112, FILE=fnameB,FORM='FORMATTED',STATUS = 'UNKNOWN')
!OPEN(UNIT=113, FILE=fnameC,FORM='FORMATTED',STATUS = 'UNKNOWN')
!OPEN(UNIT=114, FILE=fnameD,FORM='FORMATTED',STATUS = 'UNKNOWN')
!else
!OPEN(UNIT=115, FILE=fnameE,FORM='FORMATTED',STATUS = 'UNKNOWN')
!OPEN(UNIT=116, FILE=fnameF,FORM='FORMATTED',STATUS = 'UNKNOWN')
!OPEN(UNIT=117, FILE=fnameG,FORM='FORMATTED',STATUS = 'UNKNOWN')
!endif
!
!do ix=0,nxmin1
!if (readclouds) then
!write(111,*) (icloud_stats(ix,jy,1,n),jy=0,nymin1)
!write(112,*) (icloud_stats(ix,jy,2,n),jy=0,nymin1)
!write(113,*) (icloud_stats(ix,jy,3,n),jy=0,nymin1)
!write(114,*) (icloud_stats(ix,jy,4,n),jy=0,nymin1)
!else
!write(115,*) (cloudsh(ix,jy,n),jy=0,nymin1)    !integer
!write(116,*) (lsprec(ix,jy,1,n),jy=0,nymin1)   !7.83691406E-02 
!write(117,*) (convprec(ix,jy,1,n),jy=0,nymin1) !5.38330078E-02
!endif
!end do
!
!if (readclouds) then
!CLOSE(111)
!CLOSE(112)
!CLOSE(113)
!CLOSE(114)
!else
!CLOSE(115)
!CLOSE(116)
!CLOSE(117)
!endif
!
!END ********* TEST *************** END
! WRITE OUT SOME TEST VARIABLES
!END ********* TEST *************** END


! PS 2012
!      lsp=lsprec(ix,jy,1,n)
!      convp=convprec(ix,jy,1,n)
!          prec=lsp+convp
!          if (lsp.gt.convp) then !  prectype='lsp'
!            lconvectprec = .false.
!          else ! prectype='cp'
!            lconvectprec = .true.
!           endif
!      else ! windfields does not contain cloud data 
!          rhmin = 0.90 ! standard condition for presence of clouds
!PS       note that original by Sabine Eckhart was 80%
!PS       however, for T<-20 C we consider saturation over ice
!PS       so I think 90% should be enough          
!          icloudbot(ix,jy,n)=icmv
!          icloudtop=icmv ! this is just a local variable
!98        do kz=1,nz
!            pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
!            rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
!ps            if (prec.gt.0.01) print*,'relhum',prec,kz,rh,height(kz)
!            if (rh .gt. rhmin) then
!              if (icloudbot(ix,jy,n) .eq. icmv) then
!                icloudbot(ix,jy,n)=nint(height(kz))
!              endif
!              icloudtop=nint(height(kz)) ! use int to save memory
!            endif
    ! end do
!PS try to get a cloud thicker than 50 m 
!PS if there is at least .01 mm/h  - changed to 0.002 and put into
!PS parameter precpmin        
!          if ((icloudbot(ix,jy,n) .eq. icmv .or. &
!              icloudtop-icloudbot(ix,jy,n) .lt. 50) .and. &
!              prec .gt. precmin) then
!            rhmin = rhmin - 0.05
!            if (rhmin .ge. 0.30) goto 98 ! give up for <= 25% rel.hum.
!   end if

!PS is based on looking at a limited set of comparison data
!          if (lconvectprec .and. icloudtop .lt. 6000 .and. &
!             prec .gt. precmin) then 
!
!            if (convp .lt. 0.1) then
!              icloudbot(ix,jy,n) = 500
!              icloudtop =         8000
!            else
!              icloudbot(ix,jy,n) = 0
!              icloudtop =      10000
!            endif
!          endif
!          if (icloudtop .ne. icmv) then
!            icloudthck(ix,jy,n) = icloudtop-icloudbot(ix,jy,n)
!          else
!            icloudthck(ix,jy,n) = icmv
!          endif
!PS  get rid of too thin clouds      
!          if (icloudthck(ix,jy,n) .lt. 50) then
!            icloudbot(ix,jy,n)=icmv
!            icloudthck(ix,jy,n)=icmv
!          endif
!hg__________________________________
!           rcw(ix,jy)=2E-7*prec**0.36
!           rpc(ix,jy)=prec
!hg end______________________________

!      endif !hg read clouds




!eso measure CPU time
!  call mpif_mtime('verttransform',1)

!eso print out the same measure as in Leo's routine
    ! write(*,*) 'verttransform: ', &
    !      sum(tt(:,:,:,n)*tt(:,:,:,n)), &
    !      sum(uu(:,:,:,n)*uu(:,:,:,n)),sum(vv(:,:,:,n)*vv(:,:,:,n)),&
    !      sum(qv(:,:,:,n)*qv(:,:,:,n)),sum(pv(:,:,:,n)*pv(:,:,:,n)),&
    !      sum(rho(:,:,:,n)*rho(:,:,:,n)),sum(drhodz(:,:,:,n)*drhodz(:,:,:,n)),&
    !      sum(ww(:,:,:,n)*ww(:,:,:,n)), &
    !      sum(clouds(:,:,:,n)), sum(cloudsh(:,:,n)),sum(idx),sum(pinmconv)
end subroutine verttransform

