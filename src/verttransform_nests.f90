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

subroutine verttransform_nests(n,uuhn,vvhn,wwhn,pvhn)
!                            i   i    i    i   i
!*****************************************************************************
!                                                                            *
!     This subroutine transforms temperature, dew point temperature and      *
!     wind components from eta to meter coordinates.                         *
!     The vertical wind component is transformed from Pa/s to m/s using      *
!     the conversion factor pinmconv.                                        *
!     In addition, this routine calculates vertical density gradients        *
!     needed for the parameterization of the turbulent velocities.           *
!     It is similar to verttransform, but makes the transformations for      *
!     the nested grids.                                                      *
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
!  Changes, Bernd C. Krueger, Feb. 2001:       (marked "C-cv")
!   Variables tthn and qvhn (on eta coordinates) from common block
!*****************************************************************************
! Sabine Eckhardt, March 2007
! add the variable cloud for use with scavenging - descr. in com_mod
!*****************************************************************************
! ESO, 2016
! -note that divide-by-zero occurs when nxmaxn,nymaxn etc. are larger than 
!  the actual field dimensions
!*****************************************************************************
! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction    *
! uun                             wind components in x-direction [m/s]       *
! vvn                             wind components in y-direction [m/s]       *
! wwn                             wind components in z-direction [deltaeta/s]*
! ttn                             temperature [K]                            *
! pvn                             potential vorticity (pvu)                  *
! psn                             surface pressure [Pa]                      *
!                                                                            *
!*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests) :: uuhn,vvhn,pvhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests) :: wwhn

  real,dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: rhohn,uvzlev,wzlev
  real,dimension(0:nxmaxn-1,0:nymaxn-1,nzmax) :: pinmconv
  real,dimension(0:nxmaxn-1,0:nymaxn-1) :: tvold,pold,pint,tv
  real,dimension(0:nymaxn-1) :: cosf

  integer,dimension(0:nxmaxn-1,0:nymaxn-1) :: rain_cloud_above, idx

  integer :: ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp,kz_inv
  real :: f_qvsat,pressure,rh,lsp,convp,cloudh_min,prec

  real :: ew,dz1,dz2,dz
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real,parameter :: const=r_air/ga
  real :: tot_cloud_h
  integer :: nxm1, nym1

!  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

! Loop over all nests
!********************

  do l=1,numbnests
    nxm1=nxn(l)-1 
    nym1=nyn(l)-1 

! Loop over the whole grid
!*************************

    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
        tvold(ix,jy)=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/ &
             psn(ix,jy,1,n,l))
      end do
    end do

    pold(0:nxm1,0:nym1)=psn(0:nxm1,0:nym1,1,n,l)
    uvzlev(0:nxm1,0:nym1,1)=0.
    wzlev(0:nxm1,0:nym1,1)=0.
    rhohn(0:nxm1,0:nym1,1)=pold(0:nxm1,0:nym1)/(r_air*tvold(0:nxm1,0:nym1))

! Compute heights of eta levels
!******************************

    do kz=2,nuvz
      pint(0:nxm1,0:nym1)=akz(kz)+bkz(kz)*psn(0:nxm1,0:nym1,1,n,l)
      tv(0:nxm1,0:nym1)=tthn(0:nxm1,0:nym1,kz,n,l)*(1.+0.608*qvhn(0:nxm1,0:nym1,kz,n,l))
      rhohn(0:nxm1,0:nym1,kz)=pint(0:nxm1,0:nym1)/(r_air*tv(0:nxm1,0:nym1))

      where (abs(tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1)).gt.0.2) 
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))* &
             (tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1))/&
             &log(tv(0:nxm1,0:nym1)/tvold(0:nxm1,0:nym1))
      elsewhere
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))*tv(0:nxm1,0:nym1)
      endwhere

      tvold(0:nxm1,0:nym1)=tv(0:nxm1,0:nym1)
      pold(0:nxm1,0:nym1)=pint(0:nxm1,0:nym1)

    end do

    do kz=2,nwz-1
      wzlev(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)+uvzlev(0:nxm1,0:nym1,kz))/2.
    end do
    wzlev(0:nxm1,0:nym1,nwz)=wzlev(0:nxm1,0:nym1,nwz-1)+ &
         uvzlev(0:nxm1,0:nym1,nuvz)-uvzlev(0:nxm1,0:nym1,nuvz-1)


    pinmconv(0:nxm1,0:nym1,1)=(uvzlev(0:nxm1,0:nym1,2))/ &
         ((aknew(2)+bknew(2)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(1)+bknew(1)*psn(0:nxm1,0:nym1,1,n,l)))
    do kz=2,nz-1
      pinmconv(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)-uvzlev(0:nxm1,0:nym1,kz-1))/ &
           ((aknew(kz+1)+bknew(kz+1)*psn(0:nxm1,0:nym1,1,n,l))- &
           (aknew(kz-1)+bknew(kz-1)*psn(0:nxm1,0:nym1,1,n,l)))
    end do
    pinmconv(0:nxm1,0:nym1,nz)=(uvzlev(0:nxm1,0:nym1,nz)-uvzlev(0:nxm1,0:nym1,nz-1))/ &
         ((aknew(nz)+bknew(nz)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(nz-1)+bknew(nz-1)*psn(0:nxm1,0:nym1,1,n,l)))

! Levels, where u,v,t and q are given
!************************************

    uun(0:nxm1,0:nym1,1,n,l)=uuhn(0:nxm1,0:nym1,1,l)
    vvn(0:nxm1,0:nym1,1,n,l)=vvhn(0:nxm1,0:nym1,1,l)
    ttn(0:nxm1,0:nym1,1,n,l)=tthn(0:nxm1,0:nym1,1,n,l)
    qvn(0:nxm1,0:nym1,1,n,l)=qvhn(0:nxm1,0:nym1,1,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,1,n,l)=clwchn(0:nxm1,0:nym1,1,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,1,n,l)=ciwchn(0:nxm1,0:nym1,1,n,l)
    end if
    pvn(0:nxm1,0:nym1,1,n,l)=pvhn(0:nxm1,0:nym1,1,l)
    rhon(0:nxm1,0:nym1,1,n,l)=rhohn(0:nxm1,0:nym1,1)

    uun(0:nxm1,0:nym1,nz,n,l)=uuhn(0:nxm1,0:nym1,nuvz,l)
    vvn(0:nxm1,0:nym1,nz,n,l)=vvhn(0:nxm1,0:nym1,nuvz,l)
    ttn(0:nxm1,0:nym1,nz,n,l)=tthn(0:nxm1,0:nym1,nuvz,n,l)
    qvn(0:nxm1,0:nym1,nz,n,l)=qvhn(0:nxm1,0:nym1,nuvz,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,nz,n,l)=clwchn(0:nxm1,0:nym1,nuvz,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,nz,n,l)=ciwchn(0:nxm1,0:nym1,nuvz,n,l)
    end if
    pvn(0:nxm1,0:nym1,nz,n,l)=pvhn(0:nxm1,0:nym1,nuvz,l)
    rhon(0:nxm1,0:nym1,nz,n,l)=rhohn(0:nxm1,0:nym1,nuvz)


    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).gt.uvzlev(ix,jy,nuvz)) then
            uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
            vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
            ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
            qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
!hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=clwcn(ix,jy,nz,n,l)
              if (.not.sumclouds_nest(l)) ciwcn(ix,jy,iz,n,l)=ciwcn(ix,jy,nz,n,l)
            end if
!hg
            pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
            rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
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
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).le.uvzlev(ix,jy,nuvz)) then
            kz=idx(ix,jy)
            dz1=height(iz)-uvzlev(ix,jy,kz-1)
            dz2=uvzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+uuhn(ix,jy,kz,l)*dz1)/dz
            vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+vvhn(ix,jy,kz,l)*dz1)/dz
            ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2 &
                 +tthn(ix,jy,kz,n,l)*dz1)/dz
            qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2 &
                 +qvhn(ix,jy,kz,n,l)*dz1)/dz
!hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=(clwchn(ix,jy,kz-1,n,l)*dz2+clwchn(ix,jy,kz,n,l)*dz1)/dz
              if (.not.sumclouds_nest(l)) &
                   &ciwcn(ix,jy,iz,n,l)=(ciwchn(ix,jy,kz-1,n,l)*dz2+ciwchn(ix,jy,kz,n,l)*dz1)/dz
            end if
!hg
            pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+pvhn(ix,jy,kz,l)*dz1)/dz
            rhon(ix,jy,iz,n,l)=(rhohn(ix,jy,kz-1)*dz2+rhohn(ix,jy,kz)*dz1)/dz
          endif
        enddo
      enddo
    enddo



!         do kz=kmin,nuvz
!           if(height(iz).gt.uvwzlev(ix,jy,nuvz)) then
!             uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
!             vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
!             ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
!             qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
!             pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
!             rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
!             goto 30
!           endif
!           if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
!           (height(iz).le.uvwzlev(ix,jy,kz))) then
!            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
!            dz2=uvwzlev(ix,jy,kz)-height(iz)
!            dz=dz1+dz2
!            uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+ &
!            uuhn(ix,jy,kz,l)*dz1)/dz
!            vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+ &
!            vvhn(ix,jy,kz,l)*dz1)/dz
!            ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2+ &
!            tthn(ix,jy,kz,n,l)*dz1)/dz
!            qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2+ &
!            qvhn(ix,jy,kz,n,l)*dz1)/dz
!            pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+ &
!            pvhn(ix,jy,kz,l)*dz1)/dz
!            rhon(ix,jy,iz,n,l)=(rhohn(kz-1)*dz2+rhohn(kz)*dz1)/dz
!            kmin=kz
!            goto 30
!           endif
!         end do
! 30      continue
!       end do


! Levels, where w is given
!*************************

    wwn(0:nxm1,0:nym1,1,n,l)=wwhn(0:nxm1,0:nym1,1,l)*pinmconv(0:nxm1,0:nym1,1)
    wwn(0:nxm1,0:nym1,nz,n,l)=wwhn(0:nxm1,0:nym1,nwz,l)*pinmconv(0:nxm1,0:nym1,nz)
    kmin=2
    idx=kmin
    do iz=2,nz
      do jy=0,nym1
        do ix=0,nxm1
          inn:         do kz=idx(ix,jy),nwz
            if(idx(ix,jy) .le. kz .and. height(iz).gt.wzlev(ix,jy,kz-1).and. &
                 height(iz).le.wzlev(ix,jy,kz)) then
              idx(ix,jy)=kz
              exit inn
            endif
          enddo inn
        enddo
      enddo
      do jy=0,nym1
        do ix=0,nxm1
          kz=idx(ix,jy)
          dz1=height(iz)-wzlev(ix,jy,kz-1)
          dz2=wzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(ix,jy,kz-1)*dz2 &
               +wwhn(ix,jy,kz,l)*pinmconv(ix,jy,kz)*dz1)/dz
        enddo
      enddo
    enddo

!       wwn(ix,jy,1,n,l)=wwhn(ix,jy,1,l)*pinmconv(1)
!       wwn(ix,jy,nz,n,l)=wwhn(ix,jy,nwz,l)*pinmconv(nz)
!       kmin=2
!       do iz=2,nz
!         do kz=kmin,nwz
!           if ((height(iz).gt.wzlev(kz-1)).and. &
!           (height(iz).le.wzlev(kz))) then
!             dz1=height(iz)-wzlev(kz-1)
!             dz2=wzlev(kz)-height(iz)
!             dz=dz1+dz2
!             wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(kz-1)*dz2 &
!             +wwhn(ix,jy,kz,l)*pinmconv(kz)*dz1)/dz
!             kmin=kz
!             goto 40
!           endif
!         end do
! 40      continue
!       end do

! Compute density gradients at intermediate levels
!*************************************************

    drhodzn(0:nxm1,0:nym1,1,n,l)=(rhon(0:nxm1,0:nym1,2,n,l)-rhon(0:nxm1,0:nym1,1,n,l))/ &
         (height(2)-height(1))
    do kz=2,nz-1
      drhodzn(0:nxm1,0:nym1,kz,n,l)=(rhon(0:nxm1,0:nym1,kz+1,n,l)-rhon(0:nxm1,0:nym1,kz-1,n,l))/ &
           (height(kz+1)-height(kz-1))
    end do
    drhodzn(0:nxm1,0:nym1,nz,n,l)=drhodzn(0:nxm1,0:nym1,nz-1,n,l)


! drhodzn(ix,jy,1,n,l)=(rhon(ix,jy,2,n,l)-rhon(ix,jy,1,n,l))/ &
! (height(2)-height(1))
! do kz=2,nz-1
!   drhodzn(ix,jy,kz,n,l)=(rhon(ix,jy,kz+1,n,l)- &
!   rhon(ix,jy,kz-1,n,l))/(height(kz+1)-height(kz-1))
! end do
! drhodzn(ix,jy,nz,n,l)=drhodzn(ix,jy,nz-1,n,l)



!   end do
! end do


!****************************************************************
! Compute slope of eta levels in windward direction and resulting
! vertical wind correction
!****************************************************************

    do jy=1,nyn(l)-2
!    cosf=cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
      cosf(jy)=1./cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
    end do

    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2

          inneta: do kz=idx(ix,jy),nz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                 (height(iz).le.uvzlev(ix,jy,kz))) then
              idx(ix,jy)=kz
              exit inneta
            endif
          enddo inneta
        enddo
      enddo

      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2
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

          wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*uun(ix,jy,iz,n,l)*dxconst*xresoln(l)*cosf(jy)+&
               &dzdy*vvn(ix,jy,iz,n,l)*dyconst*yresoln(l))

        end do

      end do
    end do


!   do jy=1,nyn(l)-2
!     do ix=1,nxn(l)-2
!       kmin=2
!       do iz=2,nz-1

!         ui=uun(ix,jy,iz,n,l)*dxconst*xresoln(l)/cosf(jy)
!         vi=vvn(ix,jy,iz,n,l)*dyconst*yresoln(l)

!         do kz=kmin,nz
!           if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
!                (height(iz).le.uvwzlev(ix,jy,kz))) then
!             dz1=height(iz)-uvwzlev(ix,jy,kz-1)
!             dz2=uvwzlev(ix,jy,kz)-height(iz)
!             dz=dz1+dz2
!             kl=kz-1
!             klp=kz
!             kmin=kz
!             goto 47
!           endif
!         end do

! 47      ix1=ix-1
!         jy1=jy-1
!         ixp=ix+1
!         jyp=jy+1

!         dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
!         dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
!         dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

!         dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
!         dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
!         dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

!         wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*ui+dzdy*vi)

!       end do

!     end do
!   end do


!write (*,*) 'initializing nested cloudsn, n:',n
!   create a cloud and rainout/washout field, cloudsn occur where rh>80%


!***********************************************************************************  
    if (readclouds_nest(l)) then !HG METHOD
! The method is loops all grids vertically and constructs the 3D matrix for clouds
! Cloud top and cloud bottom gid cells are assigned as well as the total column 
! cloud water. For precipitating grids, the type and whether it is in or below 
! cloud scavenging are assigned with numbers 2-5 (following the old metod).
! Distinction is done for lsp and convp though they are treated the same in regards
! to scavenging. Also clouds that are not precipitating are defined which may be 
! to include future cloud processing by non-precipitating-clouds. 
!***********************************************************************************
      write(*,*) 'Nested ECMWF fields: using cloud water'
      clwn(0:nxm1,0:nym1,:,n,l)=0.0
!    icloud_stats(0:nxm1,0:nym1,:,n)=0.0
      ctwcn(0:nxm1,0:nym1,n,l)=0.0
      cloudsn(0:nxm1,0:nym1,:,n,l)=0
! If water/ice are read separately into clwc and ciwc, store sum in clwcn
      if (.not.sumclouds_nest(l)) then 
        clwcn(0:nxm1,0:nym1,:,n,l) = clwcn(0:nxm1,0:nym1,:,n,l) + ciwcn(0:nxm1,0:nym1,:,n,l)
      end if
      do jy=0,nym1
        do ix=0,nxm1
          lsp=lsprecn(ix,jy,1,n,l)
          convp=convprecn(ix,jy,1,n,l)
          prec=lsp+convp
          tot_cloud_h=0
! Find clouds in the vertical
          do kz=1, nz-1 !go from top to bottom
            if (clwcn(ix,jy,kz,n,l).gt.0) then      
! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3 
              clwn(ix,jy,kz,n,l)=(clwcn(ix,jy,kz,n,l)*rhon(ix,jy,kz,n,l))*(height(kz+1)-height(kz))
              tot_cloud_h=tot_cloud_h+(height(kz+1)-height(kz)) 
              ctwcn(ix,jy,n,l) = ctwcn(ix,jy,n,l)+clwn(ix,jy,kz,n,l)
!            icloud_stats(ix,jy,4,n)= icloud_stats(ix,jy,4,n)+clw(ix,jy,kz,n)          ! Column cloud water [m3/m3]
!           icloud_stats(ix,jy,3,n)= min(height(kz+1),height(kz))                     ! Cloud BOT height stats      [m]
              cloudh_min=min(height(kz+1),height(kz))
!ZHG 2015 extra for testing
!            clh(ix,jy,kz,n)=height(kz+1)-height(kz)
!            icloud_stats(ix,jy,1,n)=icloud_stats(ix,jy,1,n)+(height(kz+1)-height(kz)) ! Cloud total vertical extent [m]
!            icloud_stats(ix,jy,2,n)= max(icloud_stats(ix,jy,2,n),height(kz))          ! Cloud TOP height            [m]
!ZHG
            endif
          end do

! If Precipitation. Define removal type in the vertical
          if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

            do kz=nz,1,-1 !go Bottom up!
              if (clwn(ix,jy,kz,n,l).gt.0.0) then ! is in cloud
                cloudshn(ix,jy,n,l)=cloudshn(ix,jy,n,l)+height(kz)-height(kz-1) 
                cloudsn(ix,jy,kz,n,l)=1                               ! is a cloud
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=3                            ! lsp in-cloud
                else
                  cloudsn(ix,jy,kz,n,l)=2                             ! convp in-cloud
                endif                                              ! convective or large scale
              else if((clwn(ix,jy,kz,n,l).le.0.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=5                             ! lsp dominated washout
                else
                  cloudsn(ix,jy,kz,n,l)=4                             ! convp dominated washout
                endif                                              ! convective or large scale 
              endif

              if (height(kz).ge. 19000) then                        ! set a max height for removal
                cloudsn(ix,jy,kz,n,l)=0
              endif !clw>0
            end do !nz
          endif ! precipitation
        end do
      end do
!********************************************************************
    else ! old method:
!********************************************************************
      write(*,*) 'Nested fields: using cloud water from Parameterization'
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          rain_cloud_above(ix,jy)=0
          lsp=lsprecn(ix,jy,1,n,l)
          convp=convprecn(ix,jy,1,n,l)
          cloudshn(ix,jy,n,l)=0
          do kz_inv=1,nz-1
            kz=nz-kz_inv+1
            pressure=rhon(ix,jy,kz,n,l)*r_air*ttn(ix,jy,kz,n,l)
            rh=qvn(ix,jy,kz,n,l)/f_qvsat(pressure,ttn(ix,jy,kz,n,l))
            cloudsn(ix,jy,kz,n,l)=0
            if (rh.gt.0.8) then ! in cloud
              if ((lsp.gt.0.01).or.(convp.gt.0.01)) then
                rain_cloud_above(ix,jy)=1
                cloudshn(ix,jy,n,l)=cloudshn(ix,jy,n,l)+ &
                     height(kz)-height(kz-1)
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=3 ! lsp dominated rainout
                else
                  cloudsn(ix,jy,kz,n,l)=2 ! convp dominated rainout
                endif
              else ! no precipitation
                cloudsn(ix,jy,kz,n,l)=1 ! cloud
              endif
            else ! no cloud
              if (rain_cloud_above(ix,jy).eq.1) then ! scavenging
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=5 ! lsp dominated washout
                else
                  cloudsn(ix,jy,kz,n,l)=4 ! convp dominated washout
                endif
              endif
            endif
          end do
        end do
      end do
    end if

  end do ! end loop over nests

end subroutine verttransform_nests
