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

  integer :: ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp
  integer :: rain_cloud_above,kz_inv
  real :: f_qvsat,pressure,rh,lsp,convp
  real :: wzlev(nwzmax),rhoh(nuvzmax),pinmconv(nzmax)
  real :: uvwzlev(0:nxmaxn-1,0:nymaxn-1,nzmax)
  real :: ew,pint,tv,tvold,pold,dz1,dz2,dz,ui,vi,cosf
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  real,parameter :: const=r_air/ga


  ! Loop over all nests
  !********************

  do l=1,numbnests

  ! Loop over the whole grid
  !*************************

  do jy=0,nyn(l)-1
    do ix=0,nxn(l)-1

      tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/ &
      psn(ix,jy,1,n,l))
      pold=psn(ix,jy,1,n,l)
      uvwzlev(ix,jy,1)=0.
      wzlev(1)=0.
      rhoh(1)=pold/(r_air*tvold)


  ! Compute heights of eta levels
  !******************************

      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)
        tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))
        rhoh(kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)* &
          (tv-tvold)/log(tv/tvold)
        else
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)*tv
        endif

        tvold=tv
        pold=pint
      end do


      do kz=2,nwz-1
        wzlev(kz)=(uvwzlev(ix,jy,kz+1)+uvwzlev(ix,jy,kz))/2.
      end do
      wzlev(nwz)=wzlev(nwz-1)+uvwzlev(ix,jy,nuvz)-uvwzlev(ix,jy,nuvz-1)

  ! pinmconv=(h2-h1)/(p2-p1)

      pinmconv(1)=(uvwzlev(ix,jy,2)-uvwzlev(ix,jy,1))/ &
      ((aknew(2)+bknew(2)*psn(ix,jy,1,n,l))- &
      (aknew(1)+bknew(1)*psn(ix,jy,1,n,l)))
      do kz=2,nz-1
        pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/ &
        ((aknew(kz+1)+bknew(kz+1)*psn(ix,jy,1,n,l))- &
        (aknew(kz-1)+bknew(kz-1)*psn(ix,jy,1,n,l)))
      end do
      pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/ &
      ((aknew(nz)+bknew(nz)*psn(ix,jy,1,n,l))- &
      (aknew(nz-1)+bknew(nz-1)*psn(ix,jy,1,n,l)))


  ! Levels, where u,v,t and q are given
  !************************************

      uun(ix,jy,1,n,l)=uuhn(ix,jy,1,l)
      vvn(ix,jy,1,n,l)=vvhn(ix,jy,1,l)
      ttn(ix,jy,1,n,l)=tthn(ix,jy,1,n,l)
      qvn(ix,jy,1,n,l)=qvhn(ix,jy,1,n,l)
      pvn(ix,jy,1,n,l)=pvhn(ix,jy,1,l)
      rhon(ix,jy,1,n,l)=rhoh(1)
      uun(ix,jy,nz,n,l)=uuhn(ix,jy,nuvz,l)
      vvn(ix,jy,nz,n,l)=vvhn(ix,jy,nuvz,l)
      ttn(ix,jy,nz,n,l)=tthn(ix,jy,nuvz,n,l)
      qvn(ix,jy,nz,n,l)=qvhn(ix,jy,nuvz,n,l)
      pvn(ix,jy,nz,n,l)=pvhn(ix,jy,nuvz,l)
      rhon(ix,jy,nz,n,l)=rhoh(nuvz)
      kmin=2
      do iz=2,nz-1
        do kz=kmin,nuvz
          if(height(iz).gt.uvwzlev(ix,jy,nuvz)) then
            uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
            vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
            ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
            qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
            pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
            rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
            goto 30
          endif
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
           dz1=height(iz)-uvwzlev(ix,jy,kz-1)
           dz2=uvwzlev(ix,jy,kz)-height(iz)
           dz=dz1+dz2
           uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+ &
           uuhn(ix,jy,kz,l)*dz1)/dz
           vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+ &
           vvhn(ix,jy,kz,l)*dz1)/dz
           ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2+ &
           tthn(ix,jy,kz,n,l)*dz1)/dz
           qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2+ &
           qvhn(ix,jy,kz,n,l)*dz1)/dz
           pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+ &
           pvhn(ix,jy,kz,l)*dz1)/dz
           rhon(ix,jy,iz,n,l)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
           kmin=kz
           goto 30
          endif
        end do
30      continue
      end do


  ! Levels, where w is given
  !*************************

      wwn(ix,jy,1,n,l)=wwhn(ix,jy,1,l)*pinmconv(1)
      wwn(ix,jy,nz,n,l)=wwhn(ix,jy,nwz,l)*pinmconv(nz)
      kmin=2
      do iz=2,nz
        do kz=kmin,nwz
          if ((height(iz).gt.wzlev(kz-1)).and. &
          (height(iz).le.wzlev(kz))) then
            dz1=height(iz)-wzlev(kz-1)
            dz2=wzlev(kz)-height(iz)
            dz=dz1+dz2
            wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(kz-1)*dz2 &
            +wwhn(ix,jy,kz,l)*pinmconv(kz)*dz1)/dz
            kmin=kz
            goto 40
          endif
        end do
40      continue
      end do

  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodzn(ix,jy,1,n,l)=(rhon(ix,jy,2,n,l)-rhon(ix,jy,1,n,l))/ &
      (height(2)-height(1))
      do kz=2,nz-1
        drhodzn(ix,jy,kz,n,l)=(rhon(ix,jy,kz+1,n,l)- &
        rhon(ix,jy,kz-1,n,l))/(height(kz+1)-height(kz-1))
      end do
      drhodzn(ix,jy,nz,n,l)=drhodzn(ix,jy,nz-1,n,l)

    end do
  end do


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,nyn(l)-2
    cosf=cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
    do ix=1,nxn(l)-2

      kmin=2
      do iz=2,nz-1

        ui=uun(ix,jy,iz,n,l)*dxconst*xresoln(l)/cosf
        vi=vvn(ix,jy,iz,n,l)*dyconst*yresoln(l)

        do kz=kmin,nz
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
               (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            kmin=kz
            goto 47
          endif
        end do

47      ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*ui+dzdy*vi)

      end do

    end do
  end do


  !write (*,*) 'initializing nested cloudsn, n:',n
  !   create a cloud and rainout/washout field, cloudsn occur where rh>80%
  do jy=0,nyn(l)-1
    do ix=0,nxn(l)-1
      rain_cloud_above=0
      lsp=lsprecn(ix,jy,1,n,l)
      convp=convprecn(ix,jy,1,n,l)
      cloudsnh(ix,jy,n,l)=0
      do kz_inv=1,nz-1
         kz=nz-kz_inv+1
         pressure=rhon(ix,jy,kz,n,l)*r_air*ttn(ix,jy,kz,n,l)
         rh=qvn(ix,jy,kz,n,l)/f_qvsat(pressure,ttn(ix,jy,kz,n,l))
         cloudsn(ix,jy,kz,n,l)=0
         if (rh.gt.0.8) then ! in cloud
            if ((lsp.gt.0.01).or.(convp.gt.0.01)) then
               rain_cloud_above=1
               cloudsnh(ix,jy,n,l)=cloudsnh(ix,jy,n,l)+ &
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
            if (rain_cloud_above.eq.1) then ! scavenging
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

  end do

end subroutine verttransform_nests
