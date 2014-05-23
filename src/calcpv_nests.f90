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

subroutine calcpv_nests(l,n,uuhn,vvhn,pvhn)
  !                     i i  i    i    o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of potential vorticity on 3-d nested grid                     *
  !                                                                            *
  !     Author: P. James                                                       *
  !     22 February 2000                                                       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 2)       *
  ! l                  index of current nest                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: n,ix,jy,i,j,k,kl,ii,jj,klvrp,klvrm,klpt,kup,kdn,kch
  integer :: jyvp,jyvm,ixvp,ixvm,jumpx,jumpy,jux,juy,ivrm,ivrp,ivr
  integer :: nlck,l
  real :: vx(2),uy(2),phi,tanphi,cosphi,dvdx,dudy,f
  real :: theta,thetap,thetam,dthetadp,dt1,dt2,dt
  real :: ppml(0:nxmaxn-1,0:nymaxn-1,nuvzmax),ppmk(0:nxmaxn-1,0:nymaxn-1,nuvzmax)
  real :: thup,thdn
  real,parameter :: eps=1.e-5,p0=101325
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)

  ! Set number of levels to check for adjacent theta
  nlck=nuvz/3
  !
  ! Loop over entire grid
  !**********************
  do kl=1,nuvz
    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
         ppml(ix,jy,kl)=akz(kl)+bkz(kl)*psn(ix,jy,1,n,l)
      enddo
    enddo
  enddo
  ppmk=(100000./ppml)**kappa

  do jy=0,nyn(l)-1
    phi = (ylat0n(l) + jy * dyn(l)) * pi / 180.
    f = 0.00014585 * sin(phi)
    tanphi = tan(phi)
    cosphi = cos(phi)
  ! Provide a virtual jy+1 and jy-1 in case we are on domain edge (Lat)
      jyvp=jy+1
      jyvm=jy-1
      if (jy.eq.0) jyvm=0
      if (jy.eq.nyn(l)-1) jyvp=nyn(l)-1
  ! Define absolute gap length
      jumpy=2
      if (jy.eq.0.or.jy.eq.nyn(l)-1) jumpy=1
      juy=jumpy
  !
    do ix=0,nxn(l)-1
  ! Provide a virtual ix+1 and ix-1 in case we are on domain edge (Long)
      ixvp=ix+1
      ixvm=ix-1
      jumpx=2
      if (ix.eq.0) ixvm=0
      if (ix.eq.nxn(l)-1) ixvp=nxn(l)-1
      ivrp=ixvp
      ivrm=ixvm
  ! Define absolute gap length
      if (ix.eq.0.or.ix.eq.nxn(l)-1) jumpx=1
      jux=jumpx
  !
  ! Loop over the vertical
  !***********************

      do kl=1,nuvz
        theta=tthn(ix,jy,kl,n,l)*ppmk(ix,jy,kl)
        klvrp=kl+1
        klvrm=kl-1
        klpt=kl
  ! If top or bottom level, dthetadp is evaluated between the current
  ! level and the level inside, otherwise between level+1 and level-1
  !
        if (klvrp.gt.nuvz) klvrp=nuvz
        if (klvrm.lt.1) klvrm=1
        thetap=tthn(ix,jy,klvrp,n,l)*ppmk(ix,jy,klvrp)
        thetam=tthn(ix,jy,klvrm,n,l)*ppmk(ix,jy,klvrm)
        dthetadp=(thetap-thetam)/(ppml(ix,jy,klvrp)-ppml(ix,jy,klvrm))

  ! Compute vertical position at pot. temperature surface on subgrid
  ! and the wind at that position
  !*****************************************************************
  ! a) in x direction
        ii=0
        do i=ixvm,ixvp,jumpx
          ivr=i
          ii=ii+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
40        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 21     ! No more levels to check,
  !                                       ! and no values found
          if (kup.ge.nuvz) goto 41
          kch=kch+1
          k=kup
          thdn=tthn(ivr,jy,k,n,l)*ppmk(ivr,jy,k)
          thup=tthn(ivr,jy,k+1,n,l)*ppmk(ivr,jy,k+1)

      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
    vx(ii)=(vvhn(ivr,jy,k,l)*dt2+vvhn(ivr,jy,k+1,l)*dt1)/dt
              goto 20
            endif
41        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 40
          kch=kch+1
          k=kdn
          thdn=tthn(ivr,jy,k,n,l)*ppmk(ivr,jy,k)
          thup=tthn(ivr,jy,k+1,n,l)*ppmk(ivr,jy,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
    vx(ii)=(vvhn(ivr,jy,k,l)*dt2+vvhn(ivr,jy,k+1,l)*dt1)/dt
              goto 20
            endif
            goto 40
  ! This section used when no values were found
21      continue
  ! Must use vv at current level and long. jux becomes smaller by 1
        vx(ii)=vvhn(ix,jy,kl,l)
        jux=jux-1
  ! Otherwise OK
20        continue
        end do
      if (jux.gt.0) then
      dvdx=(vx(2)-vx(1))/real(jux)/(dxn(l)*pi/180.)
      else
      dvdx=vvhn(ivrp,jy,kl,l)-vvhn(ivrm,jy,kl,l)
      dvdx=dvdx/real(jumpx)/(dxn(l)*pi/180.)
  ! Only happens if no equivalent theta value
  ! can be found on either side, hence must use values
  ! from either side, same pressure level.
      end if

  ! b) in y direction

        jj=0
        do j=jyvm,jyvp,jumpy
          jj=jj+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
70        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 51     ! No more levels to check,
  !                                     ! and no values found
          if (kup.ge.nuvz) goto 71
          kch=kch+1
          k=kup
          thdn=tthn(ix,j,k,n,l)*ppmk(ix,j,k)
          thup=tthn(ix,j,k+1,n,l)*ppmk(ix,j,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
        uy(jj)=(uuhn(ix,j,k,l)*dt2+uuhn(ix,j,k+1,l)*dt1)/dt
              goto 50
            endif
71        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 70
          kch=kch+1
          k=kdn
          thdn=tthn(ix,j,k,n,l)*ppmk(ix,j,k)
          thup=tthn(ix,j,k+1,n,l)*ppmk(ix,j,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
        uy(jj)=(uuhn(ix,j,k,l)*dt2+uuhn(ix,j,k+1,l)*dt1)/dt
              goto 50
            endif
            goto 70
  ! This section used when no values were found
51      continue
  ! Must use uu at current level and lat. juy becomes smaller by 1
        uy(jj)=uuhn(ix,jy,kl,l)
        juy=juy-1
  ! Otherwise OK
50        continue
        end do
      if (juy.gt.0) then
      dudy=(uy(2)-uy(1))/real(juy)/(dyn(l)*pi/180.)
      else
      dudy=uuhn(ix,jyvp,kl,l)-uuhn(ix,jyvm,kl,l)
      dudy=dudy/real(jumpy)/(dyn(l)*pi/180.)
      end if
  !
      pvhn(ix,jy,kl,l)=dthetadp*(f+(dvdx/cosphi-dudy &
           +uuhn(ix,jy,kl,l)*tanphi)/r_earth)*(-1.e6)*9.81

  !
  ! Resest jux and juy
      jux=jumpx
      juy=jumpy
      end do
    end do
  end do
  !
end subroutine calcpv_nests
