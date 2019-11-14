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

subroutine interpol_all(itime,xt,yt,zt)
  !                          i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates everything that is needed for calculating the*
  !  dispersion.                                                               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    culated                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: i,m,n,indexh
  real,parameter :: eps=1.0e-30


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

  ! Calculate variables for time interpolation
  !*******************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  !*****************************************
  ! 1. Interpolate u*, w* and Obukhov length
  !*****************************************

  ! a) Bilinear horizontal interpolation

  do m=1,2
    indexh=memind(m)

    ust1(m)=p1*ustar(ix ,jy ,1,indexh) &
         + p2*ustar(ixp,jy ,1,indexh) &
         + p3*ustar(ix ,jyp,1,indexh) &
         + p4*ustar(ixp,jyp,1,indexh)
    wst1(m)=p1*wstar(ix ,jy ,1,indexh) &
         + p2*wstar(ixp,jy ,1,indexh) &
         + p3*wstar(ix ,jyp,1,indexh) &
         + p4*wstar(ixp,jyp,1,indexh)
    oli1(m)=p1*oli(ix ,jy ,1,indexh) &
         + p2*oli(ixp,jy ,1,indexh) &
         + p3*oli(ix ,jyp,1,indexh) &
         + p4*oli(ixp,jyp,1,indexh)
  end do

  ! b) Temporal interpolation

  ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
  wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
  oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

  if (oliaux.ne.0.) then
    ol=1./oliaux
  else
    ol=99999.
  endif


  !*****************************************************
  ! 2. Interpolate vertical profiles of u,v,w,rho,drhodz
  !*****************************************************


  ! Determine the level below the current position
  !***********************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      indzp=i
      goto 6
    endif
  end do
6   continue

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps and indz levels
  !***************************************

  do n=indz,indzp
    usl=0.
    vsl=0.
    wsl=0.
    usq=0.
    vsq=0.
    wsq=0.
    do m=1,2
      indexh=memind(m)
      if (ngrid.lt.0) then
        y1(m)=p1*uupol(ix ,jy ,n,indexh) &
             +p2*uupol(ixp,jy ,n,indexh) &
             +p3*uupol(ix ,jyp,n,indexh) &
             +p4*uupol(ixp,jyp,n,indexh)
        y2(m)=p1*vvpol(ix ,jy ,n,indexh) &
             +p2*vvpol(ixp,jy ,n,indexh) &
             +p3*vvpol(ix ,jyp,n,indexh) &
             +p4*vvpol(ixp,jyp,n,indexh)
        usl=usl+uupol(ix ,jy ,n,indexh)+uupol(ixp,jy ,n,indexh) &
             +uupol(ix ,jyp,n,indexh)+uupol(ixp,jyp,n,indexh)
        vsl=vsl+vvpol(ix ,jy ,n,indexh)+vvpol(ixp,jy ,n,indexh) &
             +vvpol(ix ,jyp,n,indexh)+vvpol(ixp,jyp,n,indexh)

        usq=usq+uupol(ix ,jy ,n,indexh)*uupol(ix ,jy ,n,indexh)+ &
             uupol(ixp,jy ,n,indexh)*uupol(ixp,jy ,n,indexh)+ &
             uupol(ix ,jyp,n,indexh)*uupol(ix ,jyp,n,indexh)+ &
             uupol(ixp,jyp,n,indexh)*uupol(ixp,jyp,n,indexh)
        vsq=vsq+vvpol(ix ,jy ,n,indexh)*vvpol(ix ,jy ,n,indexh)+ &
             vvpol(ixp,jy ,n,indexh)*vvpol(ixp,jy ,n,indexh)+ &
             vvpol(ix ,jyp,n,indexh)*vvpol(ix ,jyp,n,indexh)+ &
             vvpol(ixp,jyp,n,indexh)*vvpol(ixp,jyp,n,indexh)
      else
        y1(m)=p1*uu(ix ,jy ,n,indexh) &
             +p2*uu(ixp,jy ,n,indexh) &
             +p3*uu(ix ,jyp,n,indexh) &
             +p4*uu(ixp,jyp,n,indexh)
        y2(m)=p1*vv(ix ,jy ,n,indexh) &
             +p2*vv(ixp,jy ,n,indexh) &
             +p3*vv(ix ,jyp,n,indexh) &
             +p4*vv(ixp,jyp,n,indexh)
        usl=usl+uu(ix ,jy ,n,indexh)+uu(ixp,jy ,n,indexh) &
             +uu(ix ,jyp,n,indexh)+uu(ixp,jyp,n,indexh)
        vsl=vsl+vv(ix ,jy ,n,indexh)+vv(ixp,jy ,n,indexh) &
             +vv(ix ,jyp,n,indexh)+vv(ixp,jyp,n,indexh)

        usq=usq+uu(ix ,jy ,n,indexh)*uu(ix ,jy ,n,indexh)+ &
             uu(ixp,jy ,n,indexh)*uu(ixp,jy ,n,indexh)+ &
             uu(ix ,jyp,n,indexh)*uu(ix ,jyp,n,indexh)+ &
             uu(ixp,jyp,n,indexh)*uu(ixp,jyp,n,indexh)
        vsq=vsq+vv(ix ,jy ,n,indexh)*vv(ix ,jy ,n,indexh)+ &
             vv(ixp,jy ,n,indexh)*vv(ixp,jy ,n,indexh)+ &
             vv(ix ,jyp,n,indexh)*vv(ix ,jyp,n,indexh)+ &
             vv(ixp,jyp,n,indexh)*vv(ixp,jyp,n,indexh)
      endif
      y3(m)=p1*ww(ix ,jy ,n,indexh) &
           +p2*ww(ixp,jy ,n,indexh) &
           +p3*ww(ix ,jyp,n,indexh) &
           +p4*ww(ixp,jyp,n,indexh)
      rhograd1(m)=p1*drhodz(ix ,jy ,n,indexh) &
           +p2*drhodz(ixp,jy ,n,indexh) &
           +p3*drhodz(ix ,jyp,n,indexh) &
           +p4*drhodz(ixp,jyp,n,indexh)
      rho1(m)=p1*rho(ix ,jy ,n,indexh) &
           +p2*rho(ixp,jy ,n,indexh) &
           +p3*rho(ix ,jyp,n,indexh) &
           +p4*rho(ixp,jyp,n,indexh)
      wsl=wsl+ww(ix ,jy ,n,indexh)+ww(ixp,jy ,n,indexh) &
           +ww(ix ,jyp,n,indexh)+ww(ixp,jyp,n,indexh)
      wsq=wsq+ww(ix ,jy ,n,indexh)*ww(ix ,jy ,n,indexh)+ &
           ww(ixp,jy ,n,indexh)*ww(ixp,jy ,n,indexh)+ &
           ww(ix ,jyp,n,indexh)*ww(ix ,jyp,n,indexh)+ &
           ww(ixp,jyp,n,indexh)*ww(ixp,jyp,n,indexh)
    end do
    uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
    vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
    wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
    rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
    rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
    indzindicator(n)=.false.

  ! Compute standard deviations
  !****************************

    xaux=usq-usl*usl/8.
    if (xaux.lt.eps) then
      usigprof(n)=0.
    else
      usigprof(n)=sqrt(xaux/7.)
    endif

    xaux=vsq-vsl*vsl/8.
    if (xaux.lt.eps) then
      vsigprof(n)=0.
    else
      vsigprof(n)=sqrt(xaux/7.)
    endif


    xaux=wsq-wsl*wsl/8.
    if (xaux.lt.eps) then
      wsigprof(n)=0.
    else
      wsigprof(n)=sqrt(xaux/7.)
    endif

  end do


end subroutine interpol_all
