! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine interpol_wind(itime,xt,yt,zt)
  !                           i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
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
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: i,m,n,indexh,indzh
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

  ! Determine the level below the current position for u,v
  !*******************************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      goto 6
    endif
  end do
6   continue


  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz

  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************

  usl=0.
  vsl=0.
  wsl=0.
  usq=0.
  vsq=0.
  wsq=0.
  do m=1,2
    indexh=memind(m)
    do n=1,2
      indzh=indz+n-1

      if (ngrid.lt.0) then
        u1(n)=p1*uupol(ix ,jy ,indzh,indexh) &
             +p2*uupol(ixp,jy ,indzh,indexh) &
             +p3*uupol(ix ,jyp,indzh,indexh) &
             +p4*uupol(ixp,jyp,indzh,indexh)
        v1(n)=p1*vvpol(ix ,jy ,indzh,indexh) &
             +p2*vvpol(ixp,jy ,indzh,indexh) &
             +p3*vvpol(ix ,jyp,indzh,indexh) &
             +p4*vvpol(ixp,jyp,indzh,indexh)
        usl=usl+uupol(ix ,jy ,indzh,indexh)+ &
             uupol(ixp,jy ,indzh,indexh) &
             +uupol(ix ,jyp,indzh,indexh)+uupol(ixp,jyp,indzh,indexh)
        vsl=vsl+vvpol(ix ,jy ,indzh,indexh)+ &
             vvpol(ixp,jy ,indzh,indexh) &
             +vvpol(ix ,jyp,indzh,indexh)+vvpol(ixp,jyp,indzh,indexh)

        usq=usq+uupol(ix ,jy ,indzh,indexh)* &
             uupol(ix ,jy ,indzh,indexh)+ &
             uupol(ixp,jy ,indzh,indexh)*uupol(ixp,jy ,indzh,indexh)+ &
             uupol(ix ,jyp,indzh,indexh)*uupol(ix ,jyp,indzh,indexh)+ &
             uupol(ixp,jyp,indzh,indexh)*uupol(ixp,jyp,indzh,indexh)
        vsq=vsq+vvpol(ix ,jy ,indzh,indexh)* &
             vvpol(ix ,jy ,indzh,indexh)+ &
             vvpol(ixp,jy ,indzh,indexh)*vvpol(ixp,jy ,indzh,indexh)+ &
             vvpol(ix ,jyp,indzh,indexh)*vvpol(ix ,jyp,indzh,indexh)+ &
             vvpol(ixp,jyp,indzh,indexh)*vvpol(ixp,jyp,indzh,indexh)
      else
        u1(n)=p1*uu(ix ,jy ,indzh,indexh) &
             +p2*uu(ixp,jy ,indzh,indexh) &
             +p3*uu(ix ,jyp,indzh,indexh) &
             +p4*uu(ixp,jyp,indzh,indexh)
        v1(n)=p1*vv(ix ,jy ,indzh,indexh) &
             +p2*vv(ixp,jy ,indzh,indexh) &
             +p3*vv(ix ,jyp,indzh,indexh) &
             +p4*vv(ixp,jyp,indzh,indexh)
        usl=usl+uu(ix ,jy ,indzh,indexh)+uu(ixp,jy ,indzh,indexh) &
             +uu(ix ,jyp,indzh,indexh)+uu(ixp,jyp,indzh,indexh)
        vsl=vsl+vv(ix ,jy ,indzh,indexh)+vv(ixp,jy ,indzh,indexh) &
             +vv(ix ,jyp,indzh,indexh)+vv(ixp,jyp,indzh,indexh)

        usq=usq+uu(ix ,jy ,indzh,indexh)*uu(ix ,jy ,indzh,indexh)+ &
             uu(ixp,jy ,indzh,indexh)*uu(ixp,jy ,indzh,indexh)+ &
             uu(ix ,jyp,indzh,indexh)*uu(ix ,jyp,indzh,indexh)+ &
             uu(ixp,jyp,indzh,indexh)*uu(ixp,jyp,indzh,indexh)
        vsq=vsq+vv(ix ,jy ,indzh,indexh)*vv(ix ,jy ,indzh,indexh)+ &
             vv(ixp,jy ,indzh,indexh)*vv(ixp,jy ,indzh,indexh)+ &
             vv(ix ,jyp,indzh,indexh)*vv(ix ,jyp,indzh,indexh)+ &
             vv(ixp,jyp,indzh,indexh)*vv(ixp,jyp,indzh,indexh)
      endif
      w1(n)=p1*ww(ix ,jy ,indzh,indexh) &
           +p2*ww(ixp,jy ,indzh,indexh) &
           +p3*ww(ix ,jyp,indzh,indexh) &
           +p4*ww(ixp,jyp,indzh,indexh)
      wsl=wsl+ww(ix ,jy ,indzh,indexh)+ww(ixp,jy ,indzh,indexh) &
           +ww(ix ,jyp,indzh,indexh)+ww(ixp,jyp,indzh,indexh)
      wsq=wsq+ww(ix ,jy ,indzh,indexh)*ww(ix ,jy ,indzh,indexh)+ &
           ww(ixp,jy ,indzh,indexh)*ww(ixp,jy ,indzh,indexh)+ &
           ww(ix ,jyp,indzh,indexh)*ww(ix ,jyp,indzh,indexh)+ &
           ww(ixp,jyp,indzh,indexh)*ww(ixp,jyp,indzh,indexh)
    end do


  !**********************************
  ! 2.) Linear vertical interpolation
  !**********************************

    uh(m)=dz2*u1(1)+dz1*u1(2)
    vh(m)=dz2*v1(1)+dz1*v1(2)
    wh(m)=dz2*w1(1)+dz1*w1(2)
  end do


  !************************************
  ! 3.) Temporal interpolation (linear)
  !************************************

  u=(uh(1)*dt2+uh(2)*dt1)*dtt
  v=(vh(1)*dt2+vh(2)*dt1)*dtt
  w=(wh(1)*dt2+wh(2)*dt1)*dtt


  ! Compute standard deviations
  !****************************

  xaux=usq-usl*usl/16.
  if (xaux.lt.eps) then
    usig=0.
  else
    usig=sqrt(xaux/15.)
  endif

  xaux=vsq-vsl*vsl/16.
  if (xaux.lt.eps) then
    vsig=0.
  else
    vsig=sqrt(xaux/15.)
  endif


  xaux=wsq-wsl*wsl/16.
  if (xaux.lt.eps) then
    wsig=0.
  else
    wsig=sqrt(xaux/15.)
  endif

end subroutine interpol_wind
