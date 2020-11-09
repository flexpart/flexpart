! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine interpol_wind_nests(itime,xt,yt,zt)
  !                                 i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
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

      u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
           +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
           +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
           +p4*uun(ixp,jyp,indzh,indexh,ngrid)
      v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
      w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

      usl=usl+uun(ix ,jy ,indzh,indexh,ngrid)+ &
           uun(ixp,jy ,indzh,indexh,ngrid) &
           +uun(ix ,jyp,indzh,indexh,ngrid)+ &
           uun(ixp,jyp,indzh,indexh,ngrid)
      vsl=vsl+vvn(ix ,jy ,indzh,indexh,ngrid)+ &
           vvn(ixp,jy ,indzh,indexh,ngrid) &
           +vvn(ix ,jyp,indzh,indexh,ngrid)+ &
           vvn(ixp,jyp,indzh,indexh,ngrid)
      wsl=wsl+wwn(ix ,jy ,indzh,indexh,ngrid)+ &
           wwn(ixp,jy ,indzh,indexh,ngrid) &
           +wwn(ix ,jyp,indzh,indexh,ngrid)+ &
           wwn(ixp,jyp,indzh,indexh,ngrid)

      usq=usq+uun(ix ,jy ,indzh,indexh,ngrid)* &
           uun(ix ,jy ,indzh,indexh,ngrid)+ &
           uun(ixp,jy ,indzh,indexh,ngrid)*uun(ixp,jy ,indzh,indexh,ngrid)+ &
           uun(ix ,jyp,indzh,indexh,ngrid)*uun(ix ,jyp,indzh,indexh,ngrid)+ &
           uun(ixp,jyp,indzh,indexh,ngrid)*uun(ixp,jyp,indzh,indexh,ngrid)
      vsq=vsq+vvn(ix ,jy ,indzh,indexh,ngrid)* &
           vvn(ix ,jy ,indzh,indexh,ngrid)+ &
           vvn(ixp,jy ,indzh,indexh,ngrid)*vvn(ixp,jy ,indzh,indexh,ngrid)+ &
           vvn(ix ,jyp,indzh,indexh,ngrid)*vvn(ix ,jyp,indzh,indexh,ngrid)+ &
           vvn(ixp,jyp,indzh,indexh,ngrid)*vvn(ixp,jyp,indzh,indexh,ngrid)
      wsq=wsq+wwn(ix ,jy ,indzh,indexh,ngrid)* &
           wwn(ix ,jy ,indzh,indexh,ngrid)+ &
           wwn(ixp,jy ,indzh,indexh,ngrid)*wwn(ixp,jy ,indzh,indexh,ngrid)+ &
           wwn(ix ,jyp,indzh,indexh,ngrid)*wwn(ix ,jyp,indzh,indexh,ngrid)+ &
           wwn(ixp,jyp,indzh,indexh,ngrid)*wwn(ixp,jyp,indzh,indexh,ngrid)
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

end subroutine interpol_wind_nests
