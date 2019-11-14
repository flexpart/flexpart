! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine interpol_misslev(n)
  !                            i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates u,v,w, density and density gradients.        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !    Update: 2 March 1999                                                    *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  level                                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod

  implicit none

  ! Auxiliary variables needed for interpolation
  real :: y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
  real :: usl,vsl,wsl,usq,vsq,wsq,xaux
  integer :: m,n,indexh
  real,parameter :: eps=1.0e-30


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************


  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************

  ! Loop over 2 time steps
  !***********************

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


end subroutine interpol_misslev
