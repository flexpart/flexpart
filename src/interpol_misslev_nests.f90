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

subroutine interpol_misslev_nests(n)
  !                                  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates u,v,w, density and density gradients.        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
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
    y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid) &
         +p2*uun(ixp,jy ,n,indexh,ngrid) &
         +p3*uun(ix ,jyp,n,indexh,ngrid) &
         +p4*uun(ixp,jyp,n,indexh,ngrid)
    y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid) &
         +p2*vvn(ixp,jy ,n,indexh,ngrid) &
         +p3*vvn(ix ,jyp,n,indexh,ngrid) &
         +p4*vvn(ixp,jyp,n,indexh,ngrid)
    y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid) &
         +p2*wwn(ixp,jy ,n,indexh,ngrid) &
         +p3*wwn(ix ,jyp,n,indexh,ngrid) &
         +p4*wwn(ixp,jyp,n,indexh,ngrid)
    rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
         +p2*rhon(ixp,jy ,n,indexh,ngrid) &
         +p3*rhon(ix ,jyp,n,indexh,ngrid) &
         +p4*rhon(ixp,jyp,n,indexh,ngrid)
    rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
         +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
         +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
         +p4*drhodzn(ixp,jyp,n,indexh,ngrid)

     usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid) &
          +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
     vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid) &
          +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
     wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid) &
          +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

    usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+ &
         uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+ &
         uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+ &
         uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
    vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+ &
         vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+ &
         vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+ &
         vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
    wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+ &
         wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+ &
         wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+ &
         wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
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

end subroutine interpol_misslev_nests
