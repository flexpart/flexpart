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

subroutine interpol_all_nests(itime,xt,yt,zt)
  !                                i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates everything that is needed for calculating the*
  !  dispersion.                                                               *
  !  Version for interpolating nested grids.                                   *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    9 February 1999                                                         *
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
  !                    calculated                                              *
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

    ust1(m)=p1*ustarn(ix ,jy ,1,indexh,ngrid) &
         + p2*ustarn(ixp,jy ,1,indexh,ngrid) &
         + p3*ustarn(ix ,jyp,1,indexh,ngrid) &
         + p4*ustarn(ixp,jyp,1,indexh,ngrid)
    wst1(m)=p1*wstarn(ix ,jy ,1,indexh,ngrid) &
         + p2*wstarn(ixp,jy ,1,indexh,ngrid) &
         + p3*wstarn(ix ,jyp,1,indexh,ngrid) &
         + p4*wstarn(ixp,jyp,1,indexh,ngrid)
    oli1(m)=p1*olin(ix ,jy ,1,indexh,ngrid) &
         + p2*olin(ixp,jy ,1,indexh,ngrid) &
         + p3*olin(ix ,jyp,1,indexh,ngrid) &
         + p4*olin(ixp,jyp,1,indexh,ngrid)
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

  do n=indz,indz+1
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
      rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid) &
           +p2*drhodzn(ixp,jy ,n,indexh,ngrid) &
           +p3*drhodzn(ix ,jyp,n,indexh,ngrid) &
           +p4*drhodzn(ixp,jyp,n,indexh,ngrid)
      rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid) &
           +p2*rhon(ixp,jy ,n,indexh,ngrid) &
           +p3*rhon(ix ,jyp,n,indexh,ngrid) &
           +p4*rhon(ixp,jyp,n,indexh,ngrid)

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

  end do

end subroutine interpol_all_nests
