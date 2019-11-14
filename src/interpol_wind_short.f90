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

subroutine interpol_wind_short(itime,xt,yt,zt)
!                                 i   i  i  i
!*****************************************************************************
!                                                                            *
!  This subroutine interpolates the wind data to current trajectory position.*
!                                                                            *
!    Author: A. Stohl                                                        *
!                                                                            *
!    16 December 1997                                                        *
!                                                                            *
!  Revision March 2005 by AST : all output variables in common block         *
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

  integer, intent(in) :: itime
  real, intent(in) :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
  integer :: i,m,n,indexh,indzh


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

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
      else
        u1(n)=p1*uu(ix ,jy ,indzh,indexh) &
             +p2*uu(ixp,jy ,indzh,indexh) &
             +p3*uu(ix ,jyp,indzh,indexh) &
             +p4*uu(ixp,jyp,indzh,indexh)
        v1(n)=p1*vv(ix ,jy ,indzh,indexh) &
             +p2*vv(ixp,jy ,indzh,indexh) &
             +p3*vv(ix ,jyp,indzh,indexh) &
             +p4*vv(ixp,jyp,indzh,indexh)
      endif
      w1(n)=p1*ww(ix ,jy ,indzh,indexh) &
           +p2*ww(ixp,jy ,indzh,indexh) &
           +p3*ww(ix ,jyp,indzh,indexh) &
           +p4*ww(ixp,jyp,indzh,indexh)
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

end subroutine interpol_wind_short
