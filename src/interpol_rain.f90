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

subroutine interpol_rain(yy1,yy2,yy3,nxmax,nymax,nzmax,nx, &
       ny,memind,xt,yt,level,itime1,itime2,itime,yint1,yint2,yint3)
  !                          i   i   i    i    i     i   i
  !i    i    i  i    i     i      i      i     o     o     o
  !****************************************************************************
  !                                                                           *
  !  Interpolation of meteorological fields on 2-d model layers.              *
  !  In horizontal direction bilinear interpolation interpolation is used.    *
  !  Temporally a linear interpolation is used.                               *
  !  Three fields are interpolated at the same time.                          *
  !                                                                           *
  !  This is a special version of levlininterpol to save CPU time.            *
  !                                                                           *
  !  1 first time                                                             *
  !  2 second time                                                            *
  !                                                                           *
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     30 August 1996                                                        *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  !                                                                           *
  ! dt1,dt2              time differences between fields and current position *
  ! dz1,dz2              z distance between levels and current position       *
  ! height(nzmax)        heights of the model levels                          *
  ! indexh               help variable                                        *
  ! indz                 the level closest to the current trajectory position *
  ! indzh                help variable                                        *
  ! itime                current time                                         *
  ! itime1               time of the first wind field                         *
  ! itime2               time of the second wind field                        *
  ! ix,jy                x,y coordinates of lower left subgrid point          *
  ! level                level at which interpolation shall be done           *
  ! memind(3)            points to the places of the wind fields              *
  ! nx,ny                actual field dimensions in x,y and z direction       *
  ! nxmax,nymax,nzmax    maximum field dimensions in x,y and z direction      *
  ! xt                   current x coordinate                                 *
  ! yint                 the final interpolated value                         *
  ! yt                   current y coordinate                                 *
  ! yy(0:nxmax,0:nymax,nzmax,3) meteorological field used for interpolation   *
  ! zt                   current z coordinate                                 *
  !                                                                           *
  !****************************************************************************

  implicit none

  integer :: nx,ny,nxmax,nymax,nzmax,memind(2),m,ix,jy,ixp,jyp
  integer :: itime,itime1,itime2,level,indexh
  real :: yy1(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: yy2(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: yy3(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: ddx,ddy,rddx,rddy,dt1,dt2,dt,y1(2),y2(2),y3(2)
  real :: xt,yt,yint1,yint2,yint3,p1,p2,p3,p4



  ! If point at border of grid -> small displacement into grid
  !***********************************************************

  if (xt.ge.real(nx-1)) xt=real(nx-1)-0.00001
  if (yt.ge.real(ny-1)) yt=real(ny-1)-0.00001



  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 2 fields (Temporal)
  !*******************************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ix=int(xt)
  jy=int(yt)
  ixp=ix+1
  jyp=jy+1
  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy


  ! Loop over 2 time steps
  !***********************

  do m=1,2
    indexh=memind(m)

    y1(m)=p1*yy1(ix ,jy ,level,indexh) &
         + p2*yy1(ixp,jy ,level,indexh) &
         + p3*yy1(ix ,jyp,level,indexh) &
         + p4*yy1(ixp,jyp,level,indexh)
    y2(m)=p1*yy2(ix ,jy ,level,indexh) &
         + p2*yy2(ixp,jy ,level,indexh) &
         + p3*yy2(ix ,jyp,level,indexh) &
         + p4*yy2(ixp,jyp,level,indexh)
    y3(m)=p1*yy3(ix ,jy ,level,indexh) &
         + p2*yy3(ixp,jy ,level,indexh) &
         + p3*yy3(ix ,jyp,level,indexh) &
         + p4*yy3(ixp,jyp,level,indexh)
  end do


  !************************************
  ! 2.) Temporal interpolation (linear)
  !************************************

  dt1=real(itime-itime1)
  dt2=real(itime2-itime)
  dt=dt1+dt2

  yint1=(y1(1)*dt2+y1(2)*dt1)/dt
  yint2=(y2(1)*dt2+y2(2)*dt1)/dt
  yint3=(y3(1)*dt2+y3(2)*dt1)/dt


end subroutine interpol_rain
