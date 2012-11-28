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

subroutine partoutput_short(itime)
  !                              i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss,numshortout,numshortall
  integer :: ix,jy,ixp,jyp
  real :: xlon,ylat,zlim,dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,topo
  character :: adate*8,atime*6

  integer(kind=2) :: idump(3,maxpart)
  integer :: i4dump(maxpart)


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! Some variables needed for temporal interpolation
  !*************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  ! Loop about all particles
  !*************************

  numshortout=0
  numshortall=0
  do i=1,numpart

  ! Take only valid particles
  !**************************

    if (itra1(i).eq.itime) then
      xlon=xlon0+xtra1(i)*dx
      ylat=ylat0+ytra1(i)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

      ix=xtra1(i)
      jy=ytra1(i)
      ixp=ix+1
      jyp=jy+1
      ddx=xtra1(i)-real(ix)
      ddy=ytra1(i)-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

  ! Topography
  !***********

      topo=p1*oro(ix ,jy) &
           + p2*oro(ixp,jy) &
           + p3*oro(ix ,jyp) &
           + p4*oro(ixp,jyp)


  ! Convert positions to integer*2 variables (from -32768 to 32767)
  ! Do this only for region of main interest, i.e. extended North Atlantic region,
  ! and for the tracer of interest, i.e. the North American one
  !*****************************************************************************

      if (xlon.gt.180.) xlon=xlon-360.
      if (xlon.lt.-180.) xlon=xlon+360.

      numshortall=numshortall+1
      if ((xlon.gt.-140).and.(xlon.lt.60).and.(ylat.gt.10).and. &
           (xmass1(i,1).gt.0.)) then
        numshortout=numshortout+1
        idump(1,numshortout)=nint(xlon*180.)
        idump(2,numshortout)=nint(ylat*360.)
        zlim=min(ztra1(i)+topo,32766.)
        idump(3,numshortout)=nint(zlim)
        i4dump(numshortout)=npoint(i)
      endif

    endif
  end do


  ! Open output file and write the output
  !**************************************

  open(unitshortpart,file=path(2)(1:length(2))//'shortposit_'//adate// &
       atime,form='unformatted')

  ! Write current time to file
  !***************************

  write(unitshortpart) itime
  write(unitshortpart) numshortout
  write(unitshortpart) &
       (i4dump(i),(idump(j,i),j=1,3),i=1,numshortout)


  write(*,*) numshortout,numshortall

  close(unitshortpart)

end subroutine partoutput_short
