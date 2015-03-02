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

subroutine partoutput(itime)
  !                           i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !     12/2014 eso: Version for MPI                                           *
  !                  processes sequentially access and append data to file     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use mpi_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: xlon,ylat
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,hm(2),hmixi,pv1(2),pvprof(2),pvi,qv1(2),qvprof(2),qvi
  real :: tt1(2),ttprof(2),tti,rho1(2),rhoprof(2),rhoi
  real :: tr(2),tri
  character :: adate*8,atime*6
  character(LEN=8) :: file_stat='OLD'

  ! MPI root process creates the file, other processes append to it
  if (lroot) file_stat='REPLACE'

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

  ! Open output file and write the output
  !**************************************

  if (ipout.eq.1) then
    open(unitpartout,file=path(2)(1:length(2))//'partposit_'//adate// &
         atime,form='unformatted',status=file_stat,position='append')
  else
    open(unitpartout,file=path(2)(1:length(2))//'partposit_end', &
         form='unformatted',status=file_stat,position='append')
  endif

  ! Write current time to file
  !***************************

  if (lroot) write(unitpartout) itime ! MPI root process only
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

  ! Potential vorticity, specific humidity, temperature, and density
  !*****************************************************************

      do il=2,nz
        if (height(il).gt.ztra1(i)) then
          indz=il-1
          indzp=il
          goto 6
        endif
      end do
6     continue

      dz1=ztra1(i)-height(indz)
      dz2=height(indzp)-ztra1(i)
      dz=1./(dz1+dz2)


      do ind=indz,indzp
        do m=1,2
          indexh=memind(m)

  ! Potential vorticity
          pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
               +p2*pv(ixp,jy ,ind,indexh) &
               +p3*pv(ix ,jyp,ind,indexh) &
               +p4*pv(ixp,jyp,ind,indexh)
  ! Specific humidity
          qv1(m)=p1*qv(ix ,jy ,ind,indexh) &
               +p2*qv(ixp,jy ,ind,indexh) &
               +p3*qv(ix ,jyp,ind,indexh) &
               +p4*qv(ixp,jyp,ind,indexh)
  ! Temperature
          tt1(m)=p1*tt(ix ,jy ,ind,indexh) &
               +p2*tt(ixp,jy ,ind,indexh) &
               +p3*tt(ix ,jyp,ind,indexh) &
               +p4*tt(ixp,jyp,ind,indexh)
  ! Density
          rho1(m)=p1*rho(ix ,jy ,ind,indexh) &
               +p2*rho(ixp,jy ,ind,indexh) &
               +p3*rho(ix ,jyp,ind,indexh) &
               +p4*rho(ixp,jyp,ind,indexh)
        end do
        pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
        qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
        ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
        rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
      end do
      pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
      qvi=(dz1*qvprof(2)+dz2*qvprof(1))*dz
      tti=(dz1*ttprof(2)+dz2*ttprof(1))*dz
      rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz

  ! Tropopause and PBL height
  !**************************

      do m=1,2
        indexh=memind(m)

  ! Tropopause
        tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
             + p2*tropopause(ixp,jy ,1,indexh) &
             + p3*tropopause(ix ,jyp,1,indexh) &
             + p4*tropopause(ixp,jyp,1,indexh)

  ! PBL height
        hm(m)=p1*hmix(ix ,jy ,1,indexh) &
             + p2*hmix(ixp,jy ,1,indexh) &
             + p3*hmix(ix ,jyp,1,indexh) &
             + p4*hmix(ixp,jyp,1,indexh)
      end do

      hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
      tri=(tr(1)*dt2+tr(2)*dt1)*dtt


  ! Write the output
  !*****************
      if (mp_measure_time) call mpif_mtime('iotime',0)

      write(unitpartout) npoint(i),xlon,ylat,ztra1(i), &
           itramem(i),topo,pvi,qvi,rhoi,hmixi,tri,tti, &
           (xmass1(i,j),j=1,nspec)

      if (mp_measure_time) call mpif_mtime('iotime',1)
    endif
  end do

  ! Only last MPI process writes EOF info
  if (mp_partid.eq.mp_partgroup_np-1) then
    write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
         -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
         (-9999.9,j=1,nspec)
  end if


  close(unitpartout)

end subroutine partoutput
