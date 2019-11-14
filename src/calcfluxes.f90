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

subroutine calcfluxes(nage,jpart,xold,yold,zold)
  !                       i     i    i    i    i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the gross fluxes across horizontal, eastward and        *
  !     northward facing surfaces. The routine calculates the mass flux        *
  !     due to the motion of only one particle. The fluxes of subsequent calls *
  !     to this subroutine are accumulated until the next output is due.       *
  !     Upon output, flux fields are re-set to zero in subroutine fluxoutput.f.*
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     04 April 2000                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nage                  Age class of the particle considered                 *
  ! jpart                 Index of the particle considered                     *
  ! xold,yold,zold        "Memorized" old positions of the particle            *
  !                                                                            *
  !*****************************************************************************

  use flux_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jpart,nage,ixave,jyave,kz,kzave,kp
  integer :: k,k1,k2,ix,ix1,ix2,ixs,jy,jy1,jy2
  real :: xold,yold,zold,xmean,ymean


  ! Determine average positions
  !****************************

  if ((ioutputforeachrelease.eq.1).and.(mdomainfill.eq.0)) then
     kp=npoint(jpart)
  else
     kp=1
  endif

  xmean=(xold+xtra1(jpart))/2.
  ymean=(yold+ytra1(jpart))/2.

  ixave=int((xmean*dx+xoutshift)/dxout)
  jyave=int((ymean*dy+youtshift)/dyout)
  do kz=1,numzgrid                ! determine height of cell
    if (outheight(kz).gt.ztra1(jpart)) goto 16
  end do
16   kzave=kz


  ! Determine vertical fluxes
  !**************************

  if ((ixave.ge.0).and.(jyave.ge.0).and.(ixave.le.numxgrid-1).and. &
       (jyave.le.numygrid-1)) then
    do kz=1,numzgrid                ! determine height of cell
      if (outheighthalf(kz).gt.zold) goto 11
    end do
11   k1=min(numzgrid,kz)
    do kz=1,numzgrid                ! determine height of cell
      if (outheighthalf(kz).gt.ztra1(jpart)) goto 21
    end do
21   k2=min(numzgrid,kz)

    do k=1,nspec
      do kz=k1,k2-1
        flux(5,ixave,jyave,kz,k,kp,nage)= &
             flux(5,ixave,jyave,kz,k,kp,nage)+ &
             xmass1(jpart,k)
      end do
      do kz=k2,k1-1
        flux(6,ixave,jyave,kz,k,kp,nage)= &
             flux(6,ixave,jyave,kz,k,kp,nage)+ &
             xmass1(jpart,k)
      end do
    end do
  endif


  ! Determine west-east fluxes (fluxw) and east-west fluxes (fluxe)
  !****************************************************************

  if ((kzave.le.numzgrid).and.(jyave.ge.0).and. &
       (jyave.le.numygrid-1)) then

  ! 1) Particle does not cross domain boundary

    if (abs(xold-xtra1(jpart)).lt.real(nx)/2.) then
      ix1=int((xold*dx+xoutshift)/dxout+0.5)
      ix2=int((xtra1(jpart)*dx+xoutshift)/dxout+0.5)
      do k=1,nspec
        do ix=ix1,ix2-1
          if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
            flux(1,ix,jyave,kzave,k,kp,nage)= &
                 flux(1,ix,jyave,kzave,k,kp,nage) &
                 +xmass1(jpart,k)
          endif
        end do
        do ix=ix2,ix1-1
          if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
            flux(2,ix,jyave,kzave,k,kp,nage)= &
                 flux(2,ix,jyave,kzave,k,kp,nage) &
                 +xmass1(jpart,k)
          endif
        end do
      end do

  ! 2) Particle crosses domain boundary: use cyclic boundary condition
  !    and attribute flux to easternmost grid row only (approximation valid
  !    for relatively slow motions compared to output grid cell size)

    else
      ixs=int(((real(nxmin1)-1.e5)*dx+xoutshift)/dxout)
      if ((ixs.ge.0).and.(ixs.le.numxgrid-1)) then
        if (xold.gt.xtra1(jpart)) then       ! west-east flux
          do k=1,nspec
            flux(1,ixs,jyave,kzave,k,kp,nage)= &
                 flux(1,ixs,jyave,kzave,k,kp,nage) &
                 +xmass1(jpart,k)
          end do
        else                                 ! east-west flux
          do k=1,nspec
            flux(2,ixs,jyave,kzave,k,kp,nage)= &
                 flux(2,ixs,jyave,kzave,k,kp,nage) &
                 +xmass1(jpart,k)
          end do
        endif
      endif
    endif
  endif


  ! Determine south-north fluxes (fluxs) and north-south fluxes (fluxn)
  !********************************************************************

  if ((kzave.le.numzgrid).and.(ixave.ge.0).and. &
       (ixave.le.numxgrid-1)) then
    jy1=int((yold*dy+youtshift)/dyout+0.5)
    jy2=int((ytra1(jpart)*dy+youtshift)/dyout+0.5)

    do k=1,nspec
      do jy=jy1,jy2-1
        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          flux(3,ixave,jy,kzave,k,kp,nage)= &
               flux(3,ixave,jy,kzave,k,kp,nage) &
               +xmass1(jpart,k)
        endif
      end do
      do jy=jy2,jy1-1
        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          flux(4,ixave,jy,kzave,k,kp,nage)= &
               flux(4,ixave,jy,kzave,k,kp,nage) &
               +xmass1(jpart,k)
        endif
      end do
    end do
  endif

end subroutine calcfluxes

