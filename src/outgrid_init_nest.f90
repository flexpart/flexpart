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

subroutine outgrid_init_nest

!*****************************************************************************
!                                                                            *
!  This routine calculates, for each grid cell of the output nest, the       *
!  volume and the surface area.                                              *
!                                                                            *
!     Author: A. Stohl                                                       *
!                                                                            *
!    30 August 2004                                                          *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
!                                                                            *
! arean              surface area of all output nest cells                   *
! volumen            volumes of all output nest cells                        *
!                                                                            *
!*****************************************************************************

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: ix,jy,kz,ks,kp,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
  integer :: stat
  real :: ylat,gridarea,ylatp,ylatm,hzone,cosfactm,cosfactp
  real :: xlon,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  real,parameter :: eps=nxmax/3.e5



! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(griduncn(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'

  if (ldirect.gt.0) then
    allocate(wetgriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
    allocate(drygriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  endif

! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
    ! allocate(griduncn0(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
    !      maxpointspec_act,nclassunc,maxageclass),stat=stat)
    ! if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'

    if (ldirect.gt.0) then
      allocate(wetgriduncn0(0:numxgridn-1,0:numygridn-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
      allocate(drygriduncn0(0:numxgridn-1,0:numygridn-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
    endif
! allocate a dummy to avoid compilator complaints
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(wetgriduncn0(1,1,1,1,1,1),stat=stat)
    allocate(drygriduncn0(1,1,1,1,1,1),stat=stat)
  end if

! Compute surface area and volume of each grid cell: area, volume;
! and the areas of the northward and eastward facing walls: areaeast, areanorth
!***********************************************************************

  do jy=0,numygridn-1
    ylat=outlat0n+(real(jy)+0.5)*dyoutn
    ylatp=ylat+0.5*dyoutn
    ylatm=ylat-0.5*dyoutn
    if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
      hzone=dyoutn*r_earth*pi180
    else

! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
!************************************************************

      cosfactp=cos(ylatp*pi180)
      cosfactm=cos(ylatm*pi180)
      if (cosfactp.lt.cosfactm) then
        hzone=sqrt(1-cosfactp**2)- &
             sqrt(1-cosfactm**2)
        hzone=hzone*r_earth
      else
        hzone=sqrt(1-cosfactm**2)- &
             sqrt(1-cosfactp**2)
        hzone=hzone*r_earth
      endif
    endif



! Surface are of a grid cell at a latitude ylat
!**********************************************

    gridarea=2.*pi*r_earth*hzone*dxoutn/360.

    do ix=0,numxgridn-1
      arean(ix,jy)=gridarea

! Volume = area x box height
!***************************

      volumen(ix,jy,1)=arean(ix,jy)*outheight(1)
      do kz=2,numzgrid
        volumen(ix,jy,kz)=arean(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do


!**************************************************************************
! Determine average height of model topography in nesteed output grid cells
!**************************************************************************

! Loop over all output grid cells
!********************************

  do jjy=0,numygridn-1
    do iix=0,numxgridn-1
      oroh=0.

! Take 100 samples of the topography in every grid cell
!******************************************************

      do j1=1,10
        ylat=outlat0n+(real(jjy)+real(j1)/10.-0.05)*dyoutn
        yl=(ylat-ylat0)/dy
        do i1=1,10
          xlon=outlon0n+(real(iix)+real(i1)/10.-0.05)*dxoutn
          xl=(xlon-xlon0)/dx

! Determine the nest we are in
!*****************************

          ngrid=0
          do j=numbnests,1,-1
            if ((xl.gt.xln(j)+eps).and.(xl.lt.xrn(j)-eps).and. &
                 (yl.gt.yln(j)+eps).and.(yl.lt.yrn(j)-eps)) then
              ngrid=j
              goto 43
            endif
          end do
43        continue

! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
!*****************************************************************************

          if (ngrid.gt.0) then
            xtn=(xl-xln(ngrid))*xresoln(ngrid)
            ytn=(yl-yln(ngrid))*yresoln(ngrid)
            ix=int(xtn)
            jy=int(ytn)
            ddy=ytn-real(jy)
            ddx=xtn-real(ix)
          else
            ix=int(xl)
            jy=int(yl)
            ddy=yl-real(jy)
            ddx=xl-real(ix)
          endif
          ixp=ix+1
          jyp=jy+1
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

          if (ngrid.gt.0) then
            oroh=oroh+p1*oron(ix ,jy ,ngrid) &
                 + p2*oron(ixp,jy ,ngrid) &
                 + p3*oron(ix ,jyp,ngrid) &
                 + p4*oron(ixp,jyp,ngrid)
          else
            oroh=oroh+p1*oro(ix ,jy) &
                 + p2*oro(ixp,jy) &
                 + p3*oro(ix ,jyp) &
                 + p4*oro(ixp,jyp)
          endif
        end do
      end do

! Divide by the number of samples taken
!**************************************

      orooutn(iix,jjy)=oroh/100.
    end do
  end do



!*******************************
! Initialization of output grids
!*******************************

  do kp=1,maxpointspec_act
    do ks=1,nspec
      do nage=1,nageclass
        do jy=0,numygridn-1
          do ix=0,numxgridn-1
            do l=1,nclassunc
! Deposition fields
              if (ldirect.gt.0) then
                wetgriduncn(ix,jy,ks,kp,l,nage)=0.
                drygriduncn(ix,jy,ks,kp,l,nage)=0.
              endif
! Concentration fields
              do kz=1,numzgrid
                griduncn(ix,jy,kz,ks,kp,l,nage)=0.
              end do
            end do
          end do
        end do
      end do
    end do
  end do


end subroutine outgrid_init_nest
