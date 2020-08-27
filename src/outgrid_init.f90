! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

! DJM - 2017-05-09 - added #ifdef USE_MPIINPLACE cpp directive to     *
! enable allocation of a gridunc0 array if required by MPI code in    *
! mpi_mod.f90                                                         *
!                                                                     *
!**********************************************************************


subroutine outgrid_init
  !
  !*****************************************************************************
  !                                                                            *
  !  This routine initializes the output grids                                 *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! area               surface area of all output grid cells                   *
  ! areaeast           eastward facing wall area of all output grid cells      *
  ! areanorth          northward facing wall area of all output grid cells     *
  ! volume             volumes of all output grid cells                        *
  !                                                                            *
  !*****************************************************************************

  use flux_mod
  use oh_mod
  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: ix,jy,kz,i,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
  integer :: ks,kp,stat
  real :: ylat,gridarea,ylatp,ylatm,hzone,cosfactm,cosfactp
  real :: xlon,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  real,parameter :: eps=nxmax/3.e5


  ! Compute surface area and volume of each grid cell: area, volume;
  ! and the areas of the northward and eastward facing walls: areaeast, areanorth
  !***********************************************************************
  do jy=0,numygrid-1
    ylat=outlat0+(real(jy)+0.5)*dyout
    ylatp=ylat+0.5*dyout
    ylatm=ylat-0.5*dyout
    if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
      hzone=dyout*r_earth*pi180
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

    gridarea=2.*pi*r_earth*hzone*dxout/360.

    do ix=0,numxgrid-1
      area(ix,jy)=gridarea

  ! Volume = area x box height
  !***************************

      volume(ix,jy,1)=area(ix,jy)*outheight(1)
      areaeast(ix,jy,1)=dyout*r_earth*pi180*outheight(1)
      areanorth(ix,jy,1)=cos(ylat*pi180)*dxout*r_earth*pi180* &
           outheight(1)
      do kz=2,numzgrid
        areaeast(ix,jy,kz)=dyout*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
        areanorth(ix,jy,kz)=cos(ylat*pi180)*dxout*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
        volume(ix,jy,kz)=area(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do




  !******************************************************************
  ! Determine average height of model topography in output grid cells
  !******************************************************************

  ! Loop over all output grid cells
  !********************************

  do jjy=0,numygrid-1
    do iix=0,numxgrid-1
      oroh=0.

  ! Take 100 samples of the topography in every grid cell
  !******************************************************

      do j1=1,10
        ylat=outlat0+(real(jjy)+real(j1)/10.-0.05)*dyout
        yl=(ylat-ylat0)/dy
        do i1=1,10
          xlon=outlon0+(real(iix)+real(i1)/10.-0.05)*dxout
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

      oroout(iix,jjy)=oroh/100.
    end do
  end do

  ! if necessary allocate flux fields
  if (iflux.eq.1) then
    allocate(flux(6,0:numxgrid-1,0:numygrid-1,numzgrid, &
         1:nspec,1:maxpointspec_act,1:nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate flux array '
  endif

  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(gridunc(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  if (ldirect.gt.0) then
    allocate(wetgridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc'
    allocate(drygridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc'
  endif

#ifdef USE_MPIINPLACE
#else
! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
! then an aux array is needed for parallel grid reduction
    allocate(gridunc0(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(gridunc0(1,1,1,1,1,1,1),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  end if
#endif
  if (ldirect.gt.0) then
    if (lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(0:numxgrid-1,0:numygrid-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc0'
      allocate(drygridunc0(0:numxgrid-1,0:numygrid-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc0'

! allocate a dummy to avoid compilator complaints
    else if (.not.lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(1,1,1,1,1,1),stat=stat)
      allocate(drygridunc0(1,1,1,1,1,1),stat=stat)
    end if
  end if

  !write (*,*) 'Dimensions for fields', numxgrid,numygrid, &
  !     maxspec,maxpointspec_act,nclassunc,maxageclass

  if (lroot) then
    write (*,*) 'Allocating fields for global output (x,y): ', numxgrid,numygrid
    write (*,*) 'Allocating fields for nested output (x,y): ', numxgridn,numygridn 
  end if

  ! allocate fields for concoutput with maximum dimension of outgrid
  ! and outgrid_nest

  allocate(gridsigma(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(grid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(densityoutgrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
! RLT
  allocate(densitydrygrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(factor_drygrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  allocate(factor3d(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(sparse_dump_r(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

   allocate(sparse_dump_u(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  allocate(sparse_dump_i(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  ! deposition fields are only allocated for forward runs
  if (ldirect.gt.0) then
     allocate(wetgridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(wetgrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  endif

  ! Initial condition field

  if (linit_cond.gt.0) then
    allocate(init_cond(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate init_cond'
    if (mpi_mode.gt.0) then
      if (lroot) then
        allocate(init_cond0(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
             maxpointspec_act),stat=stat)
      else
        allocate(init_cond0(1,1,1,1,1))
      end if
    end if
  endif

  !************************
  ! Initialize output grids
  !************************

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
  ! Receptor points
      creceptor(i,ks)=0.
    end do
    do nage=1,nageclass
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
  ! Deposition fields
            if (ldirect.gt.0) then
            wetgridunc(ix,jy,ks,kp,l,nage)=0.
            drygridunc(ix,jy,ks,kp,l,nage)=0.
            endif
            do kz=1,numzgrid
              if (iflux.eq.1) then
  ! Flux fields
                 do i=1,5
                   flux(i,ix,jy,kz,ks,kp,nage)=0.
                 end do
              endif
  ! Initial condition field
              if ((l.eq.1).and.(nage.eq.1).and.(linit_cond.gt.0)) &
                   init_cond(ix,jy,kz,ks,kp)=0.
  ! Concentration fields
              gridunc(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do



end subroutine outgrid_init
