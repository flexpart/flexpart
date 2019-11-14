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

subroutine boundcond_domainfill(itime,loutend)
!                                  i      i
!*****************************************************************************
!                                                                            *
! Particles are created by this subroutine continuously throughout the       *
! simulation at the boundaries of the domain-filling box.                    *
! All particles carry the same amount of mass which alltogether comprises the*
! mass of air within the box, which remains (more or less) constant.         *
!                                                                            *
!     Author: A. Stohl                                                       *
!                                                                            *
!     16 October 2002                                                        *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
!                                                                            *
! nx_we(2)       grid indices for western and eastern boundary of domain-    *
!                filling trajectory calculations                             *
! ny_sn(2)       grid indices for southern and northern boundary of domain-  *
!                filling trajectory calculations                             *
!                                                                            *
!*****************************************************************************
!  CHANGES                                                                   
!    08/2016 eso: MPI version:                                                
!
!   -Root process release particles and distributes to other processes.
!    Temporary arrays are used, also for the non-root (receiving) processes.
!   -The scheme can be improved by having all processes report numpart
!    (keeping track of how many particles have left the domain), so that
!    a proportional amount of new particles can be distributed (however
!    we have a separate function called from timemanager that will
!    redistribute particles among processes if there are imbalances)     
!*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use random_mod, only: ran1
  use mpi_mod

  implicit none

  real :: dz,dz1,dz2,dt1,dt2,dtt,ylat,xm,cosfact,accmasst
  integer :: itime,in,indz,indzp,i,loutend
  integer :: j,k,ix,jy,m,indzh,indexh,minpart,ipart,mmass
  integer :: numactiveparticles, numpart_total, rel_counter
  integer,allocatable,dimension(:) ::  numrel_mpi !, numactiveparticles_mpi

  real :: windl(2),rhol(2)
  real :: windhl(2),rhohl(2)
  real :: windx,rhox
  real :: deltaz,boundarea,fluxofmass

  integer :: ixm,ixp,jym,jyp,indzm,mm
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2),yh1(2)

  integer :: idummy = -11
  integer :: mtag
  logical :: first_call=.true.
! Sizes of temporary arrays are maxpartfract*maxpart. Increase maxpartfract if 
! needed.  
  real,parameter :: maxpartfract=0.1
  integer :: tmp_size = int(maxpartfract*maxpart)

! Use different seed for each process
  if (first_call) then
    idummy=idummy+mp_seed
    first_call=.false.
  end if


! If domain-filling is global, no boundary conditions are needed
!***************************************************************

  if (gdomainfill) return

  accmasst=0.
  numactiveparticles=0
! Keep track of active particles on each process
  allocate(numrel_mpi(0:mp_partgroup_np-1))
! numactiveparticles_mpi(0:mp_partgroup_np-1)

! New particles to be released on each process
  numrel_mpi(:)=0

! Terminate trajectories that have left the domain, if domain-filling
! trajectory calculation domain is not global. Done for all processes
!********************************************************************

  do i=1,numpart
    if (itra1(i).eq.itime) then
      if ((ytra1(i).gt.real(ny_sn(2))).or. &
           (ytra1(i).lt.real(ny_sn(1)))) itra1(i)=-999999999
      if (((.not.xglobal).or.(nx_we(2).ne.(nx-2))).and. &
           ((xtra1(i).lt.real(nx_we(1))).or. &
           (xtra1(i).gt.real(nx_we(2))))) itra1(i)=-999999999
    endif
    if (itra1(i).ne.-999999999) numactiveparticles= &
         numactiveparticles+1
  end do
!  numactiveparticles_mpi(mp_partid) = numactiveparticles


! Collect number of active particles from all processes
  ! call MPI_Allgather(numactiveparticles, 1, MPI_INTEGER, &
  !      &numactiveparticles_mpi, 1, MPI_INTEGER, mp_comm_used, mp_ierr)


! Total number of new releases
  numpart_total = 0


! This section only done by root process
!***************************************

  if (lroot) then

! Use separate arrays for newly released particles
!*************************************************

    allocate(itra1_tmp(tmp_size),npoint_tmp(tmp_size),nclass_tmp(tmp_size),&
         & idt_tmp(tmp_size),itramem_tmp(tmp_size),itrasplit_tmp(tmp_size),&
         & xtra1_tmp(tmp_size),ytra1_tmp(tmp_size),ztra1_tmp(tmp_size),&
         & xmass1_tmp(tmp_size, maxspec))

! Initialize all particles as non-existent    
    itra1_tmp(:)=-999999999

! Determine auxiliary variables for time interpolation
!*****************************************************

    dt1=real(itime-memtime(1))
    dt2=real(memtime(2)-itime)
    dtt=1./(dt1+dt2)

! Initialize auxiliary variable used to search for vacant storage space
!**********************************************************************

    minpart=1

!***************************************
! Western and eastern boundary condition
!***************************************

! Loop from south to north
!*************************

    do jy=ny_sn(1),ny_sn(2)

! Loop over western (index 1) and eastern (index 2) boundary
!***********************************************************

      do k=1,2

! Loop over all release locations in a column
!********************************************

        do j=1,numcolumn_we(k,jy)

! Determine, for each release location, the area of the corresponding boundary
!*****************************************************************************

          if (j.eq.1) then
            deltaz=(zcolumn_we(k,jy,2)+zcolumn_we(k,jy,1))/2.
          else if (j.eq.numcolumn_we(k,jy)) then
!        deltaz=height(nz)-(zcolumn_we(k,jy,j-1)+
!    +        zcolumn_we(k,jy,j))/2.
! In order to avoid taking a very high column for very many particles,
! use the deltaz from one particle below instead
            deltaz=(zcolumn_we(k,jy,j)-zcolumn_we(k,jy,j-2))/2.
          else
            deltaz=(zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))/2.
          endif
          if ((jy.eq.ny_sn(1)).or.(jy.eq.ny_sn(2))) then
            boundarea=deltaz*111198.5/2.*dy
          else
            boundarea=deltaz*111198.5*dy
          endif


! Interpolate the wind velocity and density to the release location
!******************************************************************

! Determine the model level below the release position
!*****************************************************

          do i=2,nz
            if (height(i).gt.zcolumn_we(k,jy,j)) then
              indz=i-1
              indzp=i
              goto 6
            endif
          end do
6         continue

! Vertical distance to the level below and above current position
!****************************************************************

          dz1=zcolumn_we(k,jy,j)-height(indz)
          dz2=height(indzp)-zcolumn_we(k,jy,j)
          dz=1./(dz1+dz2)

! Vertical and temporal interpolation
!************************************

          do m=1,2
            indexh=memind(m)
            do in=1,2
              indzh=indz+in-1
              windl(in)=uu(nx_we(k),jy,indzh,indexh)
              rhol(in)=rho(nx_we(k),jy,indzh,indexh)
            end do

            windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
            rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
          end do

          windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
          rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

! Calculate mass flux
!********************

          fluxofmass=windx*rhox*boundarea*real(lsynctime)


! If the mass flux is directed into the domain, add it to previous mass fluxes;
! if it is out of the domain, set accumulated mass flux to zero
!******************************************************************************

          if (k.eq.1) then
            if (fluxofmass.ge.0.) then
              acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+fluxofmass
            else
              acc_mass_we(k,jy,j)=0.
            endif
          else
            if (fluxofmass.le.0.) then
              acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+abs(fluxofmass)
            else
              acc_mass_we(k,jy,j)=0.
            endif
          endif
          accmasst=accmasst+acc_mass_we(k,jy,j)

! If the accumulated mass exceeds half the mass that each particle shall carry,
! one (or more) particle(s) is (are) released and the accumulated mass is
! reduced by the mass of this (these) particle(s)
!******************************************************************************

          if (acc_mass_we(k,jy,j).ge.xmassperparticle/2.) then
            mmass=int((acc_mass_we(k,jy,j)+xmassperparticle/2.)/ &
                 xmassperparticle)
            acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)- &
                 real(mmass)*xmassperparticle
          else
            mmass=0
          endif

          do m=1,mmass
            do ipart=minpart,maxpart

! If a vacant storage space is found, attribute everything to this array element
! TODO: for the MPI version this test can be removed, as all 
!       elements in _tmp arrays are initialized to zero
!*****************************************************************************

              if (itra1_tmp(ipart).ne.itime) then

! Assign particle positions
!**************************

                xtra1_tmp(ipart)=real(nx_we(k))
                if (jy.eq.ny_sn(1)) then
                  ytra1_tmp(ipart)=real(jy)+0.5*ran1(idummy)
                else if (jy.eq.ny_sn(2)) then
                  ytra1_tmp(ipart)=real(jy)-0.5*ran1(idummy)
                else
                  ytra1_tmp(ipart)=real(jy)+(ran1(idummy)-.5)
                endif
                if (j.eq.1) then
                  ztra1_tmp(ipart)=zcolumn_we(k,jy,1)+(zcolumn_we(k,jy,2)- &
                       zcolumn_we(k,jy,1))/4.
                else if (j.eq.numcolumn_we(k,jy)) then
                  ztra1_tmp(ipart)=(2.*zcolumn_we(k,jy,j)+ &
                       zcolumn_we(k,jy,j-1)+height(nz))/4.
                else
                  ztra1_tmp(ipart)=zcolumn_we(k,jy,j-1)+ran1(idummy)* &
                       (zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))
                endif

! Interpolate PV to the particle position
!****************************************
                ixm=int(xtra1_tmp(ipart))
                jym=int(ytra1_tmp(ipart))
                ixp=ixm+1
                jyp=jym+1
                ddx=xtra1_tmp(ipart)-real(ixm)
                ddy=ytra1_tmp(ipart)-real(jym)
                rddx=1.-ddx
                rddy=1.-ddy
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy
                do i=2,nz
                  if (height(i).gt.ztra1_tmp(ipart)) then
                    indzm=i-1
                    indzp=i
                    goto 26
                  endif
                end do
26              continue
                dz1=ztra1_tmp(ipart)-height(indzm)
                dz2=height(indzp)-ztra1_tmp(ipart)
                dz=1./(dz1+dz2)
                do mm=1,2
                  indexh=memind(mm)
                  do in=1,2
                    indzh=indzm+in-1
                    y1(in)=p1*pv(ixm,jym,indzh,indexh) &
                         +p2*pv(ixp,jym,indzh,indexh) &
                         +p3*pv(ixm,jyp,indzh,indexh) &
                         +p4*pv(ixp,jyp,indzh,indexh)
                  end do
                  yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
                end do
                pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
                ylat=ylat0+ytra1_tmp(ipart)*dy
                if (ylat.lt.0.) pvpart=-1.*pvpart


! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
!*****************************************************************************

                if (((ztra1_tmp(ipart).gt.3000.).and. &
                     (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                  nclass_tmp(ipart)=min(int(ran1(idummy)* &
                       real(nclassunc))+1,nclassunc)
                  numactiveparticles=numactiveparticles+1
                  numparticlecount=numparticlecount+1
                  npoint_tmp(ipart)=numparticlecount
                  idt_tmp(ipart)=mintime
                  itra1_tmp(ipart)=itime
                  itramem_tmp(ipart)=itra1_tmp(ipart)
                  itrasplit_tmp(ipart)=itra1_tmp(ipart)+ldirect*itsplit
                  xmass1_tmp(ipart,1)=xmassperparticle
                  if (mdomainfill.eq.2) xmass1_tmp(ipart,1)= &
                       xmass1_tmp(ipart,1)*pvpart*48./29.*ozonescale/10.**9
                else
                  goto 71
                endif


! Increase numpart, if necessary
!*******************************

                numpart_total=max(numpart_total,ipart)
                goto 73      ! Storage space has been found, stop searching
              endif
            end do
            if (ipart.gt.tmp_size) &
                 stop 'boundcond_domainfill_mpi.f90: too many particles required'
73          minpart=ipart+1
71          continue
          end do


        end do
      end do
    end do


!*****************************************
! Southern and northern boundary condition
!*****************************************

! Loop from west to east
!***********************

    do ix=nx_we(1),nx_we(2)

! Loop over southern (index 1) and northern (index 2) boundary
!*************************************************************

      do k=1,2
        ylat=ylat0+real(ny_sn(k))*dy
        cosfact=cos(ylat*pi180)

! Loop over all release locations in a column
!********************************************

        do j=1,numcolumn_sn(k,ix)

! Determine, for each release location, the area of the corresponding boundary
!*****************************************************************************

          if (j.eq.1) then
            deltaz=(zcolumn_sn(k,ix,2)+zcolumn_sn(k,ix,1))/2.
          else if (j.eq.numcolumn_sn(k,ix)) then
!        deltaz=height(nz)-(zcolumn_sn(k,ix,j-1)+
!    +        zcolumn_sn(k,ix,j))/2.
! In order to avoid taking a very high column for very many particles,
! use the deltaz from one particle below instead
            deltaz=(zcolumn_sn(k,ix,j)-zcolumn_sn(k,ix,j-2))/2.
          else
            deltaz=(zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))/2.
          endif
          if ((ix.eq.nx_we(1)).or.(ix.eq.nx_we(2))) then
            boundarea=deltaz*111198.5/2.*cosfact*dx
          else
            boundarea=deltaz*111198.5*cosfact*dx
          endif


! Interpolate the wind velocity and density to the release location
!******************************************************************

! Determine the model level below the release position
!*****************************************************

          do i=2,nz
            if (height(i).gt.zcolumn_sn(k,ix,j)) then
              indz=i-1
              indzp=i
              goto 16
            endif
          end do
16        continue

! Vertical distance to the level below and above current position
!****************************************************************

          dz1=zcolumn_sn(k,ix,j)-height(indz)
          dz2=height(indzp)-zcolumn_sn(k,ix,j)
          dz=1./(dz1+dz2)

! Vertical and temporal interpolation
!************************************

          do m=1,2
            indexh=memind(m)
            do in=1,2
              indzh=indz+in-1
              windl(in)=vv(ix,ny_sn(k),indzh,indexh)
              rhol(in)=rho(ix,ny_sn(k),indzh,indexh)
            end do

            windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
            rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
          end do

          windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
          rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

! Calculate mass flux
!********************

          fluxofmass=windx*rhox*boundarea*real(lsynctime)

! If the mass flux is directed into the domain, add it to previous mass fluxes;
! if it is out of the domain, set accumulated mass flux to zero
!******************************************************************************

          if (k.eq.1) then
            if (fluxofmass.ge.0.) then
              acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+fluxofmass
            else
              acc_mass_sn(k,ix,j)=0.
            endif
          else
            if (fluxofmass.le.0.) then
              acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+abs(fluxofmass)
            else
              acc_mass_sn(k,ix,j)=0.
            endif
          endif
          accmasst=accmasst+acc_mass_sn(k,ix,j)

! If the accumulated mass exceeds half the mass that each particle shall carry,
! one (or more) particle(s) is (are) released and the accumulated mass is
! reduced by the mass of this (these) particle(s)
!******************************************************************************

          if (acc_mass_sn(k,ix,j).ge.xmassperparticle/2.) then
            mmass=int((acc_mass_sn(k,ix,j)+xmassperparticle/2.)/ &
                 xmassperparticle)
            acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)- &
                 real(mmass)*xmassperparticle
          else
            mmass=0
          endif

          do m=1,mmass
            do ipart=minpart,maxpart

! If a vacant storage space is found, attribute everything to this array element
!*****************************************************************************

              if (itra1_tmp(ipart).ne.itime) then

! Assign particle positions
!**************************

                ytra1_tmp(ipart)=real(ny_sn(k))
                if (ix.eq.nx_we(1)) then
                  xtra1_tmp(ipart)=real(ix)+0.5*ran1(idummy)
                else if (ix.eq.nx_we(2)) then
                  xtra1_tmp(ipart)=real(ix)-0.5*ran1(idummy)
                else
                  xtra1_tmp(ipart)=real(ix)+(ran1(idummy)-.5)
                endif
                if (j.eq.1) then
                  ztra1_tmp(ipart)=zcolumn_sn(k,ix,1)+(zcolumn_sn(k,ix,2)- &
                       zcolumn_sn(k,ix,1))/4.
                else if (j.eq.numcolumn_sn(k,ix)) then
                  ztra1_tmp(ipart)=(2.*zcolumn_sn(k,ix,j)+ &
                       zcolumn_sn(k,ix,j-1)+height(nz))/4.
                else
                  ztra1_tmp(ipart)=zcolumn_sn(k,ix,j-1)+ran1(idummy)* &
                       (zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))
                endif


! Interpolate PV to the particle position
!****************************************
                ixm=int(xtra1_tmp(ipart))
                jym=int(ytra1_tmp(ipart))
                ixp=ixm+1
                jyp=jym+1
                ddx=xtra1_tmp(ipart)-real(ixm)
                ddy=ytra1_tmp(ipart)-real(jym)
                rddx=1.-ddx
                rddy=1.-ddy
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy
                do i=2,nz
                  if (height(i).gt.ztra1_tmp(ipart)) then
                    indzm=i-1
                    indzp=i
                    goto 126
                  endif
                end do
126             continue
                dz1=ztra1_tmp(ipart)-height(indzm)
                dz2=height(indzp)-ztra1_tmp(ipart)
                dz=1./(dz1+dz2)
                do mm=1,2
                  indexh=memind(mm)
                  do in=1,2
                    indzh=indzm+in-1
                    y1(in)=p1*pv(ixm,jym,indzh,indexh) &
                         +p2*pv(ixp,jym,indzh,indexh) &
                         +p3*pv(ixm,jyp,indzh,indexh) &
                         +p4*pv(ixp,jyp,indzh,indexh)
                  end do
                  yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
                end do
                pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
                if (ylat.lt.0.) pvpart=-1.*pvpart


! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
!*****************************************************************************

                if (((ztra1_tmp(ipart).gt.3000.).and. &
                     (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                  nclass_tmp(ipart)=min(int(ran1(idummy)* &
                       real(nclassunc))+1,nclassunc)
                  numactiveparticles=numactiveparticles+1
                  numparticlecount=numparticlecount+1
                  npoint_tmp(ipart)=numparticlecount
                  idt_tmp(ipart)=mintime
                  itra1_tmp(ipart)=itime
                  itramem_tmp(ipart)=itra1_tmp(ipart)
                  itrasplit_tmp(ipart)=itra1_tmp(ipart)+ldirect*itsplit
                  xmass1_tmp(ipart,1)=xmassperparticle
                  if (mdomainfill.eq.2) xmass1_tmp(ipart,1)= &
                       xmass1_tmp(ipart,1)*pvpart*48./29.*ozonescale/10.**9
                else
                  goto 171
                endif


! Increase numpart, if necessary
!*******************************
                numpart_total=max(numpart_total,ipart)
                goto 173      ! Storage space has been found, stop searching
              endif
            end do
            if (ipart.gt.tmp_size) &
                 stop 'boundcond_domainfill.f: too many particles required'
173         minpart=ipart+1
171         continue
          end do


        end do
      end do
    end do


! xm=0.
! do i=1,numpart_total
!   if (itra1_tmp(i).eq.itime) xm=xm+xmass1(i,1)
! end do

!write(*,*) itime,numactiveparticles,numparticlecount,numpart,
!    +xm,accmasst,xm+accmasst

  end if ! if lroot

! Distribute the number of particles to be released
! *************************************************
  call MPI_Bcast(numpart_total, 1, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)

  do i=0, mp_partgroup_np-1
    numrel_mpi(i) = numpart_total/mp_partgroup_np
    if (i.lt.mod(numpart_total,mp_partgroup_np)) numrel_mpi(i) = numrel_mpi(i) + 1
  end do

! Allocate temporary arrays for receiving processes
  if (.not.lroot) then 
    allocate(itra1_tmp(numrel_mpi(mp_partid)),&
         & npoint_tmp(numrel_mpi(mp_partid)),&
         & nclass_tmp(numrel_mpi(mp_partid)),&
         & idt_tmp(numrel_mpi(mp_partid)),&
         & itramem_tmp(numrel_mpi(mp_partid)),&
         & itrasplit_tmp(numrel_mpi(mp_partid)),&
         & xtra1_tmp(numrel_mpi(mp_partid)),&
         & ytra1_tmp(numrel_mpi(mp_partid)),&
         & ztra1_tmp(numrel_mpi(mp_partid)),&
         & xmass1_tmp(numrel_mpi(mp_partid),maxspec))
      
! Initialize all particles as non-existent    
    itra1_tmp(:)=-999999999
  end if

! Distribute particles
! Keep track of released particles so far
  rel_counter = 0
  mtag = 1000

  do i=0, mp_partgroup_np-1

! For root process, nothing to do except update release count
    if (i.eq.0) then 
      rel_counter = rel_counter + numrel_mpi(i)
      cycle
    end if
      
! Send particles from root to non-root processes
    if (lroot.and.numrel_mpi(i).gt.0) then

      call MPI_SEND(nclass_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+1*i,mp_comm_used,mp_ierr)

      call MPI_SEND(npoint_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+2*i,mp_comm_used,mp_ierr)

      call MPI_SEND(itra1_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+3*i,mp_comm_used,mp_ierr)

      call MPI_SEND(idt_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+4*i,mp_comm_used,mp_ierr)

      call MPI_SEND(itramem_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+5*i,mp_comm_used,mp_ierr)

      call MPI_SEND(itrasplit_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),MPI_INTEGER,i,mtag+6*i,mp_comm_used,mp_ierr)

      call MPI_SEND(xtra1_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),mp_dp,i,mtag+7*i,mp_comm_used,mp_ierr)

      call MPI_SEND(ytra1_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),mp_dp,i,mtag+8*i,mp_comm_used,mp_ierr)

      call MPI_SEND(ztra1_tmp(rel_counter+1:rel_counter+numrel_mpi(i)),&
           &numrel_mpi(i),mp_sp,i,mtag+9*i,mp_comm_used,mp_ierr)

      do j=1,nspec
        call MPI_SEND(xmass1_tmp(rel_counter+1:rel_counter+numrel_mpi(i),j),&
             &numrel_mpi(i),mp_sp,i,mtag+(9+j)*i,mp_comm_used,mp_ierr)
      end do

! Non-root processes issue receive requests
    else if (i.eq.mp_partid.and.numrel_mpi(i).gt.0) then 
      call MPI_RECV(nclass_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+1*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(npoint_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+2*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itra1_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+3*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(idt_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+4*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itramem_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+5*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itrasplit_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &MPI_INTEGER,id_root,mtag+6*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(xtra1_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &mp_dp,id_root,mtag+7*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(ytra1_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &mp_dp,id_root,mtag+8*i,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(ztra1_tmp(1:numrel_mpi(i)),numrel_mpi(i),&
           &mp_sp,id_root,mtag+9*i,mp_comm_used,mp_status,mp_ierr)

      do j=1,nspec
        call MPI_RECV(xmass1_tmp(1:numrel_mpi(i),j),numrel_mpi(i),&
             &mp_sp,id_root,mtag+(9+j)*i,mp_comm_used,mp_status,mp_ierr)

      end do
    end if
    rel_counter = rel_counter + numrel_mpi(i)
  end do

! Find free storage space for the new particles.
! This section is independent of the redistribution scheme used
! ********************************************************************

! Keep track of released particles so far
  minpart=1
    
! The algorithm should be correct also for root process
  do i=1, numrel_mpi(mp_partid)
    do ipart=minpart, maxpart
      if (itra1(ipart).ne.itime) then
        itra1(ipart) = itra1_tmp(i)
        npoint(ipart) = npoint_tmp(i)
        nclass(ipart) = nclass_tmp(i)
        idt(ipart) = idt_tmp(i)
        itramem(ipart) = itramem_tmp(i)
        itrasplit(ipart) = itrasplit_tmp(i)
        xtra1(ipart) = xtra1_tmp(i)
        ytra1(ipart) = ytra1_tmp(i)
        ztra1(ipart) = ztra1_tmp(i)
        xmass1(ipart,:) = xmass1_tmp(i,:)
! Increase numpart, if necessary
        numpart=max(numpart,ipart)
        goto 200 ! Storage space has been found, stop searching
      end if
    end do
200 minpart=ipart+1
  end do

! If particles shall be dumped, then accumulated masses at the domain boundaries
! must be dumped, too, to be used for later runs
!*****************************************************************************

  if ((ipout.gt.0).and.(itime.eq.loutend)) then
    if (lroot) then
      call mpif_mtime('iotime',0)
      open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
           form='unformatted')
      write(unitboundcond) numcolumn_we,numcolumn_sn, &
           zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
      close(unitboundcond)
      call mpif_mtime('iotime',1)
    end if
  endif

! Deallocate temporary arrays 
  deallocate(itra1_tmp,npoint_tmp,nclass_tmp,idt_tmp,itramem_tmp,itrasplit_tmp,&
       & xtra1_tmp,ytra1_tmp,ztra1_tmp,xmass1_tmp,numrel_mpi)
! numactiveparticles_mpi


end subroutine boundcond_domainfill
