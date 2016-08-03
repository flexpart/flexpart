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

  use point_mod
  use par_mod
  use com_mod
  use random_mod, only: ran1
  use mpi_mod

  implicit none

  real :: dz,dz1,dz2,dt1,dt2,dtt,ylat,xm,cosfact,accmasst
  integer :: itime,in,indz,indzp,i,loutend
  integer :: j,k,ix,jy,m,indzh,indexh,minpart,ipart,mmass
  integer :: numactiveparticles

  real :: windl(2),rhol(2)
  real :: windhl(2),rhohl(2)
  real :: windx,rhox
  real :: deltaz,boundarea,fluxofmass

  integer :: ixm,ixp,jym,jyp,indzm,mm
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2),yh1(2)

  integer :: idummy = -11
  logical :: first_call=.true.

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

  ! Terminate trajectories that have left the domain, if domain-filling
  ! trajectory calculation domain is not global
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
6       continue

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

  ! Calculate mass flux, divided by number of processes
  !****************************************************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)/mp_partgroup_np


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
          do ipart=minpart,maxpart_mpi

  ! If a vacant storage space is found, attribute everything to this array element
  !*****************************************************************************

            if (itra1(ipart).ne.itime) then

  ! Assign particle positions
  !**************************

              xtra1(ipart)=real(nx_we(k))
              if (jy.eq.ny_sn(1)) then
                ytra1(ipart)=real(jy)+0.5*ran1(idummy)
              else if (jy.eq.ny_sn(2)) then
                ytra1(ipart)=real(jy)-0.5*ran1(idummy)
              else
                ytra1(ipart)=real(jy)+(ran1(idummy)-.5)
              endif
              if (j.eq.1) then
                ztra1(ipart)=zcolumn_we(k,jy,1)+(zcolumn_we(k,jy,2)- &
                     zcolumn_we(k,jy,1))/4.
              else if (j.eq.numcolumn_we(k,jy)) then
                ztra1(ipart)=(2.*zcolumn_we(k,jy,j)+ &
                     zcolumn_we(k,jy,j-1)+height(nz))/4.
              else
                ztra1(ipart)=zcolumn_we(k,jy,j-1)+ran1(idummy)* &
                     (zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))
              endif

  ! Interpolate PV to the particle position
  !****************************************
              ixm=int(xtra1(ipart))
              jym=int(ytra1(ipart))
              ixp=ixm+1
              jyp=jym+1
              ddx=xtra1(ipart)-real(ixm)
              ddy=ytra1(ipart)-real(jym)
              rddx=1.-ddx
              rddy=1.-ddy
              p1=rddx*rddy
              p2=ddx*rddy
              p3=rddx*ddy
              p4=ddx*ddy
              do i=2,nz
                if (height(i).gt.ztra1(ipart)) then
                  indzm=i-1
                  indzp=i
                  goto 26
                endif
              end do
26            continue
              dz1=ztra1(ipart)-height(indzm)
              dz2=height(indzp)-ztra1(ipart)
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
              ylat=ylat0+ytra1(ipart)*dy
              if (ylat.lt.0.) pvpart=-1.*pvpart


  ! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
  !*****************************************************************************

              if (((ztra1(ipart).gt.3000.).and. &
                   (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                nclass(ipart)=min(int(ran1(idummy)* &
                     real(nclassunc))+1,nclassunc)
                numactiveparticles=numactiveparticles+1
                numparticlecount=numparticlecount+1
                npoint(ipart)=numparticlecount
                idt(ipart)=mintime
                itra1(ipart)=itime
                itramem(ipart)=itra1(ipart)
                itrasplit(ipart)=itra1(ipart)+ldirect*itsplit
                xmass1(ipart,1)=xmassperparticle
                if (mdomainfill.eq.2) xmass1(ipart,1)= &
                     xmass1(ipart,1)*pvpart*48./29.*ozonescale/10.**9
              else
                goto 71
              endif


  ! Increase numpart, if necessary
  !*******************************

              numpart=max(numpart,ipart)
              goto 73      ! Storage space has been found, stop searching
            endif
          end do
          if (ipart.gt.maxpart_mpi) &
               stop 'boundcond_domainfill.f: too many particles required'
73        minpart=ipart+1
71        continue
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
16      continue

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

  ! Calculate mass flux, divided by number of processes
  !****************************************************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)/mp_partgroup_np

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
          do ipart=minpart,maxpart_mpi

  ! If a vacant storage space is found, attribute everything to this array element
  !*****************************************************************************

            if (itra1(ipart).ne.itime) then

  ! Assign particle positions
  !**************************

              ytra1(ipart)=real(ny_sn(k))
              if (ix.eq.nx_we(1)) then
                xtra1(ipart)=real(ix)+0.5*ran1(idummy)
              else if (ix.eq.nx_we(2)) then
                xtra1(ipart)=real(ix)-0.5*ran1(idummy)
              else
                xtra1(ipart)=real(ix)+(ran1(idummy)-.5)
              endif
              if (j.eq.1) then
                ztra1(ipart)=zcolumn_sn(k,ix,1)+(zcolumn_sn(k,ix,2)- &
                     zcolumn_sn(k,ix,1))/4.
              else if (j.eq.numcolumn_sn(k,ix)) then
                ztra1(ipart)=(2.*zcolumn_sn(k,ix,j)+ &
                     zcolumn_sn(k,ix,j-1)+height(nz))/4.
              else
                ztra1(ipart)=zcolumn_sn(k,ix,j-1)+ran1(idummy)* &
                     (zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))
              endif


  ! Interpolate PV to the particle position
  !****************************************
              ixm=int(xtra1(ipart))
              jym=int(ytra1(ipart))
              ixp=ixm+1
              jyp=jym+1
              ddx=xtra1(ipart)-real(ixm)
              ddy=ytra1(ipart)-real(jym)
              rddx=1.-ddx
              rddy=1.-ddy
              p1=rddx*rddy
              p2=ddx*rddy
              p3=rddx*ddy
              p4=ddx*ddy
              do i=2,nz
                if (height(i).gt.ztra1(ipart)) then
                  indzm=i-1
                  indzp=i
                  goto 126
                endif
              end do
126           continue
              dz1=ztra1(ipart)-height(indzm)
              dz2=height(indzp)-ztra1(ipart)
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

              if (((ztra1(ipart).gt.3000.).and. &
                   (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                nclass(ipart)=min(int(ran1(idummy)* &
                     real(nclassunc))+1,nclassunc)
                numactiveparticles=numactiveparticles+1
                numparticlecount=numparticlecount+1
                npoint(ipart)=numparticlecount
                idt(ipart)=mintime
                itra1(ipart)=itime
                itramem(ipart)=itra1(ipart)
                itrasplit(ipart)=itra1(ipart)+ldirect*itsplit
                xmass1(ipart,1)=xmassperparticle
                if (mdomainfill.eq.2) xmass1(ipart,1)= &
                     xmass1(ipart,1)*pvpart*48./29.*ozonescale/10.**9
              else
                goto 171
              endif


  ! Increase numpart, if necessary
  !*******************************
              numpart=max(numpart,ipart)
              goto 173      ! Storage space has been found, stop searching
            endif
          end do
          if (ipart.gt.maxpart_mpi) &
               stop 'boundcond_domainfill.f: too many particles required'
173       minpart=ipart+1
171       continue
        end do


      end do
    end do
  end do


  xm=0.
  do i=1,numpart
    if (itra1(i).eq.itime) xm=xm+xmass1(i,1)
  end do

  !write(*,*) itime,numactiveparticles,numparticlecount,numpart,
  !    +xm,accmasst,xm+accmasst


  ! If particles shall be dumped, then accumulated masses at the domain boundaries
  ! must be dumped, too, to be used for later runs
  !*****************************************************************************

! :TODO: eso parallelize
  if ((ipout.gt.0).and.(itime.eq.loutend)) then
    if (lroot) then
      open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
           form='unformatted')
      write(unitboundcond) numcolumn_we,numcolumn_sn, &
           zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
      close(unitboundcond)
    end if
  endif

end subroutine boundcond_domainfill
