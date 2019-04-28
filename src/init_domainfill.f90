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

subroutine init_domainfill
!
!*****************************************************************************
!                                                                            *
! Initializes particles equally distributed over the first release location  *
! specified in file RELEASES. This box is assumed to be the domain for doing *
! domain-filling trajectory calculations.                                    *
! All particles carry the same amount of mass which alltogether comprises the*
! mass of air within the box.                                                *
!                                                                            *
!     Author: A. Stohl                                                       *
!                                                                            *
!     15 October 2002                                                        *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
!                                                                            *
! numparticlecount    consecutively counts the number of particles released  *
! nx_we(2)       grid indices for western and eastern boundary of domain-    *
!                filling trajectory calculations                             *
! ny_sn(2)       grid indices for southern and northern boundary of domain-  *
!                filling trajectory calculations                             *
!                                                                            *
!*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use random_mod

  implicit none

  integer :: j,ix,jy,kz,ncolumn,numparttot
  real :: gridarea(0:nymax-1),pp(nzmax),ylat,ylatp,ylatm,hzone
  real :: cosfactm,cosfactp,deltacol,dz1,dz2,dz,pnew,fractus
  real,parameter :: pih=pi/180.
  real :: colmass(0:nxmax-1,0:nymax-1),colmasstotal,zposition

  integer :: ixm,ixp,jym,jyp,indzm,indzp,in,indzh,i,jj
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2)

  integer :: idummy = -11


! Determine the release region (only full grid cells), over which particles
! shall be initialized
! Use 2 fields for west/east and south/north boundary
!**************************************************************************

  nx_we(1)=max(int(xpoint1(1)),0)
  nx_we(2)=min((int(xpoint2(1))+1),nxmin1)
  ny_sn(1)=max(int(ypoint1(1)),0)
  ny_sn(2)=min((int(ypoint2(1))+1),nymin1)

! For global simulations (both global wind data and global domain-filling),
! set a switch, such that no boundary conditions are used
!**************************************************************************
  if (xglobal.and.sglobal.and.nglobal) then
    if ((nx_we(1).eq.0).and.(nx_we(2).eq.nxmin1).and. &
         (ny_sn(1).eq.0).and.(ny_sn(2).eq.nymin1)) then
      gdomainfill=.true.
    else
      gdomainfill=.false.
    endif
  endif

! Do not release particles twice (i.e., not at both in the leftmost and rightmost
! grid cell) for a global domain
!*****************************************************************************
  if (xglobal) nx_we(2)=min(nx_we(2),nx-2)


! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
!************************************************************

  do jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    ylat=ylat0+real(jy)*dy
    ylatp=ylat+0.5*dy
    ylatm=ylat-0.5*dy
    if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
      hzone=1./dyconst
    else
      cosfactp=cos(ylatp*pih)*r_earth
      cosfactm=cos(ylatm*pih)*r_earth
      if (cosfactp.lt.cosfactm) then
        hzone=sqrt(r_earth**2-cosfactp**2)- &
             sqrt(r_earth**2-cosfactm**2)
      else
        hzone=sqrt(r_earth**2-cosfactm**2)- &
             sqrt(r_earth**2-cosfactp**2)
      endif
    endif
    gridarea(jy)=2.*pi*r_earth*hzone*dx/360.
  end do

! Do the same for the south pole

  if (sglobal) then
    ylat=ylat0
    ylatp=ylat+0.5*dy
    ylatm=ylat
    cosfactm=0.
    cosfactp=cos(ylatp*pih)*r_earth
    hzone=sqrt(r_earth**2-cosfactm**2)- &
         sqrt(r_earth**2-cosfactp**2)
    gridarea(0)=2.*pi*r_earth*hzone*dx/360.
  endif

! Do the same for the north pole

  if (nglobal) then
    ylat=ylat0+real(nymin1)*dy
    ylatp=ylat
    ylatm=ylat-0.5*dy
    cosfactp=0.
    cosfactm=cos(ylatm*pih)*r_earth
    hzone=sqrt(r_earth**2-cosfactp**2)- &
         sqrt(r_earth**2-cosfactm**2)
    gridarea(nymin1)=2.*pi*r_earth*hzone*dx/360.
  endif


! Calculate total mass of each grid column and of the whole atmosphere
!*********************************************************************

  colmasstotal=0.
  do jy=ny_sn(1),ny_sn(2)          ! loop about latitudes
    do ix=nx_we(1),nx_we(2)      ! loop about longitudes
      pp(1)=rho(ix,jy,1,1)*r_air*tt(ix,jy,1,1)
      pp(nz)=rho(ix,jy,nz,1)*r_air*tt(ix,jy,nz,1)
      colmass(ix,jy)=(pp(1)-pp(nz))/ga*gridarea(jy)
      colmasstotal=colmasstotal+colmass(ix,jy)
    end do
  end do

  write(*,*) 'Atm. mass: ',colmasstotal


  if (ipin.eq.0) numpart=0

! Determine the particle positions
!*********************************

  numparttot=0
  numcolumn=0
  do jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    ylat=ylat0+real(jy)*dy
    do ix=nx_we(1),nx_we(2)      ! loop about longitudes
      ncolumn=nint(0.999*real(npart(1))*colmass(ix,jy)/ &
           colmasstotal)
      if (ncolumn.eq.0) goto 30
      if (ncolumn.gt.numcolumn) numcolumn=ncolumn

! Calculate pressure at the altitudes of model surfaces, using the air density
! information, which is stored as a 3-d field
!*****************************************************************************

      do kz=1,nz
        pp(kz)=rho(ix,jy,kz,1)*r_air*tt(ix,jy,kz,1)
      end do


      deltacol=(pp(1)-pp(nz))/real(ncolumn)
      pnew=pp(1)+deltacol/2.
      jj=0
      do j=1,ncolumn
        jj=jj+1


! For columns with many particles (i.e. around the equator), distribute
! the particles equally, for columns with few particles (i.e. around the
! poles), distribute the particles randomly
!***********************************************************************


        if (ncolumn.gt.20) then
          pnew=pnew-deltacol
        else
          pnew=pp(1)-ran1(idummy)*(pp(1)-pp(nz))
        endif

        do kz=1,nz-1
          if ((pp(kz).ge.pnew).and.(pp(kz+1).lt.pnew)) then
            dz1=pp(kz)-pnew
            dz2=pnew-pp(kz+1)
            dz=1./(dz1+dz2)

! Assign particle position
!*************************
! Do the following steps only if particles are not read in from previous model run
!*****************************************************************************
            if (ipin.eq.0) then
              xtra1(numpart+jj)=real(ix)-0.5+ran1(idummy)
              if (ix.eq.0) xtra1(numpart+jj)=ran1(idummy)
              if (ix.eq.nxmin1) xtra1(numpart+jj)= &
                   real(nxmin1)-ran1(idummy)
              ytra1(numpart+jj)=real(jy)-0.5+ran1(idummy)
              ztra1(numpart+jj)=(height(kz)*dz2+height(kz+1)*dz1)*dz
              if (ztra1(numpart+jj).gt.height(nz)-0.5) &
                   ztra1(numpart+jj)=height(nz)-0.5


! Interpolate PV to the particle position
!****************************************
              ixm=int(xtra1(numpart+jj))
              jym=int(ytra1(numpart+jj))
              ixp=ixm+1
              jyp=jym+1
              ddx=xtra1(numpart+jj)-real(ixm)
              ddy=ytra1(numpart+jj)-real(jym)
              rddx=1.-ddx
              rddy=1.-ddy
              p1=rddx*rddy
              p2=ddx*rddy
              p3=rddx*ddy
              p4=ddx*ddy
              do i=2,nz
                if (height(i).gt.ztra1(numpart+jj)) then
                  indzm=i-1
                  indzp=i
                  goto 6
                endif
              end do
6             continue
              dz1=ztra1(numpart+jj)-height(indzm)
              dz2=height(indzp)-ztra1(numpart+jj)
              dz=1./(dz1+dz2)
              do in=1,2
                indzh=indzm+in-1
                y1(in)=p1*pv(ixm,jym,indzh,1) &
                     +p2*pv(ixp,jym,indzh,1) &
                     +p3*pv(ixm,jyp,indzh,1) &
                     +p4*pv(ixp,jyp,indzh,1)
              end do
              pvpart=(dz2*y1(1)+dz1*y1(2))*dz
              if (ylat.lt.0.) pvpart=-1.*pvpart


! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
!*****************************************************************************

              if (((ztra1(numpart+jj).gt.3000.).and. &
                   (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then

! Assign certain properties to the particle
!******************************************
                nclass(numpart+jj)=min(int(ran1(idummy)* &
                     real(nclassunc))+1,nclassunc)
                numparticlecount=numparticlecount+1
                npoint(numpart+jj)=numparticlecount
                idt(numpart+jj)=mintime
                itra1(numpart+jj)=0
                itramem(numpart+jj)=0
                itrasplit(numpart+jj)=itra1(numpart+jj)+ldirect* &
                     itsplit
                xmass1(numpart+jj,1)=colmass(ix,jy)/real(ncolumn)
                if (mdomainfill.eq.2) xmass1(numpart+jj,1)= &
                     xmass1(numpart+jj,1)*pvpart*48./29.*ozonescale/10.**9
              else
                jj=jj-1
              endif
            endif
          endif
        end do
      end do
      numparttot=numparttot+ncolumn
      if (ipin.eq.0) numpart=numpart+jj
30    continue
    end do
  end do


! Check whether numpart is really smaller than maxpart
!*****************************************************

! ESO :TODO: this warning need to be moved further up, else out-of-bounds error earlier
  if (numpart.gt.maxpart) then
    write(*,*) 'numpart too large: change source in init_atm_mass.f'
    write(*,*) 'numpart: ',numpart,' maxpart: ',maxpart
  endif


  xmassperparticle=colmasstotal/real(numparttot)


! Make sure that all particles are within domain
!***********************************************

  do j=1,numpart
    if ((xtra1(j).lt.0.).or.(xtra1(j).ge.real(nxmin1)).or. &
         (ytra1(j).lt.0.).or.(ytra1(j).ge.real(nymin1))) then
      itra1(j)=-999999999
    endif
  end do




! For boundary conditions, we need fewer particle release heights per column,
! because otherwise it takes too long until enough mass has accumulated to
! release a particle at the boundary (would take dx/u seconds), leading to
! relatively large position errors of the order of one grid distance.
! It's better to release fewer particles per column, but to do so more often.
! Thus, use on the order of nz starting heights per column.
! We thus repeat the above to determine fewer starting heights, that are
! used furtheron in subroutine boundcond_domainfill.f.
!****************************************************************************

  fractus=real(numcolumn)/real(nz)
  write(*,*) 'Total number of particles at model start: ',numpart
  write(*,*) 'Maximum number of particles per column: ',numcolumn
  write(*,*) 'If ',fractus,' <1, better use more particles'
  fractus=sqrt(max(fractus,1.))/2.

  do jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
    do ix=nx_we(1),nx_we(2)      ! loop about longitudes
      ncolumn=nint(0.999/fractus*real(npart(1))*colmass(ix,jy) &
           /colmasstotal)
      if (ncolumn.gt.maxcolumn) stop 'maxcolumn too small'
      if (ncolumn.eq.0) goto 80


! Memorize how many particles per column shall be used for all boundaries
! This is further used in subroutine boundcond_domainfill.f
! Use 2 fields for west/east and south/north boundary
!************************************************************************

      if (ix.eq.nx_we(1)) numcolumn_we(1,jy)=ncolumn
      if (ix.eq.nx_we(2)) numcolumn_we(2,jy)=ncolumn
      if (jy.eq.ny_sn(1)) numcolumn_sn(1,ix)=ncolumn
      if (jy.eq.ny_sn(2)) numcolumn_sn(2,ix)=ncolumn

! Calculate pressure at the altitudes of model surfaces, using the air density
! information, which is stored as a 3-d field
!*****************************************************************************

      do kz=1,nz
        pp(kz)=rho(ix,jy,kz,1)*r_air*tt(ix,jy,kz,1)
      end do

! Determine the reference starting altitudes
!*******************************************

      deltacol=(pp(1)-pp(nz))/real(ncolumn)
      pnew=pp(1)+deltacol/2.
      do j=1,ncolumn
        pnew=pnew-deltacol
        do kz=1,nz-1
          if ((pp(kz).ge.pnew).and.(pp(kz+1).lt.pnew)) then
            dz1=pp(kz)-pnew
            dz2=pnew-pp(kz+1)
            dz=1./(dz1+dz2)
            zposition=(height(kz)*dz2+height(kz+1)*dz1)*dz
            if (zposition.gt.height(nz)-0.5) zposition=height(nz)-0.5

! Memorize vertical positions where particles are introduced
! This is further used in subroutine boundcond_domainfill.f
!***********************************************************

            if (ix.eq.nx_we(1)) zcolumn_we(1,jy,j)=zposition
            if (ix.eq.nx_we(2)) zcolumn_we(2,jy,j)=zposition
            if (jy.eq.ny_sn(1)) zcolumn_sn(1,ix,j)=zposition
            if (jy.eq.ny_sn(2)) zcolumn_sn(2,ix,j)=zposition

! Initialize mass that has accumulated at boundary to zero
!*********************************************************

            acc_mass_we(1,jy,j)=0.
            acc_mass_we(2,jy,j)=0.
            acc_mass_sn(1,jy,j)=0.
            acc_mass_sn(2,jy,j)=0.
          endif
        end do
      end do
80    continue
    end do
  end do

! Reduce numpart if invalid particles at end of arrays
  do i=numpart, 1, -1
      if (itra1(i).eq.-999999999) then
        numpart=numpart-1
      else
        exit
      end if
    end do

! If particles shall be read in to continue an existing run,
! then the accumulated masses at the domain boundaries must be read in, too.
! This overrides any previous calculations.
!***************************************************************************

  if ((ipin.eq.1).and.(.not.gdomainfill)) then
    open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
         form='unformatted')
    read(unitboundcond) numcolumn_we,numcolumn_sn, &
         zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
    close(unitboundcond)
  endif


end subroutine init_domainfill
