!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann, Ignacio Pisso                                      *
!                                                                     *
! This file is part of FLEXPART-NorESM                                *
!                                                                     *
! FLEXPART-NorESM is free software: you can redistribute it           *
! and/or modify                                                       *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART-NorESM is distributed in the hope that it will be useful,  *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART-NorESM.                                         *
!  If not, see <http://www.gnu.org/licenses/>.                        * 
!**********************************************************************

program flexpart_NorESM

  !*****************************************************************************
  !                                                                            *
  !     This is the Lagrangian Particle Dispersion Model FLEXPART-NorESM.      *
  !     based on FLEXPART version 9.0.1                                        *
  !     The main program manages the reading of model run specifications, etc. *
  !     All actual computing is done within subroutine timemanager.            *
  !     Note: FLEXPRT-NorESM uses meteo files from the NorESM(CAM4 based)      *
  !     climate model in netcdf format.                                        * 
  !     Note: some lines of code to check netcdf files content have been       *
  !     copied/modfied from routines in FLEXPART-WRF                           *
  !                                                                            *
  !     Author:                                                                *
  !     A. Stohl 18 May 1996                                                   *   
  !                                                                            *
  !     Modified by:                                                           *
  !     M. Cassiani  2016                                                      *
  !     Nested input not allowed                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use conv_mod

  implicit none

  integer :: i,j,ix,jy,inest
  integer :: idummy = -320

  ! Generate a large number of random numbers
  !******************************************

  do i=1,maxrand-1,2
    call gasdev1(idummy,rannumb(i),rannumb(i+1))
  end do
  call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))

  ! Print the GPL License statement
  !*******************************************************
  print*,'Welcome to FLEXPART-NorESM Version 1.0'
  print*,'FLEXPART-NorESM is free software released under the GNU'// &
       'General Public License.'

  ! Read the pathnames where input/output files are stored
  !*******************************************************

  call readpaths

  ! Read the user specifications for the current model run
  !*******************************************************

  call readcommand


  ! Read the age classes to be used
  !********************************

  call readageclasses


  ! Read, which wind fields are available within the modelling period
  !******************************************************************

  call readavailable


  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************

  call gridcheck
  !call gridcheck_nests ! NORESM VERSION NESTING OF INPUT NOT ACTIVATED : comment by  mc


  ! Read the output grid specifications
  !************************************

  call readoutgrid
  if (nested_output.eq.1) call readoutgrid_nest


  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************

  call readreceptors


  ! Read the physico-chemical species property table
  !*************************************************
  !SEC: now only needed SPECIES are read in readreleases.f
  !call readspecies


  ! Read the landuse inventory
  !***************************

  call readlanduse


  ! Assign fractional cover of landuse classes to each grid point
  !********************************************************************

  call assignland



  ! Read the coordinates of the release locations
  !**********************************************

  call readreleases

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************

  call readdepo


  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************

  call coordtrafo


  ! Initialize all particles to non-existent
  !*****************************************

  do j=1,maxpart
    itra1(j)=-999999999
  end do

  ! For continuation of previous run, read in particle positions
  !*************************************************************

  if (ipin.eq.1) then
    call readpartpositions
  else
    numpart=0
    numparticlecount=0
  endif


  ! Calculate volume, surface area, etc., of all output grid cells
  ! Allocate fluxes and OHfield if necessary
  !***************************************************************

  call outgrid_init
  if (nested_output.eq.1) call outgrid_init_nest


  ! Read the OH field
  !******************

  if (OHREA.eqv..TRUE.) &
       call readOHfield

  ! Write basic information on the simulation to a file "header"
  ! and open files that are to be kept open throughout the simulation
  !******************************************************************

  call writeheader
  if (nested_output.eq.1) call writeheader_nest
  open(unitdates,file=path(2)(1:length(2))//'dates')
  call openreceptors
  if ((iout.eq.4).or.(iout.eq.5)) call openouttraj


  ! Releases can only start and end at discrete times (multiples of lsynctime)
  !***************************************************************************

  do i=1,numpoint
    ireleasestart(i)=nint(real(ireleasestart(i))/ &
         real(lsynctime))*lsynctime
    ireleaseend(i)=nint(real(ireleaseend(i))/ &
         real(lsynctime))*lsynctime
  end do


  ! Initialize cloud-base mass fluxes for the convection scheme
  !************************************************************

  do jy=0,nymin1
    do ix=0,nxmin1
      cbaseflux(ix,jy)=0.
    end do
  end do
  do inest=1,numbnests
    do jy=0,nyn(inest)-1
      do ix=0,nxn(inest)-1
        cbasefluxn(ix,jy,inest)=0.
      end do
    end do
  end do


  ! Calculate particle trajectories
  !********************************

  call timemanager


  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLE&
       &XPART MODEL RUN!'

end program flexpart_NorESM
