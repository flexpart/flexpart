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

program flexpart

  !*****************************************************************************
  !                                                                            *
  !     This is the Lagrangian Particle Dispersion Model FLEXPART.             *
  !     The main program manages the reading of model run specifications, etc. *
  !     All actual computing is done within subroutine timemanager.            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
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
  character(len=256) :: inline_options  !pathfile, flexversion, arg2
  integer :: index_v

  ! Generate a large number of random numbers
  !******************************************

  do i=1,maxrand-1,2
    call gasdev1(idummy,rannumb(i),rannumb(i+1))
  end do
  call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))

  ! FLEXPART version string
  ! flexversion='Version 9.2 beta (2014-07-01)'
  !flexversion='Version 9.2.0.1 (2015-01-27)'
  flexversion='Version 9.2.0.2 (2015-03-01)'
  ! default inlide options
  inline_options='none'
  !verbosity flags  defined in com_mod.f90
 
  ! Read the pathnames where input/output files are stored
  !*******************************************************

  select case (iargc())
  case (2) !2 parameters: pathfile and inline options
    call getarg(1,arg1)
    pathfile=arg1
    call getarg(2,arg2)
    inline_options=arg2
  case (1) !1 parameter pathfiel or inline options
    call getarg(1,arg1)
    pathfile=arg1
    if (arg1(1:1).eq.'-') then
      write(pathfile,'(a11)') './pathnames'
      inline_options=arg1 
    endif
  case (0) !default behavior
    write(pathfile,'(a11)') './pathnames'
  end select
  
  ! Print the GPL License statement
  !*******************************************************
  print*,'Welcome to FLEXPART ', trim(flexversion)
  print*,'FLEXPART is free software released under the GNU General Public License.'

  ! inline options allow to fine tune the verbosity during run time
  ! e.g.: show compilation parameters or input variables, time execution     
  if (inline_options(1:1).eq.'-') then
   ! if (index(inline_options,'v').gt.0) then
   !    print*, 'verbose mode'
   !    verbosity=1
   !    index_v=index(inline_options,'v')
   !    if (inline_options(index_v+1:index_v+1).eq.'2') then
   !    verbosity=2
   !    endif         
   ! endif   
 
    !if (trim(inline_options).eq.'-v'.or.trim(inline_options).eq.'-v1') then
    if (index(inline_options,'v').gt.0) then
       index_v=index(inline_options,'v')
       print*, 'Verbose mode: display  additional information during run'
       verbosity=1
       if (inline_options(index_v+1:index_v+1).eq.'2') then
       verbosity=2
       endif
       print*, 'verbosity level=', verbosity !inline_options(index_v+1:index_v+1)
    endif
    !iif (trim(inline_options).eq.'-v2') then
    !   print*, 'Verbose mode 2: display more detailed information during run'
    !   verbosity=2
    !endif

    if (index(inline_options,'i').gt.0) then   
       index_v=index(inline_options,'i')
       print*, 'Info mode: provide compile and run specific information, then stop'
       verbosity=1
       info_flag=1
       if (inline_options(index_v+1:index_v+1).eq.'2') then
       info_flag=2
       endif
    endif
    if (index(inline_options,'t').gt.0) then
       time_flag=1
       print*, 'timing execution activated'
       !stop
    endif
    if (index(inline_options,'d').gt.0) then
       debug_flag=1
       print*, 'debug messages activated'
       print*, 'debug_flag=', debug_flag
       !these messages give additional info on top on verbose mode
    endif 
  endif
           
  if (verbosity.gt.0) then
    print*, 'FLEXPART>******************************'
    print*, 'FLEXPART>* verbosity level:', verbosity
    print*, 'FLEXPART>* info only:      ', info_flag
    print*, 'FLEXPART>* time execution: ', time_flag
    print*, 'FLEXPART>******************************'
    
    print*, 'FLEXPART> parameters from par_mod'    
    print*, 'FLEXPART> nxmax=  ', nxmax
    print*, 'FLEXPART> nymax=  ', nymax
    print*, 'FLEXPART> nuvzmax=', nuvzmax
    print*, 'FLEXPART> nwzmax= ', nwzmax
    print*, 'FLEXPART> nzmax=  ', nzmax
    print*, 'FLEXPART> nxshift=', nxshift
    print*, 'FLEXPART> maxpart=', maxpart
    print*, 'FLEXPART> maxspec=', maxspec 

    if (info_flag.eq.1) stop
    write(*,*) 'call readpaths'
  endif 
  call readpaths(pathfile)
 
  !if (time_flag.gt.1) then !show clock info 
     !count=0,count_rate=1000
  CALL SYSTEM_CLOCK(count_clock0, count_rate, count_max)
     !WRITE(*,*) 'SYSTEM_CLOCK',count, count_rate, count_max
     !WRITE(*,*) 'SYSTEM_CLOCK, count_clock0', count_clock0
     !WRITE(*,*) 'SYSTEM_CLOCK, count_rate', count_rate
     !WRITE(*,*) 'SYSTEM_CLOCK, count_max', count_max
  !endif

  ! Read the user specifications for the current model run
  !*******************************************************

  if (verbosity.gt.0) then
    write(*,*) 'FLEXPART> call readcommand'
  endif
  call readcommand
  if (verbosity.gt.0) then
    write(*,*) '    ldirect      =', ldirect
    write(*,*) '    ibdate,ibtime=', ibdate,ibtime
    write(*,*) '    iedate,ietime=', iedate,ietime
  endif
    if (time_flag.gt.0) then   
      CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
      write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
    endif     

  ! Read the age classes to be used
  !********************************
  if (verbosity.gt.0) then
    write(*,*) 'FLEXPART> call readageclasses'
  endif
  call readageclasses

  if (time_flag.gt.1) then   
    CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
    write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
  endif     

  ! Read, which wind fields are available within the modelling period
  !******************************************************************

  if (verbosity.gt.0) then
    write(*,*) 'FLEXPART> call readavailable'
  endif  
  call readavailable

  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************
 
  if (verbosity.gt.0) then
     write(*,*) 'FLEXPART> call gridcheck'
  endif
  call gridcheck

  if (time_flag.gt.0) then   
    CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
    write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
  endif      

  if (verbosity.gt.0) then
    write(*,*) 'FLEXPART> call gridcheck_nests'
  endif  
  call gridcheck_nests

  ! Read the output grid specifications
  !************************************

  if (verbosity.gt.0) then
    write(*,*) 'FLEXPART> call readoutgrid'
  endif

  call readoutgrid

  if (nested_output.eq.1) then
    call readoutgrid_nest
    if (verbosity.gt.0) then
      write(*,*) 'FLEXPART>  readoutgrid_nest'
    endif
  endif

  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************

  if (verbosity.eq.1) then
     print*,'FLEXPART> call readreceptors'
  endif
  call readreceptors

  ! Read the physico-chemical species property table
  !*************************************************
  !SEC: now only needed SPECIES are read in readreleases.f
  !call readspecies


  ! Read the landuse inventory
  !***************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call readlanduse'
  endif
  call readlanduse

  ! Assign fractional cover of landuse classes to each ECMWF grid point
  !********************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call assignland'
  endif
  call assignland

  ! Read the coordinates of the release locations
  !**********************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call readreleases'
  endif
  call readreleases

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call readdepo'
  endif
  call readdepo

  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************

  call coordtrafo  
  if (verbosity.gt.0) then
    print*,'FLEXPART> call coordtrafo'
  endif

  ! Initialize all particles to non-existent
  !*****************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> Initialize all particles to non-existent'
  endif
  do j=1,maxpart
    itra1(j)=-999999999
  end do

  ! For continuation of previous run, read in particle positions
  !*************************************************************

  if (ipin.eq.1) then
    if (verbosity.gt.0) then
      print*,'FLEXPART> call readpartpositions'
    endif
    call readpartpositions
  else
    if (verbosity.gt.0) then
      print*,'FLEXPART> set numpart=0, numparticlecount=0'
    endif    
    numpart=0
    numparticlecount=0
  endif

  ! Calculate volume, surface area, etc., of all output grid cells
  ! Allocate fluxes and OHfield if necessary
  !***************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call outgrid_init'
  endif
  call outgrid_init
  if (nested_output.eq.1) call outgrid_init_nest

  ! Read the OH field
  !******************

  if (OHREA.eqv..TRUE.) then
    if (verbosity.gt.0) then
      print*,'FLEXPART> call readOHfield'
    endif
    call readOHfield
  endif

  ! Write basic information on the simulation to a file "header"
  ! and open files that are to be kept open throughout the simulation
  !******************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> call variuos writeheader routines'
  endif

  call writeheader
  ! write header in ASCII format 
  call writeheader_txt
  !if (nested_output.eq.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.ne.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.eq.1) call writeheader_nest_surf
  if (nested_output.ne.1.and.surf_only.eq.1) call writeheader_surf

  !open(unitdates,file=path(2)(1:length(2))//'dates')

  if (verbosity.gt.0) then
    print*,'FLEXPART> call openreceptors'
  endif
  call openreceptors
  if ((iout.eq.4).or.(iout.eq.5)) call openouttraj

  ! Releases can only start and end at discrete times (multiples of lsynctime)
  !***************************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> discretize release times'
  endif
  do i=1,numpoint
    ireleasestart(i)=nint(real(ireleasestart(i))/real(lsynctime))*lsynctime
    ireleaseend(i)=nint(real(ireleaseend(i))/real(lsynctime))*lsynctime
  end do

  ! Initialize cloud-base mass fluxes for the convection scheme
  !************************************************************

  if (verbosity.gt.0) then
    print*,'FLEXPART> Initialize cloud-base mass fluxes for the convection scheme'
  endif

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

  if (time_flag.gt.0) then   
    CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
    write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
  endif
  if (info_flag.eq.2) then
    print*, 'FLEXPART> info only mode (stop before call timemanager)'
    stop
  endif
  if (verbosity.gt.0) then
     print*,'FLEXPART> call timemanager'
  endif

  call timemanager

  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLEXPART MODEL RUN!'

end program flexpart
