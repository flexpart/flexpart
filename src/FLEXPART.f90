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
  ! Changes:                                                                   *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Added detection of metdata format using gributils routines           *
  !     - Distinguished calls to ecmwf/gfs gridcheck versions based on         *
  !       detected metdata format                                              *
  !     - Passed metdata format down to timemanager                            *
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
  use netcdf_output_mod, only: writeheader_netcdf
  use random_mod, only: gasdev1
  use class_gribfile

  implicit none

  integer :: i,j,ix,jy,inest
  integer :: idummy = -320
  character(len=256) :: inline_options  !pathfile, flexversion, arg2
  integer :: metdata_format = GRIBFILE_CENTRE_UNKNOWN
  integer :: detectformat



  ! Initialize arrays in com_mod
  !*****************************
  call com_mod_allocate_part(maxpart)
  

  ! Generate a large number of random numbers
  !******************************************

  do i=1,maxrand-1,2
    call gasdev1(idummy,rannumb(i),rannumb(i+1))
  end do
  call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))


  ! FLEXPART version string
  flexversion_major = '10' ! Major version number, also used for species file names
  flexversion='Version '//trim(flexversion_major)//'.2beta (2017-08-01)'
  verbosity=0

  ! Read the pathnames where input/output files are stored
  !*******************************************************

  inline_options='none'
  select case (iargc())
  case (2)
    call getarg(1,arg1)
    pathfile=arg1
    call getarg(2,arg2)
    inline_options=arg2
  case (1)
    call getarg(1,arg1)
    pathfile=arg1
    if (arg1(1:1).eq.'-') then
      write(pathfile,'(a11)') './pathnames'
      inline_options=arg1 
    endif
  case (0)
    write(pathfile,'(a11)') './pathnames'
  end select
  
  ! Print the GPL License statement
  !*******************************************************
  print*,'Welcome to FLEXPART ', trim(flexversion)
  print*,'FLEXPART is free software released under the GNU General Public License.'
 
  if (inline_options(1:1).eq.'-') then
    if (trim(inline_options).eq.'-v'.or.trim(inline_options).eq.'-v1') then
       print*, 'Verbose mode 1: display detailed information during run'
       verbosity=1
    endif
    if (trim(inline_options).eq.'-v2') then
       print*, 'Verbose mode 2: display more detailed information during run'
       verbosity=2
    endif
    if (trim(inline_options).eq.'-i') then
       print*, 'Info mode: provide detailed run specific information and stop'
       verbosity=1
       info_flag=1
    endif
    if (trim(inline_options).eq.'-i2') then
       print*, 'Info mode: provide more detailed run specific information and stop'
       verbosity=2
       info_flag=1
    endif
  endif
           
  if (verbosity.gt.0) then
    write(*,*) 'call readpaths'
  endif 
  call readpaths(pathfile)
 
  if (verbosity.gt.1) then !show clock info 
     !print*,'length(4)',length(4)
     !count=0,count_rate=1000
     CALL SYSTEM_CLOCK(count_clock0, count_rate, count_max)
     !WRITE(*,*) 'SYSTEM_CLOCK',count, count_rate, count_max
     !WRITE(*,*) 'SYSTEM_CLOCK, count_clock0', count_clock0
     !WRITE(*,*) 'SYSTEM_CLOCK, count_rate', count_rate
     !WRITE(*,*) 'SYSTEM_CLOCK, count_max', count_max
  endif

  ! Read the user specifications for the current model run
  !*******************************************************

  if (verbosity.gt.0) then
    write(*,*) 'call readcommand'
  endif
  call readcommand
  if (verbosity.gt.0) then
    write(*,*) '    ldirect=', ldirect
    write(*,*) '    ibdate,ibtime=',ibdate,ibtime
    write(*,*) '    iedate,ietime=', iedate,ietime
    if (verbosity.gt.1) then   
      CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
      write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
    endif     
  endif

  ! Read the age classes to be used
  !********************************
  if (verbosity.gt.0) then
    write(*,*) 'call readageclasses'
  endif
  call readageclasses

  if (verbosity.gt.1) then   
    CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
    write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
  endif     

  ! Read, which wind fields are available within the modelling period
  !******************************************************************

  if (verbosity.gt.0) then
    write(*,*) 'call readavailable'
  endif  
  call readavailable

  ! Detect metdata format
  !**********************

  metdata_format = detectformat()

  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    print *,'ECMWF metdata detected'
  elseif (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
    print *,'NCEP metdata detected'
  else
    print *,'Unknown metdata format'
    return
  endif



  ! If nested wind fields are used, allocate arrays
  !************************************************

  if (verbosity.gt.0) then
    write(*,*) 'call com_mod_allocate_nests'
  endif
  call com_mod_allocate_nests

  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************
 
  if (verbosity.gt.0) then
     write(*,*) 'call gridcheck'
  endif

  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    call gridcheck_ecmwf
  else
    call gridcheck_gfs
  end if

  if (verbosity.gt.1) then   
    CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
    write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
  endif      

  if (verbosity.gt.0) then
    write(*,*) 'call gridcheck_nests'
  endif  
  call gridcheck_nests

  ! Read the output grid specifications
  !************************************

  if (verbosity.gt.0) then
    write(*,*) 'call readoutgrid'
  endif

  call readoutgrid

  if (nested_output.eq.1) then
    call readoutgrid_nest
    if (verbosity.gt.0) then
      write(*,*) '# readoutgrid_nest'
    endif
  endif

  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************

  if (verbosity.eq.1) then
     print*,'call readreceptors'
  endif
  call readreceptors

  ! Read the physico-chemical species property table
  !*************************************************
  !SEC: now only needed SPECIES are read in readreleases.f
  !call readspecies


  ! Read the landuse inventory
  !***************************

  if (verbosity.gt.0) then
    print*,'call readlanduse'
  endif
  call readlanduse

  ! Assign fractional cover of landuse classes to each ECMWF grid point
  !********************************************************************

  if (verbosity.gt.0) then
    print*,'call assignland'
  endif
  call assignland

  ! Read the coordinates of the release locations
  !**********************************************

  if (verbosity.gt.0) then
    print*,'call readreleases'
  endif
  call readreleases

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************

  if (verbosity.gt.0) then
    print*,'call readdepo'
  endif
  call readdepo

  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************

  call coordtrafo  
  if (verbosity.gt.0) then
    print*,'call coordtrafo'
  endif

  ! Initialize all particles to non-existent
  !*****************************************

  if (verbosity.gt.0) then
    print*,'Initialize all particles to non-existent'
  endif
  do j=1,maxpart
    itra1(j)=-999999999
  end do

  ! For continuation of previous run, read in particle positions
  !*************************************************************

  if (ipin.eq.1) then
    if (verbosity.gt.0) then
      print*,'call readpartpositions'
    endif
    call readpartpositions
  else
    if (verbosity.gt.0) then
      print*,'numpart=0, numparticlecount=0'
    endif    
    numpart=0
    numparticlecount=0
  endif

  ! Calculate volume, surface area, etc., of all output grid cells
  ! Allocate fluxes and OHfield if necessary
  !***************************************************************

  if (verbosity.gt.0) then
    print*,'call outgrid_init'
  endif
  call outgrid_init
  if (nested_output.eq.1) call outgrid_init_nest

  ! Read the OH field
  !******************

  if (OHREA.eqv..TRUE.) then
    if (verbosity.gt.0) then
      print*,'call readOHfield'
    endif
    call readOHfield
  endif

  ! Write basic information on the simulation to a file "header"
  ! and open files that are to be kept open throughout the simulation
  !******************************************************************

  if (lnetcdfout.eq.1) then 
    call writeheader_netcdf(lnest=.false.)
  else 
    call writeheader
  end if

  if (nested_output.eq.1) then
    if (lnetcdfout.eq.1) then
      call writeheader_netcdf(lnest=.true.)
    else
      call writeheader_nest
    endif
  endif

  if (verbosity.gt.0) then
    print*,'call writeheader'
  endif

  call writeheader
  ! FLEXPART 9.2 ticket ?? write header in ASCII format 
  call writeheader_txt
  !if (nested_output.eq.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.ne.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.eq.1) call writeheader_nest_surf
  if (nested_output.ne.1.and.surf_only.eq.1) call writeheader_surf

  !open(unitdates,file=path(2)(1:length(2))//'dates')

  if (verbosity.gt.0) then
    print*,'call openreceptors'
  endif
  call openreceptors
  if ((iout.eq.4).or.(iout.eq.5)) call openouttraj

  ! Releases can only start and end at discrete times (multiples of lsynctime)
  !***************************************************************************

  if (verbosity.gt.0) then
    print*,'discretize release times'
  endif
  do i=1,numpoint
    ireleasestart(i)=nint(real(ireleasestart(i))/real(lsynctime))*lsynctime
    ireleaseend(i)=nint(real(ireleaseend(i))/real(lsynctime))*lsynctime
  end do

  ! Initialize cloud-base mass fluxes for the convection scheme
  !************************************************************

  if (verbosity.gt.0) then
    print*,'Initialize cloud-base mass fluxes for the convection scheme'
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

  ! Inform whether output kernel is used or not
  !*********************************************
  if (lroot) then
    if (.not.lusekerneloutput) then
      write(*,*) "Concentrations are calculated without using kernel"
    else
      write(*,*) "Concentrations are calculated using kernel"
    end if
  end if


  ! Calculate particle trajectories
  !********************************

  if (verbosity.gt.0) then
     if (verbosity.gt.1) then   
       CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
       write(*,*) 'SYSTEM_CLOCK',(count_clock - count_clock0)/real(count_rate) !, count_rate, count_max
     endif
     if (info_flag.eq.1) then
       print*, 'info only mode (stop)'    
       stop
     endif
     print*,'call timemanager'
  endif

  call timemanager(metdata_format)

! NIK 16.02.2005 
  do i=1,nspec
    write(*,*) '**********************************************'
    write(*,*) 'Scavenging statistics for species ', species(i), ':'
    write(*,*) 'Total number of occurences of below-cloud scavenging', &
         & tot_blc_count(i)
    write(*,*) 'Total number of occurences of in-cloud    scavenging', &
         & tot_inc_count(i)
    write(*,*) '**********************************************'
  end do
  
  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLE&
       &XPART MODEL RUN!'

end program flexpart
