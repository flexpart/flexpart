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

subroutine readreleases

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the release point specifications for the current    *
  !     model run. Several release points can be used at the same time.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !                                                                            *
  !     Update: 29 January 2001                                                *
  !     Release altitude can be either in magl or masl                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decay               decay constant of species                              *
  ! dquer [um]          mean particle diameters                                *
  ! dsigma              e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass*
  !                     are between 0.1*dquer and 10*dquer                     *
  ! ireleasestart, ireleaseend [s] starting time and ending time of each       *
  !                     release                                                *
  ! kindz               1: zpoint is in m agl, 2: zpoint is in m asl, 3: zpoint*
  !                     is in hPa                                              *
  ! npart               number of particles to be released                     *
  ! nspec               number of species to be released                       *
  ! density [kg/m3]     density of the particles                               *
  ! rm [s/m]            Mesophyll resistance                                   *
  ! species             name of species                                        *
  ! xmass               total mass of each species                             *
  ! xpoint1,ypoint1     geograf. coordinates of lower left corner of release   *
  !                     area                                                   *
  ! xpoint2,ypoint2     geograf. coordinates of upper right corner of release  *
  !                     area                                                   *
  ! weta, wetb          parameters to determine the wet scavenging coefficient *
  ! zpoint1,zpoint2     height range, over which release takes place           *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use xmass_mod
  use par_mod
  use com_mod

  implicit none

  integer :: numpartmax,i,j,id1,it1,id2,it2,specnum_rel,idum,stat
  real :: vsh(ni),fracth(ni),schmih(ni),releaserate,xdum,cun
  real(kind=dp) :: jul1,jul2,juldate
  character(len=50) :: line
  logical :: old

  !sec, read release to find how many releasepoints should be allocated

  open(unitreleases,file=path(1)(1:length(1))//'RELEASES',status='old', &
       err=999)

  ! Check the format of the RELEASES file (either in free format,
  ! or using a formatted mask)
  ! Use of formatted mask is assumed if line 10 contains the word 'DIRECTION'
  !**************************************************************************

  call skplin(12,unitreleases)
  read (unitreleases,901) line
901   format (a)
  if (index(line,'Total') .eq. 0) then
    old = .false.
  else
    old = .true.
  endif
  rewind(unitreleases)


  ! Skip first 11 lines (file header)
  !**********************************

  call skplin(11,unitreleases)


  read(unitreleases,*,err=998) nspec
  if (old) call skplin(2,unitreleases)
  do i=1,nspec
    read(unitreleases,*,err=998) specnum_rel
    if (old) call skplin(2,unitreleases)
  end do

  numpoint=0
100   numpoint=numpoint+1
  read(unitreleases,*,end=25)
  read(unitreleases,*,err=998,end=25) idum,idum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) idum,idum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) idum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xdum
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) idum
  if (old) call skplin(2,unitreleases)
  do i=1,nspec
    read(unitreleases,*,err=998) xdum
    if (old) call skplin(2,unitreleases)
  end do
  !save compoint only for the first 1000 release points
  read(unitreleases,'(a40)',err=998) compoint(1)(1:40)
  if (old) call skplin(1,unitreleases)

  goto 100

25   numpoint=numpoint-1

  !allocate memory for numpoint releaspoint
    allocate(ireleasestart(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(ireleaseend(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(xpoint1(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(xpoint2(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(ypoint1(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(ypoint2(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(zpoint1(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(zpoint2(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(kindz(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(xmass(numpoint,maxspec) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(rho_rel(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(npart(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'
    allocate(xmasssave(numpoint) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate RELEASPOINT'

   write (*,*) ' Releasepoints allocated: ', numpoint

   do i=1,numpoint
     xmasssave(i)=0.
   end do

  !now save the information
  DEP=.false.
  DRYDEP=.false.
  WETDEP=.false.
  OHREA=.false.
  do i=1,maxspec
    DRYDEPSPEC(i)=.false.
  end do

  rewind(unitreleases)


  ! Skip first 11 lines (file header)
  !**********************************

  call skplin(11,unitreleases)


  ! Assign species-specific parameters needed for physical processes
  !*************************************************************************

  read(unitreleases,*,err=998) nspec
  if (nspec.gt.maxspec) goto 994
  if (old) call skplin(2,unitreleases)
  do i=1,nspec
    read(unitreleases,*,err=998) specnum_rel
    if (old) call skplin(2,unitreleases)
    call readspecies(specnum_rel,i)

  ! For backward runs, only 1 species is allowed
  !*********************************************

  !if ((ldirect.lt.0).and.(nspec.gt.1)) then
  !write(*,*) '#####################################################'
  !write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  !write(*,*) '#### FOR BACKWARD RUNS, ONLY 1 SPECIES IS ALLOWED####'
  !write(*,*) '#####################################################'
  !  stop
  !endif

  ! Molecular weight
  !*****************

    if (((iout.eq.2).or.(iout.eq.3)).and. &
         (weightmolar(i).lt.0.)) then
      write(*,*) 'For mixing ratio output, valid molar weight'
      write(*,*) 'must be specified for all simulated species.'
      write(*,*) 'Check table SPECIES or choose concentration'
      write(*,*) 'output instead if molar weight is not known.'
      stop
    endif


  ! Radioactive decay
  !******************

    decay(i)=0.693147/decay(i) !conversion half life to decay constant


  ! Dry deposition of gases
  !************************

    if (reldiff(i).gt.0.) &
         rm(i)=1./(henry(i)/3000.+100.*f0(i))    ! mesophyll resistance

  ! Dry deposition of particles
  !****************************

    vsetaver(i)=0.
    cunningham(i)=0.
    dquer(i)=dquer(i)*1000000.         ! Conversion m to um
    if (density(i).gt.0.) then                  ! Additional parameters
     call part0(dquer(i),dsigma(i),density(i),fracth,schmih,cun,vsh)
      do j=1,ni
        fract(i,j)=fracth(j)
        schmi(i,j)=schmih(j)
        vset(i,j)=vsh(j)
        cunningham(i)=cunningham(i)+cun*fract(i,j)
        vsetaver(i)=vsetaver(i)-vset(i,j)*fract(i,j)
      end do
      write(*,*) 'Average setting velocity: ',i,vsetaver(i)
    endif

  ! Dry deposition for constant deposition velocity
  !************************************************

    dryvel(i)=dryvel(i)*0.01         ! conversion to m/s

  ! Check if wet deposition or OH reaction shall be calculated
  !***********************************************************
    if (weta(i).gt.0.)  then
      WETDEP=.true.
      write (*,*) 'Wetdeposition switched on: ',weta(i),i
    endif
    if (ohreact(i).gt.0) then
      OHREA=.true.
      write (*,*) 'OHreaction switched on: ',ohreact(i),i
    endif


    if ((reldiff(i).gt.0.).or.(density(i).gt.0.).or. &
         (dryvel(i).gt.0.)) then
      DRYDEP=.true.
      DRYDEPSPEC(i)=.true.
    endif

  end do

    if (WETDEP.or.DRYDEP) DEP=.true.

  ! Read specifications for each release point
  !*******************************************

  numpoint=0
  numpartmax=0
  releaserate=0.
1000   numpoint=numpoint+1
  read(unitreleases,*,end=250)
  read(unitreleases,*,err=998,end=250) id1,it1
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) id2,it2
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xpoint1(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) ypoint1(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) xpoint2(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) ypoint2(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) kindz(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) zpoint1(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) zpoint2(numpoint)
  if (old) call skplin(2,unitreleases)
  read(unitreleases,*,err=998) npart(numpoint)
  if (old) call skplin(2,unitreleases)
  do i=1,nspec
    read(unitreleases,*,err=998) xmass(numpoint,i)
    if (old) call skplin(2,unitreleases)
  end do
  !save compoint only for the first 1000 release points
  if (numpoint.le.1000) then
    read(unitreleases,'(a40)',err=998) compoint(numpoint)(1:40)
  else
    read(unitreleases,'(a40)',err=998) compoint(1001)(1:40)
  endif
  if (old) call skplin(1,unitreleases)
  if (numpoint.le.1000) then
    if((xpoint1(numpoint).eq.0.).and.(ypoint1(numpoint).eq.0.).and. &
         (xpoint2(numpoint).eq.0.).and.(ypoint2(numpoint).eq.0.).and. &
         (compoint(numpoint)(1:8).eq.'        ')) goto 250
  else
    if((xpoint1(numpoint).eq.0.).and.(ypoint1(numpoint).eq.0.).and. &
         (xpoint2(numpoint).eq.0.).and.(ypoint2(numpoint).eq.0.)) goto 250
  endif


  ! If a release point contains no particles, stop and issue error message
  !***********************************************************************

  if (npart(numpoint).eq.0) then
    write(*,*) 'FLEXPART MODEL ERROR'
    write(*,*) 'RELEASES file is corrupt.'
    write(*,*) 'At least for one release point, there are zero'
    write(*,*) 'particles released. Make changes to RELEASES.'
    stop
  endif

  ! Check whether x coordinates of release point are within model domain
  !*********************************************************************

   if (xpoint1(numpoint).lt.xlon0) &
        xpoint1(numpoint)=xpoint1(numpoint)+360.
   if (xpoint1(numpoint).gt.xlon0+(nxmin1)*dx) &
        xpoint1(numpoint)=xpoint1(numpoint)-360.
   if (xpoint2(numpoint).lt.xlon0) &
        xpoint2(numpoint)=xpoint2(numpoint)+360.
   if (xpoint2(numpoint).gt.xlon0+(nxmin1)*dx) &
        xpoint2(numpoint)=xpoint2(numpoint)-360.

  ! Determine relative beginning and ending times of particle release
  !******************************************************************

  jul1=juldate(id1,it1)
  jul2=juldate(id2,it2)
  if (jul1.gt.jul2) then
    write(*,*) 'FLEXPART MODEL ERROR'
    write(*,*) 'Release stops before it begins.'
    write(*,*) 'Make changes to file RELEASES.'
    stop
  endif
  if (mdomainfill.eq.0) then   ! no domain filling
    if (ldirect.eq.1) then
      if ((jul1.lt.bdate).or.(jul2.gt.edate)) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'Release starts before simulation begins or ends'
        write(*,*) 'after simulation stops.'
        write(*,*) 'Make files COMMAND and RELEASES consistent.'
        stop
      endif
      ireleasestart(numpoint)=int((jul1-bdate)*86400.)
      ireleaseend(numpoint)=int((jul2-bdate)*86400.)
    else if (ldirect.eq.-1) then
      if ((jul1.lt.edate).or.(jul2.gt.bdate)) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'Release starts before simulation begins or ends'
        write(*,*) 'after simulation stops.'
        write(*,*) 'Make files COMMAND and RELEASES consistent.'
        stop
      endif
      ireleasestart(numpoint)=int((jul1-bdate)*86400.)
      ireleaseend(numpoint)=int((jul2-bdate)*86400.)
    endif
  endif


  ! Check, whether the total number of particles may exceed totally allowed
  ! number of particles at some time during the simulation
  !************************************************************************

  ! Determine the release rate (particles per second) and total number
  ! of particles released during the simulation
  !*******************************************************************

  if (ireleasestart(numpoint).ne.ireleaseend(numpoint)) then
    releaserate=releaserate+real(npart(numpoint))/ &
         real(ireleaseend(numpoint)-ireleasestart(numpoint))
  else
    releaserate=99999999
  endif
  numpartmax=numpartmax+npart(numpoint)
  goto 1000


250   close(unitreleases)

  write (*,*) ' Particles allocated for this run: ',maxpart, ', released in simulation: ',  numpartmax
  numpoint=numpoint-1

  if (ioutputforeachrelease.eq.1) then
    maxpointspec_act=numpoint
  else
    maxpointspec_act=1
  endif

  if (releaserate.gt. &
       0.99*real(maxpart)/real(lage(nageclass))) then
    if (numpartmax.gt.maxpart) then
  write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  write(*,*) '####                                             ####'
  write(*,*) '####WARNING - TOTAL NUMBER OF PARTICLES SPECIFIED####'
  write(*,*) '#### IN FILE "RELEASES" MAY AT SOME POINT DURING ####'
  write(*,*) '#### THE SIMULATION EXCEED THE MAXIMUM ALLOWED   ####'
  write(*,*) '#### NUMBER (MAXPART).IF RELEASES DO NOT OVERLAP,####'
  write(*,*) '#### FLEXPART CAN POSSIBLY COMPLETE SUCCESSFULLY.####'
  write(*,*) '#### HOWEVER, FLEXPART MAY HAVE TO STOP          ####'
  write(*,*) '#### AT SOME TIME DURING THE SIMULATION. PLEASE  ####'
  write(*,*) '#### MAKE SURE THAT YOUR SETTINGS ARE CORRECT.   ####'
  write(*,*) '#####################################################'
      write(*,*) 'Maximum release rate may be: ',releaserate, &
           ' particles per second'
      write(*,*) 'Maximum allowed release rate is: ', &
           real(maxpart)/real(lage(nageclass)),' particles per second'
      write(*,*) &
           'Total number of particles released during the simulation is: ', &
           numpartmax
      write(*,*) 'Maximum allowed number of particles is: ',maxpart
    endif
  endif

  return

994   write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR - MAXIMUM NUMBER OF EMITTED SPECIES IS####'
  write(*,*) '#### TOO LARGE. PLEASE REDUCE NUMBER OF SPECIES. ####'
  write(*,*) '#####################################################'
  stop

998   write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  write(*,*) '####                                             ####'
  write(*,*) '#### FATAL ERROR - FILE "RELEASES" IS            ####'
  write(*,*) '#### CORRUPT. PLEASE CHECK YOUR INPUTS FOR       ####'
  write(*,*) '#### MISTAKES OR GET A NEW "RELEASES"-           ####'
  write(*,*) '#### FILE ...                                    ####'
  write(*,*) '#####################################################'
  stop


999   write(*,*) '#####################################################'
  write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES: '
  write(*,*)
  write(*,*) 'FATAL ERROR - FILE CONTAINING PARTICLE RELEASE POINTS'
  write(*,*) 'POINTS IS NOT AVAILABLE OR YOU ARE NOT'
  write(*,*) 'PERMITTED FOR ANY ACCESS'
  write(*,*) '#####################################################'
  stop

end subroutine readreleases
