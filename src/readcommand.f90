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

subroutine readcommand

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the current model run.  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !                                                                            *
  !     modified by M. Cassiani, 2016                                          *
  !     added flag for alterntive methods to compute vertical velocity          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! bdate                beginning date as Julian date                         *
  ! ctl                  factor by which time step must be smaller than        *
  !                      Lagrangian time scale                                 *
  ! ibdate,ibtime        beginnning date and time (YYYYMMDD, HHMISS)           *
  ! ideltas [s]          modelling period                                      *
  ! iedate,ietime        ending date and time (YYYYMMDD, HHMISS)               *
  ! ifine                reduction factor for vertical wind time step          *
  ! outputforeachrel     for forward runs it is possible either to create      *
  !                      one outputfield or several for each releasepoint      *
  ! iflux                switch to turn on (1)/off (0) flux calculations       *
  ! iout                 1 for conc. (residence time for backward runs) output,*
  !                      2 for mixing ratio output, 3 both, 4 for plume        *
  !                      trajectory output, 5 = options 1 and 4                *
  ! ipin                 1 continue simulation with dumped particle data, 0 no *
  ! ipout                0 no particle dump, 1 every output time, 3 only at end*
  ! itsplit [s]          time constant for particle splitting                  *
  ! loutaver [s]         concentration output is an average over loutaver      *
  !                      seconds                                               *
  ! loutsample [s]       average is computed from samples taken every [s]      *
  !                      seconds                                               *
  ! loutstep [s]         time interval of concentration output                 *
  ! lsynctime [s]        synchronisation time interval for all particles       *
  ! lagespectra          switch to turn on (1)/off (0) calculation of age      *
  !                      spectra                                               *
  ! lconvection          value of either 0 and 1 indicating mixing by          *
  !                      convection                                            *
  !                      = 0 .. no convection                                  *
  !                      + 1 .. parameterisation of mixing by subgrid-scale    *
  !                              convection = on                               *
  ! lsubgrid             switch to turn on (1)/off (0) subgrid topography      *
  !                      parameterization                                      *
  ! method               method used to compute the particle pseudovelocities  *
  ! mdomainfill          1 use domain-filling option, 0 not, 2 use strat. O3   *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitcommand          unit connected to file COMMAND                        *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  real(kind=dp) :: juldate
  character(len=50) :: line
  logical :: old


  ! Open the command file and read user options
  !********************************************


  open(unitcommand,file=path(1)(1:length(1))//'COMMAND',status='old', &
       err=999)

  ! Check the format of the COMMAND file (either in free format,
  ! or using formatted mask)
  ! Use of formatted mask is assumed if line 10 contains the word 'DIRECTION'
  !**************************************************************************

  call skplin(9,unitcommand)
  read (unitcommand,901) line
901   format (a)
  if (index(line,'LDIRECT') .eq. 0) then
    old = .false.
  else
    old = .true.
  endif
  rewind(unitcommand)

  ! Read parameters
  !****************

  call skplin(7,unitcommand)
  if (old) call skplin(1,unitcommand)

  read(unitcommand,*) ldirect
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ibdate,ibtime
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) iedate,ietime
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) loutstep
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) loutaver
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) loutsample
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) itsplit
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) lsynctime
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ctl
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ifine
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) iout
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ipout
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) lsubgrid
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) lconvection
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) lagespectra
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ipin
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ioutputforeachrelease
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) iflux
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) mdomainfill
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ind_source
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) ind_receptor
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) mquasilag
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) nested_output
  if (old) call skplin(3,unitcommand)
  read(unitcommand,*) linit_cond
  if (old) call skplin(3,unitcommand) !added by mc
  read(unitcommand,*) method_w
  close(unitcommand)


  close(unitcommand)

  ifine=max(ifine,1)

 
  ! Determine how Markov chain is formulated (for w or for w/sigw)
  !***************************************************************

  if (ctl.ge.0.1) then
    turbswitch=.true.
  else
    turbswitch=.false.
    ifine=1
  endif
  fine=1./real(ifine)
  ctl=1./ctl

  ! Set the switches required for the various options for input/output units
  !*************************************************************************
  !AF Set the switches IND_REL and IND_SAMP for the release and sampling
  !Af switches for the releasefile:
  !Af IND_REL =  1 : xmass * rho
  !Af IND_REL =  0 : xmass * 1

  !Af switches for the conccalcfile:
  !AF IND_SAMP =  0 : xmass * 1
  !Af IND_SAMP = -1 : xmass / rho

  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af   NOTE that in backward simulations the release of computational particles
  !Af   takes place at the "receptor" and the sampling of p[articles at the "source".
  !Af          1 = mass units
  !Af          2 = mass mixing ratio units
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1 = mass units
  !Af          2 = mass mixing ratio units

  if ( ldirect .eq. 1 ) then  ! FWD-Run
  !Af set release-switch
     if (ind_source .eq. 1 ) then !mass
        ind_rel = 0
     else ! mass mix
        ind_rel = 1
     endif
  !Af set sampling switch
     if (ind_receptor .eq. 1) then !mass
        ind_samp = 0
     else ! mass mix
        ind_samp = -1
     endif
  elseif (ldirect .eq. -1 ) then !BWD-Run
  !Af set sampling switch
     if (ind_source .eq. 1 ) then !mass
        ind_samp = -1
     else ! mass mix
        ind_samp = 0
     endif
  !Af set release-switch
     if (ind_receptor .eq. 1) then !mass
        ind_rel = 1
     else ! mass mix
        ind_rel = 0
     endif
  endif

  !*************************************************************
  ! Check whether valid options have been chosen in file COMMAND
  !*************************************************************

  ! Check options for initial condition output: Switch off for forward runs
  !************************************************************************

  if (ldirect.eq.1) linit_cond=0
  if ((linit_cond.lt.0).or.(linit_cond.gt.2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INVALID OPTION    #### '
    write(*,*) ' #### FOR LINIT_COND IN FILE "COMMAND".       #### '
    stop
  endif

  ! Check input dates
  !******************

  if (iedate.lt.ibdate) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING DATE    #### '
    write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
    write(*,*) ' #### "COMMAND".                              #### '
    stop
  else if (iedate.eq.ibdate) then
    if (ietime.lt.ibtime) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### '
    write(*,*) ' #### IS LARGER THAN ENDING TIME. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
    write(*,*) ' #### "COMMAND".                              #### '
    stop
    endif
  else if (bdate.lt.yearorigin) then !modified by mc
   write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### ' !modified by mc
   write(*,*) ' #### must be larger than yearzero: Difference in time    ####'  !modified by mc
   write(*,*) ' #### computed using no-leap calendar from an ####'  !modified by mc
   write(*,*) ' #### arbitrary starting yearorigin           #### ' !modified by mc
   write(*,*) ' #### "COMMAND".                              #### ' !modified by mc
    stop
   endif

  ! Determine kind of dispersion method
  !************************************

  if (ctl.gt.0.) then
    method=1
    mintime=minstep
  else
    method=0
    mintime=lsynctime
  endif

  ! Check whether a valid option for gridded model output has been chosen
  !**********************************************************************

  if ((iout.lt.1).or.(iout.gt.5)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### IOUT MUST BE 1, 2, 3, 4, OR 5!          #### '
    stop
  endif

  !AF check consistency between units and volume mixing ratio
  if ( ((iout.eq.2).or.(iout.eq.3)).and. &
       (ind_source.gt.1 .or.ind_receptor.gt.1) ) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### VOLUME MIXING RATIO ONLY SUPPORTED      #### '
    write(*,*) ' #### FOR MASS UNITS (at the moment)          #### '
    stop
  endif



  ! For quasilag output for each release is forbidden
  !*****************************************************************************

  if ((ioutputforeachrelease.eq.1).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### OUTPUTFOREACHRELEASE AND QUASILAGRANGIAN####'
      write(*,*) '#### MODE IS NOT POSSIBLE   !                ####'
      stop
  endif


  ! For quasilag backward is forbidden
  !*****************************************************************************

  if ((ldirect.lt.0).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, QUASILAGRANGIAN MODE ####'
      write(*,*) '#### IS NOT POSSIBLE   !                     ####'
      stop
  endif


  ! For backward runs one releasefield for all releases makes no sense,
  ! For quasilag and domainfill ioutputforechrelease is forbidden
  !*****************************************************************************

  if ((ldirect.lt.0).and.(ioutputforeachrelease.eq.0)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUTPUTFOREACHRLEASE ####'
      write(*,*) '#### MUST BE SET TO ONE!                     ####'
      stop
  endif


  ! For backward runs one releasefield for all releases makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if ((mdomainfill.eq.1).and.(ioutputforeachrelease.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR DOMAIN FILLING RUNS OUTPUT FOR      ####'
      write(*,*) '#### EACH RELEASE IS FORBIDDEN !             ####'
      stop
  endif


  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! For backward runs, only residence time output (iout=1) or plume trajectories (iout=4),
  ! or both (iout=5) makes sense; other output options are "forbidden"
  !*****************************************************************************

  if (ldirect.lt.0) then
    if ((iout.eq.2).or.(iout.eq.3)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUT MUST BE 1,4,OR 5####'
      stop
    endif
  endif


  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if (mdomainfill.ge.1) then
    if ((iout.eq.4).or.(iout.eq.5)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR DOMAIN-FILLING TRAJECTORY OPTION,   ####'
      write(*,*) '#### IOUT MUST NOT BE SET TO 4 OR 5.         ####'
      stop
    endif
  endif



  ! Check whether a valid options for particle dump has been chosen
  !****************************************************************

  if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### IPOUT MUST BE 1, 2 OR 3!                #### '
    stop
  endif

  if(lsubgrid.ne.1) then
    write(*,*) '             ----------------               '
    write(*,*) ' INFORMATION: SUBGRIDSCALE TERRAIN EFFECT IS'
    write(*,*) ' NOT PARAMETERIZED DURING THIS SIMULATION.  '
    write(*,*) '             ----------------               '
  endif


  ! Check whether convection scheme is either turned on or off
  !***********************************************************

  if ((lconvection.ne.0).and.(lconvection.ne.1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### LCONVECTION MUST BE SET TO EITHER 1 OR 0#### '
    stop
  endif


  ! Check whether synchronisation interval is sufficiently short
  !*************************************************************

  if (lsynctime.gt.(idiffnorm/2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SYNCHRONISATION   #### '
    write(*,*) ' #### TIME IS TOO LONG. MAKE IT SHORTER.      #### '
    write(*,*) ' #### MINIMUM HAS TO BE: ', idiffnorm/2
    stop
  endif


  ! Check consistency of the intervals, sampling periods, etc., for model output
  !*****************************************************************************

  if (loutaver.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### ZERO.                                   #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    stop
  endif

  if (loutaver.gt.loutstep) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN INTERVAL OF OUTPUT.        #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    stop
  endif

  if (loutsample.gt.loutaver) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN TIME AVERAGE OF OUTPUT.    #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    stop
  endif

  if (mod(loutaver,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if ((loutaver/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE AT LEAST    #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    stop
  endif

  if (mod(loutstep,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE A MULTIPLE #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if ((loutstep/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE AT LEAST   #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    stop
  endif

  if (mod(loutsample,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    stop
  endif

  if (itsplit.lt.loutaver) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SPLITTING TIME FOR#### '
    write(*,*) ' #### PARTICLES IS TOO SHORT. PLEASE INCREASE #### '
    write(*,*) ' #### SPLITTING TIME CONSTANT.                #### '
    stop
  endif

  if ((mquasilag.eq.1).and.(iout.ge.4)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! CONFLICTING       #### '
    write(*,*) ' #### OPTIONS: IF MQUASILAG=1, PLUME          #### '
    write(*,*) ' #### TRAJECTORY OUTPUT IS IMPOSSIBLE.        #### '
    stop
  endif

  ! Compute modeling time in seconds and beginning date in Julian date
  !*******************************************************************

  outstep=real(abs(loutstep))
  if (ldirect.eq.1) then
    bdate=juldate(ibdate,ibtime)
    edate=juldate(iedate,ietime)
    ideltas=nint((edate-bdate)*86400.)
  else if (ldirect.eq.-1) then
    loutaver=-1*loutaver
    loutstep=-1*loutstep
    loutsample=-1*loutsample
    lsynctime=-1*lsynctime
    bdate=juldate(iedate,ietime)
    edate=juldate(ibdate,ibtime)
    ideltas=nint((edate-bdate)*86400.)
  else
    write(*,*) ' #### FLEXPART MODEL ERROR! DIRECTION IN      #### '
    write(*,*) ' #### FILE "COMMAND" MUST BE EITHER -1 OR 1.  #### '
    stop
  endif

  return

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "COMMAND"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

end subroutine readcommand
