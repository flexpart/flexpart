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

subroutine readspecies(id_spec,pos_spec)

  !*****************************************************************************
  !                                                                            *
  !     This routine reads names and physical constants of chemical species/   *
  !     radionuclides given in the parameter pos_spec                          *
  !                                                                            *
  !   Author: A. Stohl                                                         *
  !                                                                            *
  !   11 July 1996                                                             *
  !                                                                            *
  !   Changes:                                                                 *
  !   N. Kristiansen, 31.01.2013: Including parameters for in-cloud scavenging *
  !                                                                            *
  !   HSO, 13 August 2013
  !   added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decaytime(maxtable)   half time for radiological decay                     *
  ! specname(maxtable)    names of chemical species, radionuclides             *
  ! weta_gas, wetb_gas    Parameters for below-cloud scavenging of gasses      *
  ! crain_aero,csnow_aero Parameters for below-cloud scavenging of aerosols    *
  ! ccn_aero,in_aero      Parameters for in-cloud scavenging of aerosols       *
  ! ohcconst              OH reaction rate constant C                          *
  ! ohdconst              OH reaction rate constant D                          *
  ! ohnconst              OH reaction rate constant n                          *
  ! id_spec               SPECIES number as referenced in RELEASE file         *
  ! id_pos                position where SPECIES data shall be stored          *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: i, pos_spec,j
  integer :: idow,ihour,id_spec
  character(len=3) :: aspecnumb
  logical :: spec_found

  character(len=16) :: pspecies
  real :: pdecay, pweta_gas, pwetb_gas, preldiff, phenry, pf0, pdensity, pdquer
  real :: pdsigma, pdryvel, pweightmolar, pohcconst, pohdconst, pohnconst
  real :: pcrain_aero, pcsnow_aero, pccn_aero, pin_aero
  real :: parea_dow(7), parea_hour(24), ppoint_dow(7), ppoint_hour(24)
  integer :: readerror

! declare namelist
  namelist /species_params/ &
       pspecies, pdecay, pweta_gas, pwetb_gas, &
       pcrain_aero, pcsnow_aero, pccn_aero, pin_aero, &
       preldiff, phenry, pf0, pdensity, pdquer, &
       pdsigma, pdryvel, pweightmolar, pohcconst, pohdconst, pohnconst, &
       parea_dow, parea_hour, ppoint_dow, ppoint_hour

  pspecies="" ! read failure indicator value
  pdecay=-999.9
  pweta_gas=-9.9E-09
  pwetb_gas=0.0
  pcrain_aero=-9.9E-09
  pcsnow_aero=-9.9E-09
  pccn_aero=-9.9E-09
  pin_aero=-9.9E-09
  preldiff=-9.9
  phenry=0.0
  pf0=0.0
  pdensity=-9.9E09
  pdquer=0.0
  pdsigma=0.0
  pdryvel=-9.99
  pohcconst=-9.99
  pohdconst=-9.9E-09
  pohnconst=2.0
  pweightmolar=-999.9
  parea_dow=-999.9
  parea_hour=-999.9
  ppoint_dow=-999.9
  ppoint_hour=-999.9


  do j=1,24           ! initialize everything to no variation
    parea_hour(j)=1.
    ppoint_hour(j)=1.
    area_hour(pos_spec,j)=1.
    point_hour(pos_spec,j)=1.
  end do
  do j=1,7
    parea_dow(j)=1.
    ppoint_dow(j)=1.
    area_dow(pos_spec,j)=1.
    point_dow(pos_spec,j)=1.
  end do

! Open the SPECIES file and read species names and properties
!************************************************************
  specnum(pos_spec)=id_spec
  write(aspecnumb,'(i3.3)') specnum(pos_spec)
  open(unitspecies,file=path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb,status='old',form='formatted',err=998)
!write(*,*) 'reading SPECIES',specnum(pos_spec)

  ASSSPEC=.FALSE.

! try namelist input
  read(unitspecies,species_params,iostat=readerror)
  close(unitspecies)

  if ((len(trim(pspecies)).eq.0).or.(readerror.ne.0)) then ! no namelist found
    if (lroot) write(*,*) "SPECIES file not in NAMELIST format, attempting to &
         &read as fixed format"

    readerror=1

    open(unitspecies,file=path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb,status='old',err=998)

    do i=1,6
      read(unitspecies,*)
    end do

    read(unitspecies,'(a10)',end=22) species(pos_spec)
!  write(*,*) species(pos_spec)
    read(unitspecies,'(f18.1)',end=22) decay(pos_spec)
!  write(*,*) decay(pos_spec)
    read(unitspecies,'(e18.1)',end=22) weta_gas(pos_spec)
!  write(*,*) weta_gas(pos_spec)
    read(unitspecies,'(f18.2)',end=22) wetb_gas(pos_spec)
!  write(*,*) wetb_gas(pos_spec)
    read(unitspecies,'(e18.1)',end=22) crain_aero(pos_spec)
!  write(*,*) crain_aero(pos_spec)
    read(unitspecies,'(f18.2)',end=22) csnow_aero(pos_spec)
!  write(*,*) csnow_aero(pos_spec)
!*** NIK 31.01.2013: including in-cloud scavening parameters
    read(unitspecies,'(e18.1)',end=22) ccn_aero(pos_spec)
!  write(*,*) ccn_aero(pos_spec)
    read(unitspecies,'(f18.2)',end=22) in_aero(pos_spec)
!  write(*,*) in_aero(pos_spec)
    read(unitspecies,'(f18.1)',end=22) reldiff(pos_spec)
!  write(*,*) reldiff(pos_spec)
    read(unitspecies,'(e18.1)',end=22) henry(pos_spec)
!  write(*,*) henry(pos_spec)
    read(unitspecies,'(f18.1)',end=22) f0(pos_spec)
!  write(*,*) f0(pos_spec)
    read(unitspecies,'(e18.1)',end=22) density(pos_spec)
!  write(*,*) density(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dquer(pos_spec)
!  write(*,*) 'dquer(pos_spec):', dquer(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dsigma(pos_spec)
!  write(*,*) dsigma(pos_spec)
    read(unitspecies,'(f18.2)',end=22) dryvel(pos_spec)
!  write(*,*) dryvel(pos_spec)
    read(unitspecies,'(f18.2)',end=22) weightmolar(pos_spec)
!  write(*,*) weightmolar(pos_spec)
    read(unitspecies,'(e18.2)',end=22) ohcconst(pos_spec)
!  write(*,*) ohcconst(pos_spec)
    read(unitspecies,'(f8.2)',end=22) ohdconst(pos_spec)
!  write(*,*) ohdconst(pos_spec)
    read(unitspecies,'(f8.2)',end=22) ohnconst(pos_spec)
!  write(*,*) ohnconst(pos_spec)

! Read in daily and day-of-week variation of emissions, if available
!*******************************************************************

    read(unitspecies,*,end=22)
    do j=1,24     ! 24 hours, starting with 0-1 local time
      read(unitspecies,*) ihour,area_hour(pos_spec,j),point_hour(pos_spec,j)
    end do
    read(unitspecies,*)
    do j=1,7      ! 7 days of the week, starting with Monday
      read(unitspecies,*) idow,area_dow(pos_spec,j),point_dow(pos_spec,j)
    end do

    pspecies=species(pos_spec)
    pdecay=decay(pos_spec)
    pweta_gas=weta_gas(pos_spec)
    pwetb_gas=wetb_gas(pos_spec)
    pcrain_aero=crain_aero(pos_spec)
    pcsnow_aero=csnow_aero(pos_spec)
    pccn_aero=ccn_aero(pos_spec)
    pin_aero=in_aero(pos_spec)
    preldiff=reldiff(pos_spec)
    phenry=henry(pos_spec)
    pf0=f0(pos_spec)
    pdensity=density(pos_spec)
    pdquer=dquer(pos_spec)
    pdsigma=dsigma(pos_spec)
    pdryvel=dryvel(pos_spec)
    pweightmolar=weightmolar(pos_spec)
    pohcconst=ohcconst(pos_spec)
    pohdconst=ohdconst(pos_spec)
    pohnconst=ohnconst(pos_spec)


    do j=1,24     ! 24 hours, starting with 0-1 local time
      parea_hour(j)=area_hour(pos_spec,j)
      ppoint_hour(j)=point_hour(pos_spec,j)
    end do
    do j=1,7      ! 7 days of the week, starting with Monday
      parea_dow(j)=area_dow(pos_spec,j)
      ppoint_dow(j)=point_dow(pos_spec,j)
    end do

  else ! namelist available

    species(pos_spec)=pspecies
    decay(pos_spec)=pdecay
    weta_gas(pos_spec)=pweta_gas
    wetb_gas(pos_spec)=pwetb_gas
    crain_aero(pos_spec)=pcrain_aero
    csnow_aero(pos_spec)=pcsnow_aero
    ccn_aero(pos_spec)=pccn_aero
    in_aero(pos_spec)=pin_aero
    reldiff(pos_spec)=preldiff
    henry(pos_spec)=phenry
    f0(pos_spec)=pf0
    density(pos_spec)=pdensity
    dquer(pos_spec)=pdquer
    dsigma(pos_spec)=pdsigma
    dryvel(pos_spec)=pdryvel
    weightmolar(pos_spec)=pweightmolar
    ohcconst(pos_spec)=pohcconst
    ohdconst(pos_spec)=pohdconst
    ohnconst(pos_spec)=pohnconst

    do j=1,24     ! 24 hours, starting with 0-1 local time
      area_hour(pos_spec,j)=parea_hour(j)
      point_hour(pos_spec,j)=ppoint_hour(j)
    end do
    do j=1,7      ! 7 days of the week, starting with Monday
      area_dow(pos_spec,j)=parea_dow(j)
      point_dow(pos_spec,j)=ppoint_dow(j)
    end do

  endif

  i=pos_spec

!NIK 16.02.2015
! Check scavenging parameters given in SPECIES file

  if (lroot) then
! ZHG 2016.04.07 Start of changes
    write(*,*) ' '
    if (dquer(pos_spec) .gt.0)  write(*,'(a,i3,a,a,a)')       ' SPECIES: ', &
         &id_spec,'  ', species(pos_spec),'  (AEROSOL) '
    if (dquer(pos_spec) .le.0)  write(*,'(a,i3,a,a,a)')       ' SPECIES: ', &
         &id_spec,'  ', species(pos_spec),'  (GAS) '

! Particles
!**********
    if (dquer(pos_spec).gt.0) then
      if (ccn_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle CCN  efficiency (CCNeff):', ccn_aero(pos_spec)
      else 
        write(*,'(a)')      '  Particle CCN  efficiency (CCNeff):   OFF'
      endif
      if (in_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle  IN  efficiency (INeff) :', in_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle  IN  efficiency (INeff) :   OFF'
      endif
      if (crain_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle Rain efficiency (Crain) :', crain_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle Rain efficiency (Crain) :   OFF'
      endif
      if (csnow_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle Snow efficiency (Csnow) :', csnow_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle Snow efficiency (Csnow) :   OFF'
      end if
      if (density(pos_spec) .gt. 0) then
        write(*,'(a)') '  Dry deposition is turned         :   ON'
        if (reldiff(pos_spec).gt.0) then
           stop 'density>0 (SPECIES is a particle) implies reldiff <=0  '
        endif
      else
        if (reldiff(pos_spec).le.0) then
           stop 'density<=0 (SPECIES is a gas) implies reldiff >0  '
        endif      
        write(*,'(a)') '  Dry deposition is (density<0)    :   OFF'
      end if
      if (crain_aero(pos_spec).gt.10.0 .or. csnow_aero(pos_spec).gt.10.0 .or. &
           &ccn_aero(pos_spec).gt.1.0 .or. in_aero(pos_spec).gt.1.0) then
        write(*,*) '*******************************************'
        write(*,*) ' WARNING: Particle Scavenging parameter likely out of range '
        write(*,*) '       Likely   range for Crain    0.0-10'
        write(*,*) '       Likely   range for Csnow    0.0-10'
        write(*,*) '       Physical range for CCNeff   0.0-1'
        write(*,*) '       Physical range for INeff    0.0-1'
        write(*,*) '*******************************************'
      end if
    else
! Gas
!****
      if (weta_gas(pos_spec) .gt. 0 .and. wetb_gas(pos_spec).gt.0) then
        write(*,*)          '  Wet removal for gases      is turned: ON'
        write(*,*)          '  Gas below-cloud scavenging parameter A  ', &
             &weta_gas(pos_spec)
        write(*,'(a,f5.2)') '  Gas below-cloud scavenging parameter B  ', &
             &wetb_gas(pos_spec)
      else
        write(*,*)          '  Wet removal for gases      is turned: OFF '
      end if
      if (reldiff(i).gt.0.) then
        write(*,*)          '  Dry deposition for gases   is turned: ON '
      else
        write(*,*)          '  Dry deposition for gases   is turned: OFF '
      end if
      if (weta_gas(pos_spec).gt.0.) then !if wet deposition is turned on
        if (weta_gas(pos_spec).gt.1E-04 .or. weta_gas(pos_spec).lt.1E-09 .or. &
             &wetb_gas(pos_spec).gt.0.8 .or. wetb_gas(pos_spec).lt.0.4) then
          write(*,*) '*******************************************'
          write(*,*) ' WARNING: Gas below-cloud scavengig is out of likely range'
          write(*,*) '          Likely range for A is 1E-04 to 1E-08'
          write(*,*) '          Likely range for B is 0.60  to 0.80 ' 
          write(*,*) '*******************************************'
        end if
      endif

      if (((weta_gas(pos_spec).gt.0).or.(wetb_gas(pos_spec).gt.0)).and.&
           &(henry(pos_spec).le.0)) then
        if (dquer(pos_spec).le.0) goto 996 ! no particle, no henry set
      endif
    end if
  end if

  !  if (dsigma(i).eq.0.) dsigma(i)=1.0001   ! avoid floating exception
  if (dquer(i).gt.0 .and. dsigma(i).le.1.) then !dsigma(i)=1.0001   ! avoid floating exception
    !write(*,*) '#### FLEXPART MODEL ERROR!                      ####'
    write(*,*) '#### FLEXPART MODEL WARNING                     ####'
    write(*,*) '#### in SPECIES_',aspecnumb, '                             ####'
    write(*,*) '#### from v10.4 dsigma has to be larger than 1  ####'  
    write(*,*) '#### to adapt older SPECIES files,              ####' 
    write(*,*) '#### if dsigma was < 1                          ####' 
    write(*,*) '#### use the reciprocal of the old dsigma       ####'  
    if (.not.debug_mode) then 
       stop
    else
       write(*,*) 'debug mode: continue'
    endif
  endif

  if ((reldiff(i).gt.0.).and.(density(i).gt.0.)) then
    write(*,*) '#### FLEXPART MODEL ERROR! FILE "SPECIES"    ####'
    write(*,*) '#### IS CORRUPT. SPECIES CANNOT BE BOTH      ####'
    write(*,*) '#### PARTICLE AND GAS.                       ####'
    write(*,*) '#### SPECIES NUMBER',aspecnumb
    stop
  endif
20 continue


22 close(unitspecies)

! namelist output if requested
  if (nmlout.and.lroot) then
    open(unitspecies,file=path(2)(1:length(2))//'SPECIES_'//aspecnumb//'.namelist',access='append',status='replace',err=1000)
    write(unitspecies,nml=species_params)
    close(unitspecies)
  endif

  return

996 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### WET DEPOSITION SWITCHED ON, BUT NO HENRYS  #### '
  write(*,*) '#### CONSTANT IS SET                            ####'
  write(*,*) '#### PLEASE MODIFY SPECIES DESCR. FILE!        #### '
  write(*,*) '#####################################################'
  stop


997 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### THE ASSSOCIATED SPECIES HAS TO BE DEFINED  #### '
  write(*,*) '#### BEFORE THE ONE WHICH POINTS AT IT          #### '
  write(*,*) '#### PLEASE CHANGE ORDER IN RELEASES OR ADD     #### '
  write(*,*) '#### THE ASSOCIATED SPECIES IN RELEASES         #### '
  write(*,*) '#####################################################'
  stop


998 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### THE SPECIES FILE FOR SPECIES ', id_spec
  write(*,*) '#### CANNOT BE FOUND: CREATE FILE'
  write(*,*) '#### ',path(1)(1:length(1)),'SPECIES/SPECIES_',aspecnumb
  write(*,*) '#####################################################'
  stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "SPECIES_',aspecnumb,'.namelist'
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop

end subroutine readspecies
