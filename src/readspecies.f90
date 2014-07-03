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
  !
  !   Changes:                                                                 *
  !   N. Kristiansen, 31.01.2013: Including parameters for in-cloud scavenging *
  !                                                                            *
  !   HSO, 13 August 2013
  !   added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decaytime(maxtable)  half time for radiological decay                      *
  ! specname(maxtable)   names of chemical species, radionuclides              *
  ! weta, wetb           Parameters for determining below-cloud scavenging     *
  ! weta_in              Parameter for determining in-cloud scavenging         *
  ! wetb_in              Parameter for determining in-cloud scavenging         *
  ! wetc_in              Parameter for determining in-cloud scavenging         *
  ! wetd_in              Parameter for determining in-cloud scavenging         *
  ! ohreact              OH reaction rate                                      *
  ! id_spec              SPECIES number as referenced in RELEASE file          *
  ! id_pos               position where SPECIES data shall be stored           *
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
  real :: pdecay, pweta, pwetb, preldiff, phenry, pf0, pdensity, pdquer
  real :: pdsigma, pdryvel, pweightmolar, pohreact, pspec_ass, pkao
  real :: pweta_in, pwetb_in, pwetc_in, pwetd_in
  integer :: readerror

  ! declare namelist
  namelist /species_params/ &
   pspecies, pdecay, pweta, pwetb, &
   pweta_in, pwetb_in, pwetc_in, pwetd_in, &
   preldiff, phenry, pf0, pdensity, pdquer, &
   pdsigma, pdryvel, pweightmolar, pohreact, pspec_ass, pkao

  pspecies=" "
  pdecay=-999.9
  pweta=-9.9E-09
  pwetb=0.0
  pweta_in=-9.9E-09
  pwetb_in=-9.9E-09
  pwetc_in=-9.9E-09
  pwetd_in=-9.9E-09
  preldiff=-9.9
  phenry=0.0
  pf0=0.0
  pdensity=-9.9E09
  pdquer=0.0
  pdsigma=0.0
  pdryvel=-9.99
  pohreact=-9.9E-09
  pspec_ass=-9
  pkao=-99.99
  pweightmolar=-789.0 ! read failure indicator value

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
   
  if ((pweightmolar.eq.-789.0).or.(readerror.ne.0)) then ! no namelist found

    readerror=1

    open(unitspecies,file=path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb,status='old',err=998)

    do i=1,6
      read(unitspecies,*)
    end do

    read(unitspecies,'(a10)',end=22) species(pos_spec)
  !  write(*,*) species(pos_spec)
    read(unitspecies,'(f18.1)',end=22) decay(pos_spec)
  !  write(*,*) decay(pos_spec)
    read(unitspecies,'(e18.1)',end=22) weta(pos_spec)
  !  write(*,*) weta(pos_spec)
    read(unitspecies,'(f18.2)',end=22) wetb(pos_spec)
  !  write(*,*) wetb(pos_spec)

  !*** NIK 31.01.2013: including in-cloud scavening parameters
   read(unitspecies,'(e18.1)',end=22) weta_in(pos_spec)
  !  write(*,*) weta_in(pos_spec)
   read(unitspecies,'(f18.2)',end=22) wetb_in(pos_spec)
  !  write(*,*) wetb_in(pos_spec)
   read(unitspecies,'(f18.2)',end=22) wetc_in(pos_spec)
  !  write(*,*) wetc_in(pos_spec)
   read(unitspecies,'(f18.2)',end=22) wetd_in(pos_spec)
  !  write(*,*) wetd_in(pos_spec)

    read(unitspecies,'(f18.1)',end=22) reldiff(pos_spec)
  !  write(*,*) reldiff(pos_spec)
    read(unitspecies,'(e18.1)',end=22) henry(pos_spec)
  !  write(*,*) henry(pos_spec)
    read(unitspecies,'(f18.1)',end=22) f0(pos_spec)
  !  write(*,*) f0(pos_spec)
    read(unitspecies,'(e18.1)',end=22) density(pos_spec)
  !  write(*,*) density(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dquer(pos_spec)
  !  write(*,*) dquer(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dsigma(pos_spec)
  !  write(*,*) dsigma(pos_spec)
    read(unitspecies,'(f18.2)',end=22) dryvel(pos_spec)
  !  write(*,*) dryvel(pos_spec)
    read(unitspecies,'(f18.2)',end=22) weightmolar(pos_spec)
  !  write(*,*) weightmolar(pos_spec)
    read(unitspecies,'(e18.1)',end=22) ohreact(pos_spec)
  !  write(*,*) ohreact(pos_spec)
    read(unitspecies,'(i18)',end=22) spec_ass(pos_spec)
  !  write(*,*) spec_ass(pos_spec)
    read(unitspecies,'(f18.2)',end=22) kao(pos_spec)
  !       write(*,*) kao(pos_spec)

    pspecies=species(pos_spec)
    pdecay=decay(pos_spec)
    pweta=weta(pos_spec)
    pwetb=wetb(pos_spec)
    pweta_in=weta_in(pos_spec)
    pwetb_in=wetb_in(pos_spec)
    pwetc_in=wetc_in(pos_spec)
    pwetd_in=wetd_in(pos_spec)
    preldiff=reldiff(pos_spec)
    phenry=henry(pos_spec)
    pf0=f0(pos_spec)
    pdensity=density(pos_spec)
    pdquer=dquer(pos_spec)
    pdsigma=dsigma(pos_spec)
    pdryvel=dryvel(pos_spec)
    pweightmolar=weightmolar(pos_spec)
    pohreact=ohreact(pos_spec)
    pspec_ass=spec_ass(pos_spec)
    pkao=kao(pos_spec)

  else

    species(pos_spec)=pspecies
    decay(pos_spec)=pdecay
    weta(pos_spec)=pweta
    wetb(pos_spec)=pwetb
    weta_in(pos_spec)=pweta_in
    wetb_in(pos_spec)=pwetb_in
    wetc_in(pos_spec)=pwetc_in
    wetd_in(pos_spec)=pwetd_in
    reldiff(pos_spec)=preldiff
    henry(pos_spec)=phenry
    f0(pos_spec)=pf0
    density(pos_spec)=pdensity
    dquer(pos_spec)=pdquer
    dsigma(pos_spec)=pdsigma
    dryvel(pos_spec)=pdryvel
    weightmolar(pos_spec)=pweightmolar
    ohreact(pos_spec)=pohreact
    spec_ass(pos_spec)=pspec_ass
    kao(pos_spec)=pkao

  endif

  i=pos_spec

  if ((weta(pos_spec).gt.0).and.(henry(pos_spec).le.0)) then
   if (dquer(pos_spec).le.0) goto 996 ! no particle, no henry set
  endif

  if (spec_ass(pos_spec).gt.0) then
    spec_found=.FALSE.
    do j=1,pos_spec-1
      if (spec_ass(pos_spec).eq.specnum(j)) then
        spec_ass(pos_spec)=j
        spec_found=.TRUE.
        ASSSPEC=.TRUE.
      endif
    end do
    if (spec_found.eqv..FALSE.) then
      goto 997
    endif
  endif

  if (dsigma(i).eq.1.) dsigma(i)=1.0001   ! avoid floating exception
  if (dsigma(i).eq.0.) dsigma(i)=1.0001   ! avoid floating exception

  if ((reldiff(i).gt.0.).and.(density(i).gt.0.)) then
    write(*,*) '#### FLEXPART MODEL ERROR! FILE "SPECIES"    ####'
    write(*,*) '#### IS CORRUPT. SPECIES CANNOT BE BOTH      ####'
    write(*,*) '#### PARTICLE AND GAS.                       ####'
    write(*,*) '#### SPECIES NUMBER',aspecnumb
    stop
  endif
20   continue

  if (readerror.ne.0) then ! text format input

    ! Read in daily and day-of-week variation of emissions, if available
    !*******************************************************************
    ! HSO: This is not yet implemented as namelist parameters

    do j=1,24           ! initialize everything to no variation
      area_hour(i,j)=1.
      point_hour(i,j)=1.
    end do
    do j=1,7
      area_dow(i,j)=1.
      point_dow(i,j)=1.
    end do

    read(unitspecies,*,end=22)
    do j=1,24     ! 24 hours, starting with 0-1 local time
      read(unitspecies,*) ihour,area_hour(i,j),point_hour(i,j)
    end do
    read(unitspecies,*)
    do j=1,7      ! 7 days of the week, starting with Monday
      read(unitspecies,*) idow,area_dow(i,j),point_dow(i,j)
    end do

  endif

22 close(unitspecies)

  ! namelist output if requested
  if (nmlout.eqv..true.) then
    open(unitspecies,file=path(2)(1:length(2))//'SPECIES_'//aspecnumb//'.namelist',access='append',status='new',err=1000)
    write(unitspecies,nml=species_params)
    close(unitspecies)
  endif

  return

996   write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### WET DEPOSITION SWITCHED ON, BUT NO HENRYS  #### '
  write(*,*) '#### CONSTANT IS SET                            ####'
  write(*,*) '#### PLEASE MODIFY SPECIES DESCR. FILE!        #### '
  write(*,*) '#####################################################'
  stop


997   write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### THE ASSSOCIATED SPECIES HAS TO BE DEFINED  #### '
  write(*,*) '#### BEFORE THE ONE WHICH POINTS AT IT          #### '
  write(*,*) '#### PLEASE CHANGE ORDER IN RELEASES OR ADD     #### '
  write(*,*) '#### THE ASSOCIATED SPECIES IN RELEASES         #### '
  write(*,*) '#####################################################'
  stop


998   write(*,*) '#####################################################'
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
