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
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decaytime(maxtable)  half time for radiological decay                      *
  ! specname(maxtable)   names of chemical species, radionuclides              *
  ! wetscava, wetscavb   Parameters for determining scavenging coefficient     *
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

  ! Open the SPECIES file and read species names and properties
  !************************************************************
  specnum(pos_spec)=id_spec
  write(aspecnumb,'(i3.3)') specnum(pos_spec)
  open(unitspecies,file= &
       path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb,status='old', &
       err=998)
  !write(*,*) 'reading SPECIES',specnum(pos_spec)

  ASSSPEC=.FALSE.

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


  ! Read in daily and day-of-week variation of emissions, if available
  !*******************************************************************

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

22   close(unitspecies)

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


end subroutine readspecies
