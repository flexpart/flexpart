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

subroutine readageclasses

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the age classes to be used for the current model    *
  !     run.                                                                   *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !     20 March 2000                                                          *
  !
  !     HSO, 1 July 2014                                                       *
  !     Added optional namelist input                                          *
  !                                                                            *
  !     PS, 6/2015-9/2018 some variable names changed as in readreleases.f90   *
  !       catch nageclass>maxage properly                                      *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: i

  ! namelist aux variables
  integer :: ios

  ! namelist declaration
  namelist /nml_ageclass/ &
    nageclass, &
    lage

  nageclass=-1 ! preset to negative value to identify failed namelist input

  ! If age spectra calculation is switched off, set number of age classes
  ! to 1 and maximum age to a large number
  !**********************************************************************

  if (lagespectra.ne.1) then
    nageclass=1
    lage(nageclass)=999999999
    return
  endif

  ! If age spectra claculation is switched on,
  ! open the AGECLASSSES file and read user options
  !************************************************

  open(unitageclasses,file=path(1)(1:length(1))//'AGECLASSES', &
    form='formatted',status='old',err=999)

  ! try to read in as a namelist
  read(unitageclasses, nml_ageclass, iostat=ios)
  close(unitageclasses)

  if (ios .ne. 0) then ! failed to read nml, assume simple text
    open(unitageclasses,file=path(1)(1:length(1))// &
      'AGECLASSES',status='old',err=999)
    call skplin(13,unitageclasses)
    read(unitageclasses,*) nageclass ! number of classes
    if (nageclass.gt.maxageclass) goto 1001
    do i=1,nageclass
      read(unitageclasses,*) lage(i) ! max age per classes
    end do
    close(unitageclasses)
  endif

  if (nageclass.lt.1) goto 1002

  ! write ageclasses file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitageclasses,file=path(2)(1:length(2))//'AGECLASSES.namelist', &
      err=1000)
    write(unitageclasses,nml=nml_ageclass)
    close(unitageclasses)
  endif

  if (lage(1).le.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AGE OF FIRST      #### '
    write(*,*) ' #### CLASS MUST BE GREATER THAN ZERO. CHANGE #### '
    write(*,*) ' #### SETTINGS IN FILE AGECLASSES.            #### '
    stop
  endif

  do i=2,nageclass
    if (lage(i).le.lage(i-1)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! AGE CLASSES     #### '
      write(*,*) ' #### MUST BE GIVEN IN TEMPORAL ORDER.      #### '
      write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES.   #### '
      stop
    endif
  end do

  return

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "AGECLASSES" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') trim(path(1))
  stop

1000  write(*,*) ' #### FLEXPART MODEL ERROR! FILE "AGECLASSES" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop

1001 continue
    write(*,*) ' #### FLEXPART MODEL ERROR! NUMBER OF AGE     #### '
    write(*,*) ' #### CLASSES GREATER THAN MAXIMUM ALLOWED.   #### '
    write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES OR   #### '
    write(*,*) ' #### RECOMPILE WITH LARGER MAXAGECLASS IN    #### '
    write(*,*) ' #### FILE PAR_MOD.                           #### '
    stop

1002 continue
    write(*,*) ' #### FLEXPART MODEL ERROR! NUMBER OF AGE     #### '
    write(*,*) ' #### CLASSES < 1                             #### '
    write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES      #### '
    stop

end subroutine readageclasses
