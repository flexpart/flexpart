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

subroutine openreceptors

  !*****************************************************************************
  !                                                                            *
  !  This routine opens the receptor output files and writes out the receptor  *
  !  names and the receptor locations. The receptor output files are not       *
  !  closed, but kept open throughout the simulation. Concentrations are       *
  !  continuously  dumped to these files.                                      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! numreceptor            actual number of receptor points specified          *
  ! receptornames          names of the receptor points                        *
  ! xreceptor,yreceptor    coordinates of the receptor points                  *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: j

  ! Open output file for receptor points and write out a short header
  ! containing receptor names and locations
  !******************************************************************

  if (numreceptor.ge.1) then           ! do it only if receptors are specified

  ! Concentration output
  !*********************

    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      open(unitoutrecept,file=path(2)(1:length(2))//'receptor_conc', &
           form='unformatted',err=997)
      write(unitoutrecept) (receptorname(j),j=1,numreceptor)
      write(unitoutrecept) (xreceptor(j)*dx+xlon0, &
           yreceptor(j)*dy+ylat0,j=1,numreceptor)
    endif

  ! Mixing ratio output
  !********************

    if ((iout.eq.2).or.(iout.eq.3)) then
      open(unitoutreceptppt,file=path(2)(1:length(2))//'receptor_pptv', &
           form='unformatted',err=998)
      write(unitoutreceptppt) (receptorname(j),j=1,numreceptor)
      write(unitoutreceptppt) (xreceptor(j)*dx+xlon0, &
           yreceptor(j)*dy+ylat0,j=1,numreceptor)
    endif
  endif

  return


997   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
  write(*,*) ' ####              receptor_conc               #### '
  write(*,*) ' #### CANNOT BE OPENED.                        #### '
  stop

998   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
  write(*,*) ' ####              receptor_pptv               #### '
  write(*,*) ' #### CANNOT BE OPENED.                        #### '
  stop

end subroutine openreceptors
