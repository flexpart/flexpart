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

subroutine writeheader_txt

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !     modified IP 2013 for text output                                       *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

!  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  integer :: jjjjmmdd,ihmmss,i,j
  real :: xp1,yp1,xp2,yp2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(2)(1:length(2))//'header_txt', &
       form='formatted',err=998)
  open(unitheader_txt,file=path(2)(1:length(2))//'header_txt_releases', &
       form='formatted',err=998)


  ! Write the header information
  !*****************************
  
  write(unitheader,*) '# ibdate,ibtime, iedate, ietime, flexversion'
  write(unitheader,*) ibdate, ibtime, iedate, ietime, trim(flexversion) !  'FLEXPART V9.0'
  !if (ldirect.eq.1) then
  !  write(unitheader,*) ibdate,ibtime,trim(flexversion) !  'FLEXPART V9.0'
  !else
  !  write(unitheader,*) iedate,ietime,trim(flexversion) ! 'FLEXPART V9.0'
  !endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************
  
  write(unitheader,*) '# interval, averaging time, sampling time'
  write(unitheader,*) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************
  
  write(unitheader,*) '# information on grid setup    '
  write(unitheader,*) '#outlon0,outlat0,numxgrid,numygrid,dxout,dyout'
  write(unitheader,*) outlon0,outlat0,numxgrid,numygrid, &
       dxout,dyout  
  write(unitheader,*) '# numzgrid, outheight(.) '
  write(unitheader,*) numzgrid,(outheight(i),i=1,numzgrid)

  write(unitheader,*) '# jjjjmmdd,ihmmss'
  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader,*) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the vertical dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader,*) '# information on species'
  write(unitheader,*) '# 3*nspec,maxpointspec_act'
  write(unitheader,*) 3*nspec,maxpointspec_act
  write(unitheader,*) '# for nspec:'
  write(unitheader,*) '# 1, WD_ '
  write(unitheader,*) '# 1, DD_ '
  write(unitheader,*) '# numzgrid,species'
  do i=1,nspec
    write(unitheader,*) 1,'WD_'//species(i)(1:7)
    write(unitheader,*) 1,'DD_'//species(i)(1:7)
    write(unitheader,*) numzgrid,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************


  write(unitheader_txt,*) '# information on release points'
  write(unitheader_txt,*) '# numpoint'
  write(unitheader_txt,*) numpoint
  write(unitheader_txt,*) '# for numpoint:'
  do i=1,numpoint
    write(unitheader_txt,*) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader_txt,*) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader_txt,*) npart(i),1
    if (numpoint.le.1000) then
      write(unitheader_txt,*) compoint(i)
    else
      write(unitheader_txt,*) compoint(1001)
    endif
    do j=1,nspec
      write(unitheader_txt,*) xmass(i,j)
      write(unitheader_txt,*) xmass(i,j)
      write(unitheader_txt,*) xmass(i,j)
    end do
  end do

  ! Write information on model switches
  !*****************************************

  write(unitheader,*) '# information on model switches'
  write(unitheader,*) '# method,lsubgrid,lconvection, ind_source,ind_receptor'
  write(unitheader,*) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************
  
  write(unitheader,*) '# information on age class     '
  write(unitheader,*) nageclass,(lage(i),i=1,nageclass)


  !Do not write topography to text output file. Keep it on the binary one
  !********************************

  !do ix=0,numxgrid-1
  !  write(unitheader,*) (oroout(ix,jy),jy=0,numygrid-1)
  !end do


  


  close(unitheader)
  close(unitheader_txt)


!  open(unitheader,file=path(2)(1:length(2))//'header_nml', &
!        form='formatted',err=998)
!  write(unitheader,NML=COMMAND) 
!  close(unitheader)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header_txt'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader_txt
