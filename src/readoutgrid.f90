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

subroutine readoutgrid

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the output grid.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     4 June 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dxout,dyout          grid distance                                         *
  ! numxgrid,numygrid,numzgrid    grid dimensions                              *
  ! outlon0,outlat0      lower left corner of grid                             *
  ! outheight(maxzgrid)  height levels of output grid [m]                      *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitoutgrid          unit connected to file OUTGRID                        *
  !                                                                            *
  !*****************************************************************************

  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i,j,stat
  real :: outhelp,xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4



  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID',status='old', &
       err=999)


  call skplin(5,unitoutgrid)


  ! 1.  Read horizontal grid specifications
  !****************************************

  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,f11.4)') outlon0
  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,f11.4)') outlat0
  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,i5)') numxgrid
  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,i5)') numygrid
  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,f12.5)') dxout
  call skplin(3,unitoutgrid)
  read(unitoutgrid,'(4x,f12.5)') dyout


  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0+real(numxgrid)*dxout
  yr=outlat0+real(numygrid)*dyout
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0+eps.lt.xlon0).or.(outlat0+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) 'outlon0,outlat0:'
    write(*,*) outlon0,outlat0
    write(*,*) 'xr1,yr1,xlon0,ylat0,xr,yr,dxout,dyout:'
    write(*,*) xr1,yr1,xlon0,ylat0,xr,yr,dxout,dyout
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') path(1)(1:length(1))
    stop
  endif

  ! 2. Count Vertical levels of output grid
  !****************************************
  j=0
100   j=j+1
    do i=1,3
      read(unitoutgrid,*,end=99)
    end do
    read(unitoutgrid,'(4x,f7.1)',end=99) outhelp
    if (outhelp.eq.0.) goto 99
    goto 100
99   numzgrid=j-1

    allocate(outheight(numzgrid) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
    allocate(outheighthalf(numzgrid) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'


  rewind(unitoutgrid)
  call skplin(29,unitoutgrid)

  ! 2. Vertical levels of output grid
  !**********************************

  j=0
1000   j=j+1
    do i=1,3
      read(unitoutgrid,*,end=990)
    end do
    read(unitoutgrid,'(4x,f7.1)',end=990) outhelp
    if (outhelp.eq.0.) goto 99
    outheight(j)=outhelp
    goto 1000
990   numzgrid=j-1


  ! Check whether vertical levels are specified in ascending order
  !***************************************************************

  do j=2,numzgrid
    if (outheight(j).le.outheight(j-1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### '
    write(*,*) ' #### OF OUTPUT LEVELS IS CORRUPT AT LEVEL    #### '
    write(*,*) ' #### ',j,'                              #### '
    write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### '
    endif
  end do

  ! Determine the half levels, i.e. middle levels of the output grid
  !*****************************************************************

  outheighthalf(1)=outheight(1)/2.
  do j=2,numzgrid
    outheighthalf(j)=(outheight(j-1)+outheight(j))/2.
  end do


  xoutshift=xlon0-outlon0
  youtshift=ylat0-outlat0
  close(unitoutgrid)

    allocate(oroout(0:numxgrid-1,0:numygrid-1) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
    allocate(area(0:numxgrid-1,0:numygrid-1) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
    allocate(volume(0:numxgrid-1,0:numygrid-1,numzgrid) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
    allocate(areaeast(0:numxgrid-1,0:numygrid-1,numzgrid) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
    allocate(areanorth(0:numxgrid-1,0:numygrid-1,numzgrid) &
         ,stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outh'
  return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,*) ' #### xxx/flexpart/options                    #### '
  stop

end subroutine readoutgrid
