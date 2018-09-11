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
  !     HSO, 1 July 2014: Add optional namelist input                          *
  !     PS, 6/2015-9/2018: read regular input with free format                 *
  !       simplify code and rename some variables                              *
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

  integer :: i,kz,istat
  real :: xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4

  ! namelist variables
  integer, parameter :: maxoutlev=500
  integer :: ios
  real,allocatable, dimension (:) :: outheights,outaux
  logical :: lnml

  ! declare namelist
  namelist /nml_outgrid/ &
    outlon0,outlat0, &
    numxgrid,numygrid, &
    dxout,dyout, &
    outheight

! allocate outheights for nml read with max dimension
  allocate(outheights(maxoutlev),outaux(maxoutlev),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate outheights'

  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  outheight(:) = -999. ! initialise for later finding #valid levels
  open(unitoutgrid,file=trim(path(1))//'OUTGRID',status='old',&
    form='formatted',err=999)

  ! try namelist input
  read(unitoutgrid,nml_outgrid,iostat=ios)
  close(unitoutgrid)

  if (ios .eq. 0) then ! namelist works

    lnml = .true.

  else ! read as regular text

    lnml = .false.

    open(unitoutgrid,file=trim(path(1))//'OUTGRID',status='old',err=999)
    call skplin(5,unitoutgrid)

   ! Read horizontal grid specifications
    !****************************************

    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) outlon0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) outlat0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) numxgrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) numygrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) dxout
    call skplin(3,unitoutgrid)
    read(unitoutgrid,*) dyout

  endif ! read OUTGRID file

  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0+real(numxgrid)*dxout
  yr=outlat0+real(numygrid)*dyout
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0+eps.lt.xlon0).or.(outlat0+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) outlon0,outlat0
    write(*,*) xr1,yr1,xlon0,ylat0,xr,yr,dxout,dyout
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') trim(path(1))
    stop
  endif

! Read (if .not. lmnl) and count vertical levels of output grid 
!**************************************************************

  do kz = 1,maxoutlev
    if (lnml) then ! we have read them already
      if (outheight(kz) .lt. 0.) exit ! 1st nondefined level
    else
      call skplin(3,unitoutgrid)
      read(unitoutgrid,*,end=10) outheight(kz)
    endif
  end do
10 continue  

  numzgrid = kz - 1 ! number of outgrid levels

! allocate the required length only, shuffle data  
  outaux = outheight ! shuffle

  deallocate(outheights)
  allocate(outheight(numzgrid),outheighthalf(numzgrid),stat=istat)
  if (istat .ne. 0) then
    write(*,*) 'ERROR: could not allocate outheight and outheighthalf'
    stop 'readoutgrid error'
  endif

  outheight=outaux(1:numzgrid) ! shuffle back
  deallocate (outaux) 

! write outgrid file in namelist format to output directory if requested
  if (nmlout) then
    open(unitoutgrid,file=trim(path(2))//'OUTGRID.namelist',err=1000)
    write(unitoutgrid,nml=nml_outgrid)
    close(unitoutgrid)
  endif

  ! Check whether vertical levels are specified in ascending order
  !***************************************************************

  do kz=2,numzgrid
    if (outheight(kz) .le. outheight(kz-1)) goto 998
  end do

  ! Determine the half levels, i.e. middle levels of the output grid
  !*****************************************************************

  outheighthalf(1) = 0.5*outheight(1)
  do kz = 2,numzgrid
    outheighthalf(kz) = 0.5*(outheight(kz-1)+outheight(kz))
  end do

  xoutshift=xlon0-outlon0
  youtshift=ylat0-outlat0

  allocate(oroout(0:numxgrid-1,0:numygrid-1),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate oroout'
  allocate(area(0:numxgrid-1,0:numygrid-1),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate area'
  allocate(volume(0:numxgrid-1,0:numygrid-1,numzgrid),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate volume'
  allocate(areaeast(0:numxgrid-1,0:numygrid-1,numzgrid),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate areaeast'
  allocate(areanorth(0:numxgrid-1,0:numygrid-1,numzgrid),stat=istat)
  if (istat .ne. 0) write(*,*)'ERROR: could not allocate areanorth'
  
  return

998 continue
  write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### '
  write(*,*) ' #### OF OUTPUT LEVELS NOT INCREASING AT LEVEL#### '
  write(*,*) ' #### ',kz,'                                  #### '
  write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### '
  STOP 'readoutgrid error'

999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') trim(path(1))
  stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') trim(path(2))
  stop

end subroutine readoutgrid
