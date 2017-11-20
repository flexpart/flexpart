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

subroutine gridcheck_nests

  !*****************************************************************************
  !                                                                            *
  !     This routine checks the grid specification for the nested model        *
  !     domains. It is similar to subroutine gridcheck, which checks the       *
  !     mother domain.                                                         *
  !                                                                            *
  !     Authors: A. Stohl, G. Wotawa                                           *
  !                                                                            *
  !     8 February 1999                                                        *
  !                                                                            *
  !*****************************************************************************
  !  CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with ECMWF grib_api    *
  !  CHANGE: 03/12/2008, Harald Sodemann, change to f90 grib_api               *
  !                                                                            *
  !   Implementation of the Vtables approach                                   *
  !   D. Morton, D. Arnold 17.11.2017                                          *
  !     - Inclusion of specific code and usage of class_vtable                 *
  !                                                                            *
  !*****************************************************************************

  use grib_api
  use par_mod
  use com_mod
  use class_vtable

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  integer :: parID !added by mc for making it consistent with new gridcheck.f90
  integer :: gotGrib
  !HSO  end
  integer :: i,j,k,l,ifn,ifield,iumax,iwmax,numskip,nlev_ecn
  integer :: nuvzn,nwzn
  real :: akmn(nwzmax),bkmn(nwzmax),akzn(nuvzmax),bkzn(nuvzmax)
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  real :: conversion_factor !added by mc to make it consistent with new gridchek.f90

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

!!! DJM -- isec1() no longer needed --  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  integer :: isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec2(60+2*nuvzmax),zsec4(jpunp)

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheck_nests'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Vtable related variables
  !
  !  Path to Vtable - current implementation assumes it's in cwd, named
  !  "Vtable"
  CHARACTER(LEN=255), PARAMETER :: VTABLE_PATH = "Vtable"
  CHARACTER(LEN=15) :: fpname      ! stores FLEXPART name for curr grib mesg.
  TYPE(Vtable) :: my_vtable    ! unallocated
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!  DJM
  INTEGER current_grib_level   ! this was isec1(8) in previous versions


  xresoln(0)=1.       ! resolution enhancement for mother grid
  yresoln(0)=1.       ! resolution enhancement for mother grid

  ! Loop about all nesting levels
  !******************************

  do l=1,numbnests

    iumax=0
    iwmax=0

    if(ideltas.gt.0) then
      ifn=1
    else
      ifn=numbwf
    endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Vtable code
  PRINT *, 'Loading Vtable: ', VTABLE_PATH
  call vtable_load_by_name(VTABLE_PATH, my_vtable)
  !! Debugging tool
  !PRINT *, 'Dump of Vtable...'
  !call vtable_dump_records(my_vtable)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  ! This is diagnostic/debugging code, and will normally be commented out.
  ! It's purpose is to look at the provided grib file and produce an
  ! inventory of the FP-related messages, relative to the Vtable that's
  ! already been open.

  !CALL vtable_gribfile_inventory(path(3)(1:length(3)) // trim(wfname(ifn)), &
!&                                my_vtable)
!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
  ifile=0
  igrib=0
  iret=0

5   call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,ifn)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  !call grib_multi_support_on

  gotGrib=0
  ifield=0
10   ifield=ifield+1

  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE)  then
    goto 30    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  ! Get the fpname
  fpname = vtable_get_fpname(igrib, my_vtable)
  !print *, 'fpname: ', trim(fpname)


  !!!!!!!!!!!!!!!!!!!  VTABLE code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !!  DJM - get the current_grib_level (in previous code it was isec1(8))
  !!  It's the same in both GRIB1 and GRIB2
  call grib_get_int(igrib,'level',current_grib_level,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)



  !get the size and data of the values array
  !!!! -- original statement -- if (isec1(6).ne.-1) then
  ! 'NOFP' is the fpname for a grib message not recognized by FLEXPART
  IF (TRIM(fpname) .NE. 'NOFP') THEN
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  !HSO  get the required fields from section 2 in a gribex compatible manner
  if (ifield.eq.1) then
    call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
         isec2(2),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
         isec2(3),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfVerticalCoordinateValues', &
         isec2(12),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  !HSO    get the size and data of the vertical coordinate array
    call grib_get_real4_array(igrib,'pv',zsec2,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    nxn(l)=isec2(2)
    nyn(l)=isec2(3)
    nlev_ecn=isec2(12)/2-1
  endif ! ifield

  if (nxn(l).gt.nxmaxn) then
  write(*,*) 'FLEXPART error: Too many grid points in x direction.'
  write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
  write(*,*) 'for nesting level ',l
  write(*,*) 'Or change parameter settings in file par_mod.'
  write(*,*) nxn(l),nxmaxn
  stop
  endif

  if (nyn(l).gt.nymaxn) then
  write(*,*) 'FLEXPART error: Too many grid points in y direction.'
  write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
  write(*,*) 'for nesting level ',l
  write(*,*) 'Or change parameter settings in file par_mod.'
  write(*,*) nyn(l),nymaxn
  stop
  endif

  !HSO  get the second part of the grid dimensions only from GRiB1 messages
  if (TRIM(fpname) .EQ. 'T2'.and. (gotGrib.eq.0)) then
!!!!!  if (isec1(6) .eq. 167 .and. (gotGrib.eq.0)) then !added by mc to make it consistent with new gridchek.f90 note that gotGrid must be changed in gotGrib!!
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', & !comment by mc: note that this was in the (if (ifield.eq.1) ..end above in gridchek.f90 see line 257
         xaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'longitudeOfLastGridPointInDegrees', &
         xaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfFirstGridPointInDegrees', &
         yaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    xaux1=xaux1in
    xaux2=xaux2in
    yaux1=yaux1in
    yaux2=yaux2in
    if(xaux1.gt.180.) xaux1=xaux1-360.0
    if(xaux2.gt.180.) xaux2=xaux2-360.0
    if(xaux1.lt.-180.) xaux1=xaux1+360.0
    if(xaux2.lt.-180.) xaux2=xaux2+360.0
    if (xaux2.lt.xaux1) xaux2=xaux2+360.0
    xlon0n(l)=xaux1
    ylat0n(l)=yaux1
    dxn(l)=(xaux2-xaux1)/real(nxn(l)-1)
    dyn(l)=(yaux2-yaux1)/real(nyn(l)-1)
    gotGrib=1 !commetn by mc note tahthere gotGRIB is used instead of gotGrid!!!
  endif ! ifield.eq.1

!!!!!!!!! DJM - orig -  k=isec1(8)
k = current_grib_level


!!!!!! DJM - orig -   if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
  IF (TRIM(fpname) .EQ. 'UU') iumax=max(iumax,nlev_ec-k+1)


!!!!!! DJM - orig -   if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)
  IF (TRIM(fpname) .EQ. 'ETADOT') iwmax=max(iwmax,nlev_ec-k+1)

  !!!!!!! DJM - orig - if(isec1(6).eq.129) then
  IF (TRIM(fpname) .EQ. 'ORO') THEN
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        oron(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif

  !!!!!!! DJM - orig - if(isec1(6).eq.172) then
  IF (TRIM(fpname) .EQ. 'LSM') THEN
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        lsmn(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif

  !!!!!! DJM - orig - if(isec1(6).eq.160) then
  IF (TRIM(fpname) .EQ. 'EXCESSORO') THEN
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        excessoron(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif

  call grib_release(igrib)
  goto 10                 !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

30   call grib_close_file(ifile)

  !error message if no fields found with correct first longitude in it
  if (gotGrib.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  nuvzn=iumax
  nwzn=iwmax
  if(nuvzn.eq.nlev_ec) nwzn=nlev_ecn+1

  if ((nuvzn.gt.nuvzmax).or.(nwzn.gt.nwzmax)) then
  write(*,*) 'FLEXPART error: Nested wind fields have too many'// &
       'vertical levels.'
  write(*,*) 'Problem was encountered for nesting level ',l
  stop
  endif


  ! Output of grid info
  !********************

  write(*,'(a,i2,a)') ' Nested domain ',l,':'
  write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Longitude range: ', &
       xlon0n(l),' to ',xlon0n(l)+(nxn(l)-1)*dxn(l), &
       '   Grid distance: ',dxn(l)
  write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Latitude range : ', &
       ylat0n(l),' to ',ylat0n(l)+(nyn(l)-1)*dyn(l), &
       '   Grid distance: ',dyn(l)
  write(*,*)

  ! Determine, how much the resolutions in the nests are enhanced as
  ! compared to the mother grid
  !*****************************************************************

  xresoln(l)=dx/dxn(l)
  yresoln(l)=dy/dyn(l)

  ! Determine the mother grid coordinates of the corner points of the
  ! nested grids
  ! Convert first to geographical coordinates, then to grid coordinates
  !********************************************************************

  xaux1=xlon0n(l)
  xaux2=xlon0n(l)+real(nxn(l)-1)*dxn(l)
  yaux1=ylat0n(l)
  yaux2=ylat0n(l)+real(nyn(l)-1)*dyn(l)

  xln(l)=(xaux1-xlon0)/dx
  xrn(l)=(xaux2-xlon0)/dx
  yln(l)=(yaux1-ylat0)/dy
  yrn(l)=(yaux2-ylat0)/dy


  if ((xln(l).lt.0.).or.(yln(l).lt.0.).or. &
       (xrn(l).gt.real(nxmin1)).or.(yrn(l).gt.real(nymin1))) then
    write(*,*) 'Nested domain does not fit into mother domain'
    write(*,*) 'For global mother domain fields, you can shift'
    write(*,*) 'shift the mother domain into x-direction'
    write(*,*) 'by setting nxshift (file par_mod) to a'
    write(*,*) 'positive value. Execution is terminated.'
    stop
  endif


  ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
  ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

  numskip=nlev_ecn-nuvzn ! number of ecmwf model layers not used by FLEXPART
  do i=1,nwzn
    j=numskip+i
    k=nlev_ecn+1+numskip+i
    akmn(nwzn-i+1)=zsec2(j)
    bkmn(nwzn-i+1)=zsec2(k)
  end do

  !
  ! CALCULATION OF AKZ, BKZ
  ! AKZ,BKZ: model discretization parameters at the center of each model
  !     layer
  !
  ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
  ! i.e. ground level
  !*****************************************************************************

    akzn(1)=0.
    bkzn(1)=1.0
    do i=1,nuvzn
      akzn(i+1)=0.5*(akmn(i+1)+akmn(i))
      bkzn(i+1)=0.5*(bkmn(i+1)+bkmn(i))
    end do
    nuvzn=nuvzn+1

  ! Check, whether the heights of the model levels of the nested
  ! wind fields are consistent with those of the mother domain.
  ! If not, terminate model run.
  !*************************************************************

    do i=1,nuvz
      if ((akzn(i).ne.akz(i)).or.(bkzn(i).ne.bkz(i))) then
  write(*,*) 'FLEXPART error: The wind fields of nesting level',l
  write(*,*) 'are not consistent with the mother domain:'
  write(*,*) 'Differences in vertical levels detected.'
        stop
      endif
    end do

    do i=1,nwz
      if ((akmn(i).ne.akm(i)).or.(bkmn(i).ne.bkm(i))) then
  write(*,*) 'FLEXPART error: The wind fields of nesting level',l
  write(*,*) 'are not consistent with the mother domain:'
  write(*,*) 'Differences in vertical levels detected.'
        stop
      endif
    end do

  end do

  return

999   write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '                FLEXPART SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfnamen(l,ifn)
  write(*,*) ' FOR NESTING LEVEL ',k
  write(*,*) ' ###########################################'// &
       '###### '
  stop

end subroutine gridcheck_nests
