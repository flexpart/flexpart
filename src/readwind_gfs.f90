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

subroutine readwind_gfs(indj,n,uuh,vvh,wwh)

  !***********************************************************************
  !*                                                                     *
  !*             TRAJECTORY MODEL SUBROUTINE READWIND                    *
  !*                                                                     *
  !***********************************************************************
  !*                                                                     *
  !*             AUTHOR:      G. WOTAWA                                  *
  !*             DATE:        1997-08-05                                 *
  !*             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
  !*             CHANGE: 01/02/2001, Bernd C. Krueger, Variables tth and *
  !*                     qvh (on eta coordinates) in common block        *
  !*             CHANGE: 16/11/2005, Caroline Forster, GFS data          *
  !*             CHANGE: 11/01/2008, Harald Sodemann, Input of GRIB1/2   *
  !*                     data with the ECMWF grib_api library            *
  !*             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
  !*                                 ECMWF grib_api                      *
  !                                                                      *
  !   Unified ECMWF and GFS builds                                       *
  !   Marian Harustak, 12.5.2017                                         *
  !     - Renamed routine from readwind to readwind_gfs                  *
  !                                                                      *
  !   Implementation of the Vtables approach                             *
  !   D. Morton, D. Arnold 17.11.2017                                    *
  !     - Inclusion of specific code and usage of class_vtable           *
  !                                                                      *
  !*                                                                     *
  !***********************************************************************
  !*                                                                     *
  !* DESCRIPTION:                                                        *
  !*                                                                     *
  !* READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
  !* INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
  !*                                                                     *
  !* INPUT:                                                              *
  !* indj               indicates number of the wind field to be read in *
  !* n                  temporal index for meteorological fields (1 to 3)*
  !*                                                                     *
  !* IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
  !*                                                                     *
  !* wfname             File name of data to be read in                  *
  !* nx,ny,nuvz,nwz     expected field dimensions                        *
  !* nlev_ec            number of vertical levels ecmwf model            *
  !* uu,vv,ww           wind fields                                      *
  !* tt,qv              temperature and specific humidity                *
  !* ps                 surface pressure                                 *
  !*                                                                     *
  !***********************************************************************

  use grib_api
  use par_mod
  use com_mod
  use class_vtable

  implicit none

  !HSO  new parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  !HSO end edits
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  integer :: ii,indj,i,j,k,n,levdiff2,ifield,iumax,iwmax

  ! NCEP
  integer :: numpt,numpu,numpv,numpw,numprh
  real :: help, temp, ew
  real :: elev
  real :: ulev1(0:nxmax-1,0:nymax-1),vlev1(0:nxmax-1,0:nymax-1)
  real :: tlev1(0:nxmax-1,0:nymax-1)
  real :: qvh2(0:nxmax-1,0:nymax-1)

  integer :: i179,i180,i181

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING
  !HSO kept isec1, isec2 and zsec4 for consistency with gribex GRIB input



  !!!!!  DJM - ultimately I will remove isec1
  !!!!!integer :: isec1(8),isec2(3)
  integer :: isec2(3)
  real(kind=4) :: zsec4(jpunp)
  real(kind=4) :: xaux,yaux,xaux0,yaux0
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real(kind=4) :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
  real :: plev1,hlev1,ff10m,fflev1

  logical :: hflswitch,strswitch

  !HSO  for grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind_gfs'


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




  hflswitch=.false.
  strswitch=.false.
  levdiff2=nlev_ec-nwz+1
  iumax=0
  iwmax=0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Vtable code
  PRINT *, 'Loading Vtable: ', VTABLE_PATH
  call vtable_load_by_name(VTABLE_PATH, my_vtable)
  !! Debugging tool
!  PRINT *, 'Dump of Vtable...'
!  call vtable_dump_records(my_vtable)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  ! This is diagnostic/debugging code, and will normally be commented out.
  ! It's purpose is to look at the provided grib file and produce an
  ! inventory of the FP-related messages, relative to the Vtable that's
  ! already been open.

!  CALL vtable_gribfile_inventory(path(3)(1:length(3)) // trim(wfname(indj)), &
! &                                my_vtable)

  !!!!!!!!!!!!!!!!!!!  VTABLE code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  ! OPENING OF DATA FILE (GRIB CODE)

  !HSO
5   call grib_open_file(ifile,path(3)(1:length(3)) &
         //trim(wfname(indj)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  call grib_multi_support_on

  numpt=0
  numpu=0
  numpv=0
  numpw=0
  numprh=0
  ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE)  then
    goto 50    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
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
!  call grib_check(iret,gribFunction,gribErrorMsg)


  if (gribVer.eq.1) then ! GRIB Edition 1

      call grib_get_int(igrib,'level',current_grib_level,iret)
      call grib_check(iret,gribFunction,gribErrorMsg)

      !!! Added by DJM 2017-08-08 - if this is GRIB1 we assume that
      !!! level units are hPa and need to be multiplied by 100 for Pa
      !!! As currently implemented, ALL levels (even non isobaric) are transformed.
      !!! I "think" that's OK, because I don't think levels are used for any of the
      !!! 2D fields.  But, I need to double check.  I could always make this
      !!! conditional on whether the level type is isobaric (type 105 I think)
      current_grib_level = current_grib_level*100.0

  else ! GRIB Edition 2

      call grib_get_int(igrib,'scaledValueOfFirstFixedSurface', &
           current_grib_level,iret)
      call grib_check(iret,gribFunction,gribErrorMsg)

     !!! Added by DJM 2017-08-08 - if this is GRIB2 we assume that
      !!! level units are Pa and don't need to be modified

  endif ! gribVer






!!!!!!  DJM - orig -   if (isec1(6).ne.-1) then
  IF (TRIM(fpname) .NE. 'NOFP') THEN
  !  get the size and data of the values array
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  if(ifield.eq.1) then

  !get the required fields from section 2
  !store compatible to gribex input
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
!  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
!  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
       xauxin,iret)
!  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
       yauxin,iret)
!  call grib_check(iret,gribFunction,gribErrorMsg)
  xaux=xauxin+real(nxshift)*dx
  yaux=yauxin

  ! CHECK GRID SPECIFICATIONS

    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    if(xaux.eq.0.) xaux=-179.0     ! NCEP DATA
    xaux0=xlon0
    yaux0=ylat0
    if(xaux.lt.0.) xaux=xaux+360.
    if(yaux.lt.0.) yaux=yaux+360.
    if(xaux0.lt.0.) xaux0=xaux0+360.
    if(yaux0.lt.0.) yaux0=yaux0+360.
    if(abs(xaux-xaux0).gt.eps) &
         stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
    if(abs(yaux-yaux0).gt.eps) &
         stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
  endif
  !HSO end of edits

  i179=nint(179./dx)
  if (dx.lt.0.7) then
    i180=nint(180./dx)+1    ! 0.5 deg data
  else
    i180=nint(179./dx)+1    ! 1 deg data
  endif
  i181=i180+1

  !!!!!! DJM - orig - if (isec1(6).ne.-1) then
  IF (TRIM(fpname) .NE. 'NOFP') THEN

      IF (TRIM(fpname) .EQ. 'TT') THEN

          DO ii=1,nuvz
              if (current_grib_level .EQ. akz(ii)) numpt=ii
          END DO

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      tth(i179+i,j,numpt,n)=help
                  else
                      tth(i-i181,j,numpt,n)=help
                  endif
              END DO
          END DO


      ELSEIF (TRIM(fpname) .EQ. 'UU') THEN

          DO ii=1,nuvz
              if (current_grib_level .EQ. akz(ii)) numpu=ii
          END DO

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      uuh(i179+i,j,numpu)=help
                  else
                      uuh(i-i181,j,numpu)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'VV') THEN

          DO ii=1,nuvz
              if (current_grib_level .EQ. akz(ii)) numpv=ii
          END DO

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      vvh(i179+i,j,numpv)=help
                  else
                      vvh(i-i181,j,numpv)=help
                  endif
              END DO
          END DO


      ELSEIF (TRIM(fpname) .EQ. 'RH') THEN

          DO ii=1,nuvz
              if (current_grib_level .EQ. akz(ii)) numprh=ii
          END DO

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      qvh(i179+i,j,numprh,n)=help
                  else
                      qvh(i-i181,j,numprh,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'WW') THEN

          DO ii=1,nuvz
              if (current_grib_level .EQ. akz(ii)) numpw=ii
          END DO

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      wwh(i179+i,j,numpw)=help
                  else
                      wwh(i-i181,j,numpw)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'PS') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      ps(i179+i,j,1,n)=help
                  else
                      ps(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'SD') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      sd(i179+i,j,1,n)=help
                  else
                      sd(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'SLP') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      msl(i179+i,j,1,n)=help
                  else
                      msl(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'TCC') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      tcc(i179+i,j,1,n)=help
                  else
                      tcc(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'U10') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      u10(i179+i,j,1,n)=help
                  else
                      u10(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'V10') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      v10(i179+i,j,1,n)=help
                  else
                      v10(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'T2') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      tt2(i179+i,j,1,n)=help
                  else
                      tt2(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'TD2') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      td2(i179+i,j,1,n)=help
                  else
                      td2(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'LSPREC') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      lsprec(i179+i,j,1,n)=help
                  else
                      lsprec(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'CONVPREC') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      convprec(i179+i,j,1,n)=help
                  else
                      convprec(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'ORO') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      oro(i179+i,j)=help
                      excessoro(i179+i,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                  else
                      oro(i-i181,j)=help
                      excessoro(i-i181,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'LSM') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      lsm(i179+i,j)=help
                  else
                      lsm(i-i181,j)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'HMIX') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      hmix(i179+i,j,1,n)=help
                  else
                      hmix(i-i181,j,1,n)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'RH2') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      qvh2(i179+i,j)=help
                  else
                      qvh2(i-i181,j)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'TSIG1') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      tlev1(i179+i,j)=help
                  else
                      tlev1(i-i181,j)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'USIG1') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      ulev1(i179+i,j)=help
                  else
                      ulev1(i-i181,j)=help
                  endif
              END DO
          END DO

      ELSEIF (TRIM(fpname) .EQ. 'VSIG1') THEN

          DO j=0,nymin1
              DO i=0,nxfield-1
                  help=zsec4(nxfield*(ny-j-1)+i+1)
                  if(i.le.i180) then
                      vlev1(i179+i,j)=help
                  else
                      vlev1(i-i181,j)=help
                  endif
              END DO
          END DO












      END IF


!!!!!!!  BEGIN OLD LOOP STRUCTURE - WILL EVENTUALLY BE DELETED
#ifdef NO_LONGER_NEEDED
  do j=0,nymin1
    do i=0,nxfield-1

!!!!!! DJM - orig -       if((isec1(6).eq.011).and.(isec1(7).eq.100)) then


      IF (TRIM(fpname) .EQ. 'TT') THEN
  ! TEMPERATURE
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              !!!!!! DJM - orig - if ((isec1(8)*100.0).eq.akz(ii)) numpt=ii
              if (current_grib_level .EQ. akz(ii)) numpt=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tth(i179+i,j,numpt,n)=help
        else
          tth(i-i181,j,numpt,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.033).and.(isec1(7).eq.100)) then
      IF (TRIM(fpname) .EQ. 'UU') THEN
  ! U VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              !!!!!! DJM - orig if ((isec1(8)*100.0).eq.akz(ii)) numpu=ii
              if (current_grib_level .EQ. akz(ii)) numpu=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          uuh(i179+i,j,numpu)=help
        else
          uuh(i-i181,j,numpu)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.034).and.(isec1(7).eq.100)) then
      IF (TRIM(fpname) .EQ. 'VV') THEN
  ! V VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              !!!!!!  DJM - orig - if ((isec1(8)*100.0).eq.akz(ii)) numpv=ii
              if (current_grib_level .EQ. akz(ii)) numpv=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          vvh(i179+i,j,numpv)=help
        else
          vvh(i-i181,j,numpv)=help
        endif
      endif




!!!!!! DJM - orig -       if((isec1(6).eq.052).and.(isec1(7).eq.100)) then
      IF (TRIM(fpname) .EQ. 'RH') THEN
  ! RELATIVE HUMIDITY -> CONVERT TO SPECIFIC HUMIDITY LATER
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              !!!!!! DJM - orig - if ((isec1(8)*100.0).eq.akz(ii)) numprh=ii
              if (current_grib_level .EQ. akz(ii)) numprh=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          qvh(i179+i,j,numprh,n)=help
        else
          qvh(i-i181,j,numprh,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.039).and.(isec1(7).eq.100)) then
      IF (TRIM(fpname) .EQ. 'WW') THEN
  ! W VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              !!!!!! DJM - orig - if ((isec1(8)*100.0).eq.akz(ii)) numpw=ii
              if (current_grib_level .EQ. akz(ii)) numpw=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          wwh(i179+i,j,numpw)=help
        else
          wwh(i-i181,j,numpw)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.001).and.(isec1(7).eq.001)) then
      IF (TRIM(fpname) .EQ. 'PS') THEN
  ! SURFACE PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          ps(i179+i,j,1,n)=help
        else
          ps(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.066).and.(isec1(7).eq.001)) then
      IF (TRIM(fpname) .EQ. 'SD') THEN
  ! SNOW DEPTH
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          sd(i179+i,j,1,n)=help
        else
          sd(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.002).and.(isec1(7).eq.102)) then
      IF (TRIM(fpname) .EQ. 'SLP') THEN
  ! MEAN SEA LEVEL PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          msl(i179+i,j,1,n)=help
        else
          msl(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.071).and.(isec1(7).eq.244)) then
      IF (TRIM(fpname) .EQ. 'TCC') THEN
  ! TOTAL CLOUD COVER
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tcc(i179+i,j,1,n)=help
        else
          tcc(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.033).and.(isec1(7).eq.105).and. &
!!!!!!           (isec1(8).eq.10)) then
      IF (TRIM(fpname) .EQ. 'U10') THEN
  ! 10 M U VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          u10(i179+i,j,1,n)=help
        else
          u10(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.034).and.(isec1(7).eq.105).and. &
!!!!!!           (isec1(8).eq.10)) then
      IF (TRIM(fpname) .EQ. 'V10') THEN
  ! 10 M V VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          v10(i179+i,j,1,n)=help
        else
          v10(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.011).and.(isec1(7).eq.105).and. &
!!!!!!           (isec1(8).eq.02)) then
      IF (TRIM(fpname) .EQ. 'T2') THEN
  ! 2 M TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tt2(i179+i,j,1,n)=help
        else
          tt2(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.017).and.(isec1(7).eq.105).and. &
!!!!!!           (isec1(8).eq.02)) then
      IF (TRIM(fpname) .EQ. 'TD2') THEN
  ! 2 M DEW POINT TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          td2(i179+i,j,1,n)=help
        else
          td2(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.062).and.(isec1(7).eq.001)) then
     IF (TRIM(fpname) .EQ. 'LSPREC') THEN
  ! LARGE SCALE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          lsprec(i179+i,j,1,n)=help
        else
          lsprec(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.063).and.(isec1(7).eq.001)) then
     IF (TRIM(fpname) .EQ. 'CONVPREC') THEN
  ! CONVECTIVE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          convprec(i179+i,j,1,n)=help
        else
          convprec(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.007).and.(isec1(7).eq.001)) then
     IF (TRIM(fpname) .EQ. 'ORO') THEN
  ! TOPOGRAPHY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          oro(i179+i,j)=help
          excessoro(i179+i,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        else
          oro(i-i181,j)=help
          excessoro(i-i181,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.081).and.(isec1(7).eq.001)) then
     IF (TRIM(fpname) .EQ. 'LSM') THEN
  ! LAND SEA MASK
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          lsm(i179+i,j)=help
        else
          lsm(i-i181,j)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.221).and.(isec1(7).eq.001)) then
     IF (TRIM(fpname) .EQ. 'HMIX') THEN
  ! MIXING HEIGHT
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          hmix(i179+i,j,1,n)=help
        else
          hmix(i-i181,j,1,n)=help
        endif
      endif

!!!!!! DJM - orig -      if((isec1(6).eq.052).and.(isec1(7).eq.105).and. &
!!!!!!           (isec1(8).eq.02)) then
     IF (TRIM(fpname) .EQ. 'RH2') THEN
  ! 2 M RELATIVE HUMIDITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          qvh2(i179+i,j)=help
        else
          qvh2(i-i181,j)=help
        endif
      endif

!!!!!! DJM - orig -       if((isec1(6).eq.011).and.(isec1(7).eq.107)) then
     IF (TRIM(fpname) .EQ. 'TSIG1') THEN
  ! TEMPERATURE LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tlev1(i179+i,j)=help
        else
          tlev1(i-i181,j)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.033).and.(isec1(7).eq.107)) then
     IF (TRIM(fpname) .EQ. 'USIG1') THEN
  ! U VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          ulev1(i179+i,j)=help
        else
          ulev1(i-i181,j)=help
        endif
      endif

!!!!!! DJM - orig      if((isec1(6).eq.034).and.(isec1(7).eq.107)) then
     IF (TRIM(fpname) .EQ. 'VSIG1') THEN
  ! V VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          vlev1(i179+i,j)=help
        else
          vlev1(i-i181,j)=help
        endif
      endif

    end do
  end do

#endif
!!!!!!!  END OLD LOOP STRUCTURE - WILL EVENTUALLY BE DELETED




  endif !!!!! (TRIM(fpname) .NE. 'NOFP')

!!!!!! DJM - orig -   if((isec1(6).eq.33).and.(isec1(7).eq.100)) then
     IF (TRIM(fpname) .EQ. 'UU') THEN
  ! NCEP ISOBARIC LEVELS
    iumax=iumax+1
  endif

  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  !HSO close grib file
50   continue
  call grib_close_file(ifile)

  ! SENS. HEAT FLUX
  sshf(:,:,1,n)=0.0     ! not available from gfs.tccz.pgrbfxx files
  hflswitch=.false.    ! Heat flux not available
  ! SOLAR RADIATIVE FLUXES
  ssr(:,:,1,n)=0.0      ! not available from gfs.tccz.pgrbfxx files
  ! EW SURFACE STRESS
  ewss=0.0         ! not available from gfs.tccz.pgrbfxx files
  ! NS SURFACE STRESS
  nsss=0.0         ! not available from gfs.tccz.pgrbfxx files
  strswitch=.false.    ! stress not available

  ! CONVERT TP TO LSP (GRIB2 only)
  if (gribVer.eq.2) then
    do j=0,nymin1
    do i=0,nxfield-1
     if(i.le.i180) then
     if (convprec(i179+i,j,1,n).lt.lsprec(i179+i,j,1,n)) then ! neg precip would occur
         lsprec(i179+i,j,1,n)= &
              lsprec(i179+i,j,1,n)-convprec(i179+i,j,1,n)
     else
         lsprec(i179+i,j,1,n)=0
     endif
     else
     if (convprec(i-i181,j,1,n).lt.lsprec(i-i181,j,1,n)) then
          lsprec(i-i181,j,1,n)= &
               lsprec(i-i181,j,1,n)-convprec(i-i181,j,1,n)
     else
          lsprec(i-i181,j,1,n)=0
     endif
     endif
    enddo
    enddo
  endif
  !HSO end edits


  ! TRANSFORM RH TO SPECIFIC HUMIDITY

  do j=0,ny-1
    do i=0,nxfield-1
      do k=1,nuvz
        help=qvh(i,j,k,n)
        temp=tth(i,j,k,n)
        plev1=akm(k)+bkm(k)*ps(i,j,1,n)
        elev=ew(temp)*help/100.0
        qvh(i,j,k,n)=xmwml*(elev/(plev1-((1.0-xmwml)*elev)))
      end do
    end do
  end do

  ! CALCULATE 2 M DEW POINT FROM 2 M RELATIVE HUMIDITY
  ! USING BOLTON'S (1980) FORMULA
  ! BECAUSE td2 IS NOT AVAILABLE FROM NCEP GFS DATA

  do j=0,ny-1
    do i=0,nxfield-1
        help=qvh2(i,j)
        temp=tt2(i,j,1,n)
        elev=ew(temp)/100.*help/100.   !vapour pressure in hPa
        td2(i,j,1,n)=243.5/(17.67/log(elev/6.112)-1)+273.
        if (help.le.0.) td2(i,j,1,n)=tt2(i,j,1,n)
    end do
  end do

  if(levdiff2.eq.0) then
    iwmax=nlev_ec+1
    do i=0,nxmin1
      do j=0,nymin1
        wwh(i,j,nlev_ec+1)=0.
      end do
    end do
  endif


  ! For global fields, assign the leftmost data column also to the rightmost
  ! data column; if required, shift whole grid by nxshift grid points
  !*************************************************************************

  if (xglobal) then
    call shift_field_0(ewss,nxfield,ny)
    call shift_field_0(nsss,nxfield,ny)
    call shift_field_0(oro,nxfield,ny)
    call shift_field_0(excessoro,nxfield,ny)
    call shift_field_0(lsm,nxfield,ny)
    call shift_field_0(ulev1,nxfield,ny)
    call shift_field_0(vlev1,nxfield,ny)
    call shift_field_0(tlev1,nxfield,ny)
    call shift_field_0(qvh2,nxfield,ny)
    call shift_field(ps,nxfield,ny,1,1,2,n)
    call shift_field(sd,nxfield,ny,1,1,2,n)
    call shift_field(msl,nxfield,ny,1,1,2,n)
    call shift_field(tcc,nxfield,ny,1,1,2,n)
    call shift_field(u10,nxfield,ny,1,1,2,n)
    call shift_field(v10,nxfield,ny,1,1,2,n)
    call shift_field(tt2,nxfield,ny,1,1,2,n)
    call shift_field(td2,nxfield,ny,1,1,2,n)
    call shift_field(lsprec,nxfield,ny,1,1,2,n)
    call shift_field(convprec,nxfield,ny,1,1,2,n)
    call shift_field(sshf,nxfield,ny,1,1,2,n)
    call shift_field(ssr,nxfield,ny,1,1,2,n)
    call shift_field(hmix,nxfield,ny,1,1,2,n)
    call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
  endif

  do i=0,nxmin1
    do j=0,nymin1
  ! Convert precip. from mm/s -> mm/hour
      convprec(i,j,1,n)=convprec(i,j,1,n)*3600.
      lsprec(i,j,1,n)=lsprec(i,j,1,n)*3600.
      surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
  !  write(*,*) 'WARNING: No flux data contained in GRIB file ',
  !    +  wfname(indj)

  ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
  !***************************************************************************

    do i=0,nxmin1
      do j=0,nymin1
        hlev1=30.0                     ! HEIGHT OF FIRST MODEL SIGMA LAYER
        ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
        fflev1=sqrt(ulev1(i,j)**2+vlev1(i,j)**2)
        call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
             tt2(i,j,1,n),tlev1(i,j),ff10m,fflev1, &
             surfstr(i,j,1,n),sshf(i,j,1,n))
        if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
        if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
      end do
    end do
  endif

  if(iumax.ne.nuvz) stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iumax.ne.nwz)    stop 'READWIND: NWZ NOT CONSISTENT'

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  stop 'Execution terminated'
999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED !!!                    #### '
  stop 'Execution terminated'

end subroutine readwind_gfs
