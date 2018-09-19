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

subroutine readwind_ecmwf(indj,n,uuh,vvh,wwh)

!**********************************************************************
!                                                                     *
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     *
!             AUTHOR:      G. WOTAWA                                  *
!             DATE:        1997-08-05                                 *
!             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
!             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
!                                 ECMWF grib_api                      *
!             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
!                                 ECMWF grib_api                      *
!                                                                     *
!**********************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:
!   Variables tth and qvh (on eta coordinates) in common block
!
!   Unified ECMWF and GFS builds                                      
!   Marian Harustak, 12.5.2017                                        
!     - Renamed from readwind to readwind_ecmwf                     
!
!                                                                     *
!   Implementation of the Vtables approach                            *
!   D. Morton, D. Arnold 17.11.2017                                   *
!     - Inclusion of specific code and usage of class_vtable          *
!                                                                     *  !                                                                     *
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
!                                                                     *
! wfname             File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of vertical levels ecmwf model            *
! uu,vv,ww           wind fields                                      *
! tt,qv              temperature and specific humidity                *
! ps                 surface pressure                                 *
!                                                                     *
!**********************************************************************

  use grib_api
  use par_mod
  use com_mod
  use class_vtable

  implicit none

!  include 'grib_api.h'

!HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl,parId
  integer :: gotGrid
!HSO  end

  real(sp) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real(sp) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real(sp) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax

! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

! dimension of isec2 at least (22+n), where n is the number of parallels or
! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

! dimension of zsec2 at least (10+nn), where nn is the number of vertical
! coordinate parameters

  !!!!!!!! DJM - orig - integer :: isec1(56),isec2(22+nxmax+nymax)
  !!!!!!!! DJM - with Vtables, isec1() array no longer necessary
  integer :: isec2(22+nxmax+nymax)
  real(sp) :: zsec4(jpunp)
  real(sp) :: xaux,yaux
  real(dp) :: xauxin,yauxin
  real(sp),parameter :: eps=1.e-4
  real(sp) :: nsss(0:nxmax-1,0:nymax-1),ewss(0:nxmax-1,0:nymax-1)
  real(sp) :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1,conversion_factor

  logical :: hflswitch,strswitch

!HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Vtable related variables
  !
  !  Path to Vtable - current implementation assumes it's in cwd, named
  !  "Vtable"
  ! ESO: Changed to use default Vtable file in options directory
  ! CHARACTER(LEN=255), PARAMETER :: VTABLE_PATH = "Vtable"
  character(LEN=255) :: VTABLE_PATH
  character(LEN=15) :: fpname      ! stores FLEXPART name for curr grib mesg.
  type(Vtable),save :: my_vtable    ! unallocated
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!  DJM
  integer current_grib_level   ! this was isec1(8) in previous versions
  logical :: linit=.true.




  hflswitch=.false.
  strswitch=.false.
!ZHG test the grib fields that have lcwc without using them
!  readcloud=.false.

  levdiff2=nlev_ec-nwz+1
  iumax=0
  iwmax=0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Vtable code
  if (linit) then
     VTABLE_PATH = path(1)(1:length(1))//'Vtables/Vtable.ecmwf'
     PRINT *, 'Loading Vtable: ', VTABLE_PATH
     call vtable_load_by_name(VTABLE_PATH, my_vtable)
     linit=.false.
  end if
  !! Debugging tool
!!  PRINT *, 'Dump of Vtable...'
  ! call vtable_dump_records(my_vtable)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  ! This is diagnostic/debugging code, and will normally be commented out.
  ! It's purpose is to look at the provided grib file and produce an
  ! inventory of the FP-related messages, relative to the Vtable that's
  ! already been open.

!!  CALL vtable_gribfile_inventory(path(3)(1:length(3)) // trim(wfname(indj)), &
!! &                                my_vtable)

  !!!!!!!!!!!!!!!!!!!  VTABLE code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!
! OPENING OF DATA FILE (GRIB CODE)
!
5 call grib_open_file(ifile,path(3)(1:length(3)) &
       //trim(wfname(indj)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif
!turn on support for multi fields messages */
!call grib_multi_support_on

  gotGrid=0
  ifield=0
10 ifield=ifield+1
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
  call grib_check(iret,gribFunction,gribErrorMsg)


  !!  DJM - get the current_grib_level (in previous code it was isec1(8))
  !!  It's the same in both GRIB1 and GRIB2
  call grib_get_int(igrib,'level',current_grib_level,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)




!HSO  get the size and data of the values array
  !!!! -- original statement -- if (isec1(6).ne.-1) then
  ! 'NOFP' is the fpname for a grib message not recognized by FLEXPART
  IF (TRIM(fpname) .NE. 'NOFP') THEN
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  ENDIF

!HSO  get the required fields from section 2 in a gribex compatible manner
  if (ifield.eq.1) then
    call grib_get_int(igrib,'numberOfPointsAlongAParallel',isec2(2),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfPointsAlongAMeridian',isec2(3),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfVerticalCoordinateValues',isec2(12))
    call grib_check(iret,gribFunction,gribErrorMsg)
! CHECK GRID SPECIFICATIONS
    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    if(isec2(12)/2-1.ne.nlev_ec) &
         stop 'READWIND: VERTICAL DISCRETIZATION NOT CONSISTENT'
  endif ! ifield

!HSO  get the second part of the grid dimensions only from GRiB1 messages


  if (TRIM(fpname) .EQ. 'T2'.and. (gotGrid.eq.0)) then
  !!!!! DJM --- if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
         xauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    if (xauxin.gt.180.) xauxin=xauxin-360.0
    if (xauxin.lt.-180.) xauxin=xauxin+360.0

    xaux=xauxin+real(nxshift)*dx
    yaux=yauxin
    if (xaux.gt.180.) xaux=xaux-360.0
    if(abs(xaux-xlon0).gt.eps) &
         stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
    if(abs(yaux-ylat0).gt.eps) &
         stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
    gotGrid=1
  endif ! gotGrid


  k = current_grib_level

  IF(TRIM(fpname) .EQ. 'TT') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
          tth(i,j,nlev_ec-k+2,n) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO


  ELSE IF(TRIM(fpname) .EQ. 'UU') THEN
      iumax=max(iumax,nlev_ec-k+1)
      DO j=0,nymin1
        DO i=0,nxfield-1
          uuh(i,j,nlev_ec-k+2) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'VV') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
          vvh(i,j,nlev_ec-k+2) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'QV') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            qvh(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

      ! this is necessary because the gridded data may contain
      ! spurious negative values
      DO j=0,nymin1
        DO i=0,nxfield-1
            if (qvh(i,j,nlev_ec-k+2,n) .lt. 0.) qvh(i,j,nlev_ec-k+2,n) = 0.
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'ETADOT') THEN
      iwmax=max(iwmax,nlev_ec-k+1)
      DO j=0,nymin1
        DO i=0,nxfield-1
          wwh(i,j,nlev_ec-k+1) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'PS') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            ps(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'SD') THEN
      !!!!!!!! DJM - WARNING - in the previous version of this code, snow depth
      !!       had been assumed to be in m if it was GRIB1, and mm if it was GRIB2.
      !!       Hence, if the values were from a GRIB2 message, they were divided by
      !!       1000 to convert from mm to m.
      !!       This was done by Leo, based on his experience with the GRIB files.  My
      !!       experience has so far been that units in both GRIB1 and GRIB2 messages
      !!       are in m, which means no conversion would be necessary.
      !!
      !!       I have therefore set the value of "conversion_factor" to 1.0 for reading
      !!       in snow depth, but I don't feel 100% good about this just yet.  It may
      !!       need to be scrutinized more closely in the future.

! ESO: reverted conversion factor to 1000.
    conversion_factor = 1000.0
      DO j=0,nymin1
        DO i=0,nxfield-1
            sd(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'MSL') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
           msl(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'TCC') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            tcc(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'U10') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            u10(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'V10') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            v10(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'T2') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            tt2(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'TD2') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            td2(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'LSPREC') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            lsprec(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            if (lsprec(i,j,1,n).lt.0.) lsprec(i,j,1,n)=0.
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'CONVPREC') THEN
      !!!!!!!! DJM - In the previous version of this code, if convective precip was
      !!       read from a GRIB1 message, it was assumed to have units of "m."  If it
      !!       was read in from a GRIB2 message, it was assumed to have units of
      !!       kg m-2, and then a "conversion_factor" of 1000 was assigned, and incoming
      !!       values would be divided by that factor.  This had been implemented by Leo.
      !!
      !!       In GRIB1 messages, the shortname is "cp" and in GRIB2 it's "acpcp."  So,
      !!       my implementation has an fpname of CONVPREC if it's GRIB1, and ACPCP if
      !!       it's GRIB2.  And, ACPCP data is divided by 1000 as it's read in.
      !!
      !!
      DO j=0,nymin1
        DO i=0,nxfield-1
! ESO: Awaiting CTBTO clarification, conversion factor=1000 added for now
            convprec(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)/1000.
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'ACPCP') THEN
      !!!!!!!! DJM - new code for GRIB2 convective precip.  Divide values by 1000
      !!       to convert from kg m-2 to m
      DO j=0,nymin1
        DO i=0,nxfield-1
            convprec(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)/1000.0
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'SHF') THEN
      DO j=0,nymin1
        DO i=0,nxfield-1
            sshf(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            IF (zsec4(nxfield*(ny-j-1)+i+1).ne.0.) hflswitch = .TRUE.  ! Heat flux available
        END DO
      END DO


  ELSE IF(TRIM(fpname) .EQ. 'SR') THEN                      !! SOLAR RADIATION
      DO j=0,nymin1
        DO i=0,nxfield-1
            ssr(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            IF (ssr(i,j,1,n).lt.0.) ssr(i,j,1,n)=0.
        END DO
      END DO



  ELSE IF(TRIM(fpname) .EQ. 'EWSS') THEN   !! EW SURFACE STRESS
      DO j=0,nymin1
        DO i=0,nxfield-1
            ewss(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            IF (zsec4(nxfield*(ny-j-1)+i+1).ne.0.) strswitch = .TRUE.  ! stress available
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'NSSS') THEN   !! NS SURFACE STRESS
      DO j=0,nymin1
        DO i=0,nxfield-1
            nsss(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
      DO j=0,nymin1
        DO i=0,nxfield-1
            IF (zsec4(nxfield*(ny-j-1)+i+1).ne.0.) strswitch = .TRUE.  ! stress available
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'ORO') THEN   !! ECMWF OROGRAPHY
      DO j=0,nymin1
        DO i=0,nxfield-1
            !!!!!!! DJM - note - I don't know where "ga" comes from, but it was in the original code
            oro(i,j) = zsec4(nxfield*(ny-j-1)+i+1)/ga
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'EXCESSORO') THEN   !! STANDARD DEVIATION OF OROGRAPHY
      DO j=0,nymin1
        DO i=0,nxfield-1
            excessoro(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'LSM') THEN   !! ECMWF LAND SEA MASK
      DO j=0,nymin1
        DO i=0,nxfield-1
            lsm(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO


  !!!!!!!! DJM - I'm under the impression that the reading of CLWC, CIWC and QC is still under
  !!       development.  But, I've gone ahead and replaced the grib code with an fpname.  There are
  !!       additional notes in the comments section of options/Vtables/Vtable.ecmwf concerning the
  !!       grib codes for these three messages

  !ZHG READING CLOUD FIELDS ASWELL
  ! ESO TODO: add check for if one of clwc/ciwc missing (error),
  ! also if all 3 cw fields present, use qc and disregard the others
  ELSE IF(TRIM(fpname) .EQ. 'CLWC') THEN  !! CLWC  Cloud liquid water content [kg/kg]
      readclouds=.true.
      sumclouds=.false.
      DO j=0,nymin1
        DO i=0,nxfield-1
            clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'CIWC') then  !! CIWC  Cloud ice water content
      DO j=0,nymin1
        DO i=0,nxfield-1
            ciwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO
  !ZHG end

  ELSE IF(TRIM(fpname) .EQ. 'QC') then  !! QC  Cloud liquid water content [kg/kg]
      readclouds=.true.
      sumclouds=.true.
      DO j=0,nymin1
        DO i=0,nxfield-1
            clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        END DO
      END DO


  END IF










#ifdef NO_LONGER_NEEDED
  do j=0,nymin1
    do i=0,nxfield-1

      !!!!!!!! DJM - orig - k=isec1(8)
      k = current_grib_level


      !!!!!!!! DJM - orig - if(isec1(6).eq.130) tth(i,j,nlev_ec-k+2,n)= &!! TEMPERATURE
      if(TRIM(fpname) .EQ. 'TT') tth(i,j,nlev_ec-k+2,n)= &!! TEMPERATURE
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.131) uuh(i,j,nlev_ec-k+2)= &!! U VELOCITY
      if(TRIM(fpname) .EQ. 'UU') uuh(i,j,nlev_ec-k+2)= &!! U VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.132) vvh(i,j,nlev_ec-k+2)= &!! V VELOCITY
      if(TRIM(fpname) .EQ. 'VV') vvh(i,j,nlev_ec-k+2)= &!! V VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.133) then                      !! SPEC. HUMIDITY
      if(TRIM(fpname) .EQ. 'QV') then                      !! SPEC. HUMIDITY
        qvh(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (qvh(i,j,nlev_ec-k+2,n) .lt. 0.) &
             qvh(i,j,nlev_ec-k+2,n) = 0.
!        this is necessary because the gridded data may contain
!        spurious negative values
      endif
      !!!!!!!! DJM - orig - if(isec1(6).eq.135) wwh(i,j,nlev_ec-k+1)= &!! W VELOCITY
      if(TRIM(fpname) .EQ. 'ETADOT') wwh(i,j,nlev_ec-k+1)= &!! W VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)

      !!!!!!!! DJM - orig - if(isec1(6).eq.134) ps(i,j,1,n)= &!! SURF. PRESS.
      if(TRIM(fpname) .EQ. 'PS') ps(i,j,1,n)= &!! SURF. PRESS.
           zsec4(nxfield*(ny-j-1)+i+1)


      !!!!!!!! DJM - WARNING - in the previous version of this code, snow depth
      !!       had been assumed to be in m if it was GRIB1, and mm if it was GRIB2.
      !!       Hence, if the values were from a GRIB2 message, they were divided by
      !!       1000 to convert from mm to m.
      !!       This was done by Leo, based on his experience with the GRIB files.  My
      !!       experience has so far been that units in both GRIB1 and GRIB2 messages
      !!       are in m, which means no conversion would be necessary.
      !!
      !!       I have therefore set the value of "conversion_factor" to 1.0 for reading
      !!       in snow depth, but I don't feel 100% good about this just yet.  It may
      !!       need to be scrutinized more closely in the future.
      conversion_factor = 1.0
      !!!!!!!! DJM - orig - if(isec1(6).eq.141) sd(i,j,1,n)= &!! SNOW DEPTH
      if(TRIM(fpname) .EQ. 'SD') sd(i,j,1,n)= &!! SNOW DEPTH
           zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor

      !!!!!!!! DJM - orig - if(isec1(6).eq.151) msl(i,j,1,n)= &!! SEA LEVEL PRESS.
      if(TRIM(fpname) .EQ. 'MSL') msl(i,j,1,n)= &!! SEA LEVEL PRESS.
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.164) tcc(i,j,1,n)= &!! CLOUD COVER
      if(TRIM(fpname) .EQ. 'TCC') tcc(i,j,1,n)= &!! CLOUD COVER
           zsec4(nxfield*(ny-j-1)+i+1)

      !!!!!!!! DJM - orig - if(isec1(6).eq.165) u10(i,j,1,n)= &!! 10 M U VELOCITY
      if(TRIM(fpname) .EQ. 'U10') u10(i,j,1,n)= &!! 10 M U VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.166) v10(i,j,1,n)= &!! 10 M V VELOCITY
      if(TRIM(fpname) .EQ. 'V10') v10(i,j,1,n)= &!! 10 M V VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.167) tt2(i,j,1,n)= &!! 2 M TEMPERATURE
      if(TRIM(fpname) .EQ. 'T2') tt2(i,j,1,n)= &!! 2 M TEMPERATURE
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.168) td2(i,j,1,n)= &!! 2 M DEW POINT
      if(TRIM(fpname) .EQ. 'TD2') td2(i,j,1,n)= &!! 2 M DEW POINT
           zsec4(nxfield*(ny-j-1)+i+1)

      !!!!!!!! DJM -orig - if(isec1(6).eq.142) then                      !! LARGE SCALE PREC.
      if(TRIM(fpname) .EQ. 'LSPREC') then                      !! LARGE SCALE PREC.
        lsprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (lsprec(i,j,1,n).lt.0.) lsprec(i,j,1,n)=0.
      endif

      !!!!!!!! DJM - In the previous version of this code, if convective precip was
      !!       read from a GRIB1 message, it was assumed to have units of "m."  If it
      !!       was read in from a GRIB2 message, it was assumed to have units of
      !!       kg m-2, and then a "conversion_factor" of 1000 was assigned, and incoming
      !!       values would be divided by that factor.  This had been implemented by Leo.
      !!
      !!       In GRIB1 messages, the shortname is "cp" and in GRIB2 it's "acpcp."  So,
      !!       my implementation has an fpname of CONVPREC if it's GRIB1, and ACPCP if
      !!       it's GRIB2.  And, ACPCP data is divided by 1000 as it's read in.
      !!
      !!

      !!!!!!!! DJM -orig - if(isec1(6).eq.143) then                      !! CONVECTIVE PREC.
      if(TRIM(fpname) .EQ. 'CONVPREC') then                      !! CONVECTIVE PREC.
        !!!!!  DJM - orig - convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
        convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
      endif

      !!!!!!!! DJM - new code for GRIB2 convective precip.  Divide values by 1000
      !!       to convert from kg m-2 to m
      if(TRIM(fpname) .EQ. 'ACPCP') then                      !! CONVECTIVE PREC.
        convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)/1000.0
        if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
      endif

      !!!!!!!! DJM - orig - if(isec1(6).eq.146) sshf(i,j,1,n)= &!! SENS. HEAT FLUX
      if(TRIM(fpname) .EQ. 'SHF') sshf(i,j,1,n)= &!! SENS. HEAT FLUX
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if((isec1(6).eq.146).and.(zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) &
      if((TRIM(fpname) .EQ. 'SHF').and.(zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) &
           hflswitch=.true.    ! Heat flux available
      !!!!!!!! DJM - orig - if(isec1(6).eq.176) then                      !! SOLAR RADIATION
      if(TRIM(fpname) .EQ. 'SR') then                      !! SOLAR RADIATION
        ssr(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (ssr(i,j,1,n).lt.0.) ssr(i,j,1,n)=0.
      endif

      !!!!!!!! DJM - orig - if(isec1(6).eq.180) ewss(i,j)= &!! EW SURFACE STRESS
      if(TRIM(fpname) .EQ. 'EWSS') ewss(i,j)= &!! EW SURFACE STRESS
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.181) nsss(i,j)= &!! NS SURFACE STRESS
      if(TRIM(fpname) .EQ. 'NSSS') nsss(i,j)= &!! NS SURFACE STRESS
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(((isec1(6).eq.180).or.(isec1(6).eq.181)).and. &
      if(((TRIM(fpname) .EQ. 'EWSS').or.(TRIM(fpname) .EQ. 'NSSS')).and. &
           (zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) strswitch=.true.    ! stress available
!sec        strswitch=.true.

      !!!!!!!! DJM - orig - if(isec1(6).eq.129) oro(i,j)= &!! ECMWF OROGRAPHY
      if(TRIM(fpname) .EQ. 'ORO') oro(i,j)= &!! ECMWF OROGRAPHY
           zsec4(nxfield*(ny-j-1)+i+1)/ga
      !!!!!!!! DJM - orig - if(isec1(6).eq.160) excessoro(i,j)= &!! STANDARD DEVIATION OF OROGRAPHY
      if(TRIM(fpname) .EQ. 'EXCESSORO') excessoro(i,j)= &!! STANDARD DEVIATION OF OROGRAPHY
           zsec4(nxfield*(ny-j-1)+i+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.172) lsm(i,j)= &!! ECMWF LAND SEA MASK
      if(TRIM(fpname) .EQ. 'LSM') lsm(i,j)= &!! ECMWF LAND SEA MASK
           zsec4(nxfield*(ny-j-1)+i+1)

      !!!!!!!! DJM - orig - if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
      if(TRIM(fpname) .EQ. 'UU') iumax=max(iumax,nlev_ec-k+1)
      !!!!!!!! DJM - orig - if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)
      if(TRIM(fpname) .EQ. 'ETADOT') iwmax=max(iwmax,nlev_ec-k+1)


!!!!!!!! DJM - I'm under the impression that the reading of CLWC, CIWC and QC is still under
!!       development.  But, I've gone ahead and replaced the grib code with an fpname.  There are
!!       additional notes in the comments section of options/Vtables/Vtable.ecmwf concerning the
!!       grib codes for these three messages

!ZHG READING CLOUD FIELDS ASWELL
! ESO TODO: add check for if one of clwc/ciwc missing (error),
! also if all 3 cw fields present, use qc and disregard the others
      !!!!!!!! DJM - orig - if(isec1(6).eq.246) then  !! CLWC  Cloud liquid water content [kg/kg]
      if(TRIM(fpname) .EQ. 'CLWC') then  !! CLWC  Cloud liquid water content [kg/kg]
        clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        readclouds=.true.
        sumclouds=.false.
      endif
      !!!!!!!! DJM - orig - if(isec1(6).eq.247) then  !! CIWC  Cloud ice water content
      if(TRIM(fpname) .EQ. 'CICE') then  !! CIWC  Cloud ice water content
        ciwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
      endif
!ZHG end
!ESO read qc (=clwc+ciwc)
      !!!!!!!! DJM - orig - if(isec1(6).eq.201031) then  !! QC  Cloud liquid water content [kg/kg]
      if(TRIM(fpname) .EQ. 'QC') then  !! QC  Cloud liquid water content [kg/kg]
        clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        readclouds=.true.
        sumclouds=.true.
      endif


    end do
  end do

#endif


  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
!
! CLOSING OF INPUT DATA FILE
!

50 call grib_close_file(ifile)

!error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

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
    call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
!ZHG
    call shift_field(clwch,nxfield,ny,nuvzmax,nuvz,2,n)
    if (.not.sumclouds) call shift_field(ciwch,nxfield,ny,nuvzmax,nuvz,2,n)
!ZHG end

  endif

  do i=0,nxmin1
    do j=0,nymin1
      surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
    write(*,*) 'WARNING: No flux data contained in GRIB file ', &
         wfname(indj)

! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
! As ECMWF has increased the model resolution, such that now the first model
! level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
! (3rd model level in FLEXPART) for the profile method
!***************************************************************************

    do i=0,nxmin1
      do j=0,nymin1
        plev1=akz(3)+bkz(3)*ps(i,j,1,n)
        pmean=0.5*(ps(i,j,1,n)+plev1)
        tv=tth(i,j,3,n)*(1.+0.61*qvh(i,j,3,n))
        fu=-r_air*tv/ga/pmean
        hlev1=fu*(plev1-ps(i,j,1,n))   ! HEIGTH OF FIRST MODEL LAYER
        ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
        fflev1=sqrt(uuh(i,j,3)**2+vvh(i,j,3)**2)
        call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
             tt2(i,j,1,n),tth(i,j,3,n),ff10m,fflev1, &
             surfstr(i,j,1,n),sshf(i,j,1,n))
        if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
        if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
      end do
    end do
  endif


! Assign 10 m wind to model level at eta=1.0 to have one additional model
! level at the ground
! Specific humidity is taken the same as at one level above
! Temperature is taken as 2 m temperature
!**************************************************************************

  do i=0,nxmin1
    do j=0,nymin1
      uuh(i,j,1)=u10(i,j,1,n)
      vvh(i,j,1)=v10(i,j,1,n)
      qvh(i,j,1,n)=qvh(i,j,2,n)
      tth(i,j,1,n)=tt2(i,j,1,n)
    end do
  end do

  if(iumax.ne.nuvz-1) stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iwmax.ne.nwz)    stop 'READWIND: NWZ NOT CONSISTENT'

  return

888 write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  stop 'Execution terminated'

end subroutine readwind_ecmwf

