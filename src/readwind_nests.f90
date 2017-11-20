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

subroutine readwind_nests(indj,n,uuhn,vvhn,wwhn)
  !                           i   i  o    o    o
  !*****************************************************************************
  !                                                                            *
  !     This routine reads the wind fields for the nested model domains.       *
  !     It is similar to subroutine readwind, which reads the mother domain.   *
  !                                                                            *
  !     Authors: A. Stohl, G. Wotawa                                           *
  !                                                                            *
  !     8 February 1999                                                        *
  !                                                                            *
  !     Last update: 17 October 2000, A. Stohl                                 *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !        Variables tthn and qvhn (on eta coordinates) in common block        *
  !  CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with ECMWF grib_api    *
  !  CHANGE: 03/12/2008, Harald Sodemann, update to f90 with ECMWF grib_api    *
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
  integer :: parId !!added by mc for making it consistent with new readwind.f90
  integer :: gotGrid
  !HSO  end

  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,l

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters


  !!!!!!!!!!!   DJM - eventually we will remove isec1()
  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec4(jpunp)
  real(kind=4) :: xaux,yaux
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real :: ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1
  real :: conversion_factor !added by mc to make it consistent with new gridchek.f90

  logical :: hflswitch,strswitch

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind_nests'

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




  do l=1,numbnests
    hflswitch=.false.
    strswitch=.false.
    levdiff2=nlev_ec-nwz+1
    iumax=0
    iwmax=0

    ifile=0
    igrib=0
    iret=0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!  Vtable code
  PRINT *, 'Loading Vtable: ', VTABLE_PATH
  call vtable_load_by_name(VTABLE_PATH, my_vtable)
  !! Debugging tool
  PRINT *, 'Dump of Vtable...'
  call vtable_dump_records(my_vtable)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!  VTABLE code
  ! This is diagnostic/debugging code, and will normally be commented out.
  ! It's purpose is to look at the provided grib file and produce an
  ! inventory of the FP-related messages, relative to the Vtable that's
  ! already been open.

  CALL vtable_gribfile_inventory( path(numpath+2*(l-1)+1)(1:length(numpath+2*(l-1)+1))// trim(wfnamen(l, indj)), &
&                                my_vtable)

  !!!!!!!!!!!!!!!!!!!  VTABLE code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !

5   call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,indj)),'r')
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages */
  !call grib_multi_support_on

    gotGrid=0
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
  print *, 'fpname: ', trim(fpname)


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
  if(ifield.eq.1) then
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfVerticalCoordinateValues', &
       isec2(12))
  call grib_check(iret,gribFunction,gribErrorMsg)
  ! CHECK GRID SPECIFICATIONS
  if(isec2(2).ne.nxn(l)) stop &
  'READWIND: NX NOT CONSISTENT FOR A NESTING LEVEL'
  if(isec2(3).ne.nyn(l)) stop &
  'READWIND: NY NOT CONSISTENT FOR A NESTING LEVEL'
  if(isec2(12)/2-1.ne.nlev_ec) stop 'READWIND: VERTICAL DISCRET&
       &IZATION NOT CONSISTENT FOR A NESTING LEVEL'
  endif ! ifield

  !HSO  get the second part of the grid dimensions only from GRiB1 messages
  if (TRIM(fpname) .EQ. 'T2'.and. (gotGrid.eq.0)) then
  !!!!! DJM --- if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then ! !added by mc to make it consisitent with new readwind.f90
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
         xauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    if (xauxin.gt.180.) xauxin=xauxin-360.0
    if (xauxin.lt.-180.) xauxin=xauxin+360.0

    xaux=xauxin
    yaux=yauxin
    if (abs(xaux-xlon0n(l)).gt.eps) &
    stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT FOR A NESTING LEVEL'
    if (abs(yaux-ylat0n(l)).gt.eps) &
    stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT FOR A NESTING LEVEL'
    gotGrid=1
  endif





!!!!!!!!!!!   BEGIN --  THIS IS OLD STUFF I'LL GET RID OF EVENTUALLY
#ifdef NO_LONGER_NEEDED
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        k=isec1(8)



        if(isec1(6).eq.130) tthn(i,j,nlev_ec-k+2,n,l)= &!! TEMPERATURE
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.131) uuhn(i,j,nlev_ec-k+2,l)= &!! U VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.132) vvhn(i,j,nlev_ec-k+2,l)= &!! V VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.133) then                         !! SPEC. HUMIDITY
          qvhn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (qvhn(i,j,nlev_ec-k+2,n,l) .lt. 0.) &
               qvhn(i,j,nlev_ec-k+2,n,l) = 0.
  !          this is necessary because the gridded data may contain
  !          spurious negative values
        endif
        if(isec1(6).eq.135) wwhn(i,j,nlev_ec-k+1,l)= &!! W VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)

        if(isec1(6).eq.134) psn(i,j,1,n,l)= &!! SURF. PRESS.
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)

        if(isec1(6).eq.141) sdn(i,j,1,n,l)= &!! SNOW DEPTH
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/conversion_factor !added by mc to make it consisitent with new readwind.f90!
        if(isec1(6).eq.151) msln(i,j,1,n,l)= &!! SEA LEVEL PRESS.
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.164) tccn(i,j,1,n,l)= &!! CLOUD COVER
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.165) u10n(i,j,1,n,l)= &!! 10 M U VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.166) v10n(i,j,1,n,l)= &!! 10 M V VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.167) tt2n(i,j,1,n,l)= &!! 2 M TEMPERATURE
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.168) td2n(i,j,1,n,l)= &!! 2 M DEW POINT
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)

        if(isec1(6).eq.142) then                         !! LARGE SCALE PREC.
          lsprecn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (lsprecn(i,j,1,n,l).lt.0.) lsprecn(i,j,1,n,l)=0.
        endif
        if(isec1(6).eq.143) then                         !! CONVECTIVE PREC.
          convprecn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/conversion_factor !added by mc to make it consisitent with new readwind.f90
          if (convprecn(i,j,1,n,l).lt.0.) convprecn(i,j,1,n,l)=0.
        endif

        if(isec1(6).eq.146) sshfn(i,j,1,n,l)= &!! SENS. HEAT FLUX
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if((isec1(6).eq.146).and. &
             (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.)) hflswitch=.true.    ! Heat flux available
        if(isec1(6).eq.176) then                         !! SOLAR RADIATION
          ssrn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (ssrn(i,j,1,n,l).lt.0.) ssrn(i,j,1,n,l)=0.
        endif
        if(isec1(6).eq.180) ewss(i,j)= &!! EW SURFACE STRESS
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.181) nsss(i,j)= &!! NS SURFACE STRESS
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(((isec1(6).eq.180).or.(isec1(6).eq.181)).and. &
             (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.)) strswitch=.true.    ! stress available


        if(isec1(6).eq.129) oron(i,j,l)= &!! ECMWF OROGRAPHY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
        if(isec1(6).eq.160) excessoron(i,j,l)= &!! STANDARD DEVIATION OF OROGRAPHY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.172) lsmn(i,j,l)= &!! ECMWF LAND SEA MASK
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)


! ESO TODO:
! -add check for if one of clwc/ciwc missing (error),
!    also if all 3 cw fields present, use qc and disregard the others
        if(isec1(6).eq.246) then  !! CLWC  Cloud liquid water content [kg/kg]
          clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          readclouds_nest(l)=.true.
          sumclouds_nest(l)=.false.
        endif
        if(isec1(6).eq.247) then  !! CIWC  Cloud ice water content
          ciwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        endif
!ZHG end
!ESO read qc (=clwc+ciwc)
        if(isec1(6).eq.201031) then  !! QC  Cloud liquid water content [kg/kg]
          clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          readclouds_nest(l)=.true.
          sumclouds_nest(l)=.true.
        endif


        if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
        if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

      end do
    end do


#endif
!!!!!!!!!!!   END --  THIS IS OLD STUFF I'LL GET RID OF EVENTUALLY




  k = current_grib_level

  IF(TRIM(fpname) .EQ. 'TT') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
          tthn(i,j,nlev_ec-k+2,n, l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'UU') THEN
      iumax=max(iumax,nlev_ec-k+1)
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
          uuhn(i,j,nlev_ec-k+2, l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'VV') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
          vvhn(i,j,nlev_ec-k+2, l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'QV') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
          qvhn(i,j,nlev_ec-k+2, n, l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      ! this is necessary because the gridded data may contain
      ! spurious negative values
      DO j=0,nyn(l)
        DO i=0,nxn(l)-1
            if (qvhn(i,j,nlev_ec-k+2, n, l) .lt. 0.) qvhn(i,j,nlev_ec-k+2, n, l) = 0.
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'ETADOT') THEN
      iwmax=max(iwmax,nlev_ec-k+1)
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
          wwhn(i,j,nlev_ec-k+1, l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'PS') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            psn(i,j,1,n,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
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

      conversion_factor = 1.0
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            sdn(i,j,1,n,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/conversion_factor
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'MSL') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
           msln(i,j,1,n,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'TCC') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            tccn(i,j,1,n,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'U10') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            u10n(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'V10') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            v10n(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'T2') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            tt2n(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'TD2') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            td2n(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'LSPREC') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            lsprecn(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            if (lsprecn(i,j,1,n,l).lt.0.) lsprecn(i,j,1,n,l)=0.
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
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            convprecn(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            if (convprecn(i,j,1,n,l).lt.0.) convprecn(i,j,1,n,l)=0.
        END DO
      END DO

  ELSE IF (TRIM(fpname) .EQ. 'ACPCP') THEN
      !!!!!!!! DJM - new code for GRIB2 convective precip.  Divide values by 1000
      !!       to convert from kg m-2 to m
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            convprecn(i,j,1,n,l) =  zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/1000.0
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            if (convprecn(i,j,1,n,l).lt.0.) convprecn(i,j,1,n,l)=0.
        END DO
      END DO



  ELSE IF(TRIM(fpname) .EQ. 'SHF') THEN
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            sshfn(i,j,1,n,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            IF (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.) hflswitch = .TRUE.  ! Heat flux available
        END DO
      END DO


  ELSE IF(TRIM(fpname) .EQ. 'SR') THEN                      !! SOLAR RADIATION
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            ssrn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            IF (ssrn(i,j,1,n,l).lt.0.) ssrn(i,j,1,n,l)=0.
        END DO
      END DO



!!!!  DJM - note that, unlike other variables, ewss and nsss are not stores in
!!!!  "nest-specific" arrays.  Rather, they are stored in the same arrays as the
!!!!  mother nest is stored.  This was in the original FPv10.1 code (and even FPv9.2,
!!!!  and likely earlier).  So, I left it this way when implementing Vtables.
  ELSE IF(TRIM(fpname) .EQ. 'EWSS') THEN   !! EW SURFACE STRESS
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            ewss(i,j) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            IF (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.) strswitch = .TRUE.  ! stress available
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'NSSS') THEN   !! NS SURFACE STRESS
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            nsss(i,j) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            IF (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.) strswitch = .TRUE.  ! stress available
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'ORO') THEN   !! ECMWF OROGRAPHY
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            !!!!!!! DJM - note - I don't know where "ga" comes from, but it was in the original code
            oron(i,j,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'EXCESSORO') THEN   !! STANDARD DEVIATION OF OROGRAPHY
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            excessoron(i,j,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'LSM') THEN   !! ECMWF LAND SEA MASK
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            lsmn(i,j,l) = zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
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
      readclouds_nest(l)=.true.
      sumclouds_nest(l)=.false.
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO

  ELSE IF(TRIM(fpname) .EQ. 'CICE') then  !! CIWC  Cloud ice water content
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            ciwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO
  !ZHG end
!ESO read qc (=clwc+ciwc)
  ELSE IF(TRIM(fpname) .EQ. 'QC') then  !! QC  Cloud liquid water content [kg/kg]
      readclouds_nest(l)=.true.
      sumclouds_nest(l)=.true.
      DO j=0,nyn(l)-1
        DO i=0,nxn(l)-1
            clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        END DO
      END DO


  END IF





  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
50   call grib_close_file(ifile)

  !error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  if(levdiff2.eq.0) then
    iwmax=nlev_ec+1
    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        wwhn(i,j,nlev_ec+1,l)=0.
      end do
    end do
  endif

  do i=0,nxn(l)-1
    do j=0,nyn(l)-1
      surfstrn(i,j,1,n,l)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
    write(*,*) 'WARNING: No flux data contained in GRIB file ', &
         wfnamen(l,indj)

  ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
  ! As ECMWF has increased the model resolution, such that now the first model
  ! level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
  ! (3rd model level in FLEXPART) for the profile method
  !***************************************************************************

    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        plev1=akz(3)+bkz(3)*psn(i,j,1,n,l)
        pmean=0.5*(psn(i,j,1,n,l)+plev1)
        tv=tthn(i,j,3,n,l)*(1.+0.61*qvhn(i,j,3,n,l))
        fu=-r_air*tv/ga/pmean
        hlev1=fu*(plev1-psn(i,j,1,n,l))   ! HEIGTH OF FIRST MODEL LAYER
        ff10m= sqrt(u10n(i,j,1,n,l)**2+v10n(i,j,1,n,l)**2)
        fflev1=sqrt(uuhn(i,j,3,l)**2+vvhn(i,j,3,l)**2)
        call pbl_profile(psn(i,j,1,n,l),td2n(i,j,1,n,l),hlev1, &
             tt2n(i,j,1,n,l),tthn(i,j,3,n,l),ff10m,fflev1, &
             surfstrn(i,j,1,n,l),sshfn(i,j,1,n,l))
        if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
        if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
      end do
    end do
  endif


  ! Assign 10 m wind to model level at eta=1.0 to have one additional model
  ! level at the ground
  ! Specific humidity is taken the same as at one level above
  ! Temperature is taken as 2 m temperature
  !**************************************************************************

    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        uuhn(i,j,1,l)=u10n(i,j,1,n,l)
        vvhn(i,j,1,l)=v10n(i,j,1,n,l)
        qvhn(i,j,1,n,l)=qvhn(i,j,2,n,l)
        tthn(i,j,1,n,l)=tt2n(i,j,1,n,l)
      end do
    end do

    if(iumax.ne.nuvz-1) stop &
         'READWIND: NUVZ NOT CONSISTENT FOR A NESTING LEVEL'
    if(iwmax.ne.nwz) stop &
         'READWIND: NWZ NOT CONSISTENT FOR A NESTING LEVEL'

  end do

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),' FOR NESTING LEVEL  #### '
  write(*,*) ' #### ',l,' IS NOT GRIB FORMAT !!!           #### '
  stop 'Execution terminated'


999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED FOR NESTING LEVEL ',l,'####'

end subroutine readwind_nests
