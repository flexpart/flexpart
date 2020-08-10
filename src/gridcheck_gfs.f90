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

subroutine gridcheck_gfs

  !**********************************************************************
  !                                                                     *
  !             FLEXPART MODEL SUBROUTINE GRIDCHECK                     *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  !             AUTHOR:      G. WOTAWA                                  *
  !             DATE:        1997-08-06                                 *
  !             LAST UPDATE: 1997-10-10                                 *
  !                                                                     *
  !             Update:      1999-02-08, global fields allowed, A. Stohl*
  !             CHANGE: 17/11/2005, Caroline Forster, GFS data          *
  !             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
  !                                 ECMWF grib_api                      *
  !             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
  !                                 ECMWF grib_api                      *
  !                                                                     *
  !   Unified ECMWF and GFS builds                                      *
  !   Marian Harustak, 12.5.2017                                        *
  !     - Renamed routine from gridcheck to gridcheck_gfs               *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION:                                                        *
  !                                                                     *
  ! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
  ! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
  ! ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
  ! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
  ! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
  ! ANY CALL.                                                           *
  !                                                                     *
  ! XLON0                geographical longitude of lower left gridpoint *
  ! YLAT0                geographical latitude of lower left gridpoint  *
  ! NX                   number of grid points x-direction              *
  ! NY                   number of grid points y-direction              *
  ! DX                   grid distance x-direction                      *
  ! DY                   grid distance y-direction                      *
  ! NUVZ                 number of grid points for horizontal wind      *
  !                      components in z direction                      *
  ! NWZ                  number of grid points for vertical wind        *
  !                      component in z direction                       *
  ! sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
  !                      points of the polar stereographic grid):       *
  !                      used to check the CFL criterion                *
  ! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
  ! UVHEIGHT(NUVZ)       given                                          *
  ! WHEIGHT(1)-          heights of gridpoints where w is given         *
  ! WHEIGHT(NWZ)                                                        *
  !                                                                     *
  !**********************************************************************

  use eccodes
  use par_mod
  use com_mod
  use conv_mod
  use cmapf_mod, only: stlmbr,stcm2p

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  !HSO  end
  integer :: ix,jy,i,ifn,ifield,j,k,iumax,iwmax,numskip
  real :: sizesouth,sizenorth,xauxa,pint
  real :: akm_usort(nwzmax)
  real,parameter :: eps=0.0001

  ! NCEP GFS
  real :: pres(nwzmax), help

  integer :: i179,i180,i181

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  integer :: isec1(8),isec2(3)
  real(kind=4) :: zsec4(jpunp)
  character(len=1) :: opt

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheckwind_gfs'
  !
  if (numbnests.ge.1) then
  write(*,*) ' ###########################################'
  write(*,*) ' FLEXPART ERROR SUBROUTINE GRIDCHECK:'
  write(*,*) ' NO NESTED WINDFIELDAS ALLOWED FOR GFS!      '
  write(*,*) ' ###########################################'
  stop
  endif

  iumax=0
  iwmax=0

  if(ideltas.gt.0) then
    ifn=1
  else
    ifn=numbwf
  endif
  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
5   call grib_open_file(ifile,path(3)(1:length(3)) &
         //trim(wfname(ifn)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  call grib_multi_support_on

  ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE )  then
    goto 30    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !read the grib1 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'indicatorOfTypeOfLevel',isec1(7),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !get the size and data of the values array
  call grib_get_real4_array(igrib,'values',zsec4,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  else ! GRIB Edition 2

  !read the grib2 identifiers
  call grib_get_int(igrib,'discipline',discipl,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'scaledValueOfFirstFixedSurface', &
       valSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  if ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.100)) then ! U
    isec1(6)=33          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.3).and.(parNum.eq.5).and.(typSurf.eq.1)) then ! TOPO
    isec1(6)=7           ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.1) &
       .and.(discipl.eq.2)) then ! LSM
    isec1(6)=81          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  endif

  if (isec1(6).ne.-1) then
  !  get the size and data of the values array
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  endif ! gribVer

  if(ifield.eq.1) then

  !get the required fields from section 2
  !store compatible to gribex input
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
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

  ! Fix for flexpart.eu ticket #48
  if (xaux2in.lt.0) xaux2in = 359.0

  xaux1=xaux1in
  xaux2=xaux2in
  yaux1=yaux1in
  yaux2=yaux2in

  nxfield=isec2(2)
  ny=isec2(3)
  if((abs(xaux1).lt.eps).and.(xaux2.ge.359)) then ! NCEP DATA FROM 0 TO
    xaux1=-179.0                             ! 359 DEG EAST ->
    xaux2=-179.0+360.-360./real(nxfield)    ! TRANSFORMED TO -179
  endif                                      ! TO 180 DEG EAST
  if (xaux1.gt.180) xaux1=xaux1-360.0
  if (xaux2.gt.180) xaux2=xaux2-360.0
  if (xaux1.lt.-180) xaux1=xaux1+360.0
  if (xaux2.lt.-180) xaux2=xaux2+360.0
  if (xaux2.lt.xaux1) xaux2=xaux2+360.
  xlon0=xaux1
  ylat0=yaux1
  dx=(xaux2-xaux1)/real(nxfield-1)
  dy=(yaux2-yaux1)/real(ny-1)
  dxconst=180./(dx*r_earth*pi)
  dyconst=180./(dy*r_earth*pi)
  !HSO end edits


  ! Check whether fields are global
  ! If they contain the poles, specify polar stereographic map
  ! projections using the stlmbr- and stcm2p-calls
  !***********************************************************

    xauxa=abs(xaux2+dx-360.-xaux1)
    if (xauxa.lt.0.001) then
      nx=nxfield+1                 ! field is cyclic
      xglobal=.true.
      if (abs(nxshift).ge.nx) &
           stop 'nxshift in file par_mod is too large'
      xlon0=xlon0+real(nxshift)*dx
    else
      nx=nxfield
      xglobal=.false.
      if (nxshift.ne.0) &
           stop 'nxshift (par_mod) must be zero for non-global domain'
    endif
    nxmin1=nx-1
    nymin1=ny-1
    if (xlon0.gt.180.) xlon0=xlon0-360.
    xauxa=abs(yaux1+90.)
    if (xglobal.and.xauxa.lt.0.001) then
      sglobal=.true.               ! field contains south pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
      sizesouth=6.*(switchsouth+90.)/dy
      call stlmbr(southpolemap,-90.,0.)
      call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth, &
           sizesouth,switchsouth,180.)
      switchsouthg=(switchsouth-ylat0)/dy
    else
      sglobal=.false.
      switchsouthg=999999.
    endif
    xauxa=abs(yaux2-90.)
    if (xglobal.and.xauxa.lt.0.001) then
      nglobal=.true.               ! field contains north pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
      sizenorth=6.*(90.-switchnorth)/dy
      call stlmbr(northpolemap,90.,0.)
      call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth, &
           sizenorth,switchnorth,180.)
      switchnorthg=(switchnorth-ylat0)/dy
    else
      nglobal=.false.
      switchnorthg=999999.
    endif
  endif ! ifield.eq.1

  if (nxshift.lt.0) stop 'nxshift (par_mod) must not be negative'
  if (nxshift.ge.nxfield) stop 'nxshift (par_mod) too large'

  ! NCEP ISOBARIC LEVELS
  !*********************

  if((isec1(6).eq.33).and.(isec1(7).eq.100)) then ! check for U wind
    iumax=iumax+1
    pres(iumax)=real(isec1(8))*100.0
  endif


  i179=nint(179./dx)
  if (dx.lt.0.7) then
    i180=nint(180./dx)+1    ! 0.5 deg data
  else
    i180=nint(179./dx)+1    ! 1 deg data
  endif
  i181=i180+1


  ! NCEP TERRAIN
  !*************

  if((isec1(6).eq.007).and.(isec1(7).eq.001)) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        help=zsec4(nxfield*(ny-jy-1)+ix+1)
        if(ix.le.i180) then
          oro(i179+ix,jy)=help
          excessoro(i179+ix,jy)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        else
          oro(ix-i181,jy)=help
          excessoro(ix-i181,jy)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        endif
      end do
    end do
  endif

  ! NCEP LAND SEA MASK
  !*******************

  if((isec1(6).eq.081).and.(isec1(7).eq.001)) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        help=zsec4(nxfield*(ny-jy-1)+ix+1)
        if(ix.le.i180) then
          lsm(i179+ix,jy)=help
        else
          lsm(ix-i181,jy)=help
        endif
      end do
    end do
  endif

  call grib_release(igrib)

  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  ! HSO
30   continue
  call grib_close_file(ifile)
  ! HSO end edits

  nuvz=iumax
  nwz =iumax
  nlev_ec=iumax

  if (nx.gt.nxmax) then
   write(*,*) 'FLEXPART error: Too many grid points in x direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nx,nxmax
    stop
  endif

  if (ny.gt.nymax) then
   write(*,*) 'FLEXPART error: Too many grid points in y direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) ny,nymax
    stop
  endif

  if (nuvz.gt.nuvzmax) then
    write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nuvz,nuvzmax
    stop
  endif

  if (nwz.gt.nwzmax) then
    write(*,*) 'FLEXPART error: Too many w grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nwz,nwzmax
    stop
  endif

  ! If desired, shift all grids by nxshift grid cells
  !**************************************************

  if (xglobal) then
    call shift_field_0(oro,nxfield,ny)
    call shift_field_0(lsm,nxfield,ny)
    call shift_field_0(excessoro,nxfield,ny)
  endif

  ! Output of grid info
  !********************

  if (lroot) then
    write(*,*)
    write(*,*)
    write(*,'(a,2i7)') 'Vertical levels in NCEP data: ', &
         nuvz,nwz
    write(*,*)
    write(*,'(a)') 'Mother domain:'
    write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ', &
         xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
    write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range : ', &
         ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy
    write(*,*)
  end if

  ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
  ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

  numskip=nlev_ec-nuvz  ! number of ecmwf model layers not used
                        ! by trajectory model
  do i=1,nwz
    j=numskip+i
    k=nlev_ec+1+numskip+i
    akm_usort(nwz-i+1)=pres(nwz-i+1)
    bkm(nwz-i+1)=0.0
  end do

  !******************************
  ! change Sabine Eckhardt: akm should always be in descending order ... readwind adapted!
  !******************************
      do i=1,nwz
         if (akm_usort(1).gt.akm_usort(2)) then
            akm(i)=akm_usort(i)
         else
            akm(i)=akm_usort(nwz-i+1)
         endif
      end do

  !
  ! CALCULATION OF AKZ, BKZ
  ! AKZ,BKZ: model discretization parameters at the center of each model
  !     layer
  !
  ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
  ! i.e. ground level
  !*****************************************************************************

  do i=1,nuvz
     akz(i)=akm(i)
     bkz(i)=bkm(i)
  end do

  ! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
  ! upon the transformation to z levels. In order to save computer memory, this is
  ! not done anymore in the standard version. However, this option can still be
  ! switched on by replacing the following lines with those below, that are
  ! currently commented out. For this, similar changes are necessary in
  ! verttransform.f and verttranform_nests.f
  !*****************************************************************************

  nz=nuvz
  if (nz.gt.nzmax) stop 'nzmax too small'
  do i=1,nuvz
    aknew(i)=akz(i)
    bknew(i)=bkz(i)
  end do

  ! Switch on following lines to use doubled vertical resolution
  !*************************************************************
  !nz=nuvz+nwz-1
  !if (nz.gt.nzmax) stop 'nzmax too small'
  !do 100 i=1,nwz
  !  aknew(2*(i-1)+1)=akm(i)
  !00     bknew(2*(i-1)+1)=bkm(i)
  !do 110 i=2,nuvz
  !  aknew(2*(i-1))=akz(i)
  !10     bknew(2*(i-1))=bkz(i)
  ! End doubled vertical resolution


  ! Determine the uppermost level for which the convection scheme shall be applied
  ! by assuming that there is no convection above 50 hPa (for standard SLP)
  !*****************************************************************************

  do i=1,nuvz-2
    pint=akz(i)+bkz(i)*101325.
    if (pint.lt.5000.) goto 96
  end do
96   nconvlev=i
  if (nconvlev.gt.nconvlevmax-1) then
    nconvlev=nconvlevmax-1
    write(*,*) 'Attention, convection only calculated up to ', &
         akz(nconvlev)+bkz(nconvlev)*1013.25,' hPa'
  endif

  return

999   write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '       TRAJECTORY MODEL SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfname(ifn)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*)
  write(*,'(a)') '!!! PLEASE INSERT A NEW CD-ROM AND   !!!'
  write(*,'(a)') '!!! PRESS ANY KEY TO CONTINUE...     !!!'
  write(*,'(a)') '!!! ...OR TERMINATE FLEXPART PRESSING!!!'
  write(*,'(a)') '!!! THE "X" KEY...                   !!!'
  write(*,*)
  read(*,'(a)') opt
  if(opt.eq.'X') then
    stop
  else
    goto 5
  endif

end subroutine gridcheck_gfs
