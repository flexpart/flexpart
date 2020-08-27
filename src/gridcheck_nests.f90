! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
  !*****************************************************************************

  use eccodes
  use par_mod
  use com_mod

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

  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec2(60+2*nuvzmax),zsec4(jpunp)

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheck_nests'

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

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !print*,'GRiB Edition 1'
  !read the grib2 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !change code for etadot to code for omega
  if (isec1(6).eq.77) then
    isec1(6)=135
  endif

  !print*,isec1(6),isec1(8)

  else

  !print*,'GRiB Edition 2'
  !read the grib2 identifiers
  call grib_get_int(igrib,'discipline',discipl,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',valSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'paramId',parId,iret) !added by mc to make it consisitent with new grid_check.f90
  call grib_check(iret,gribFunction,gribErrorMsg) !added by mc to make it consisitent with new  grid_check.f90

  !print*,discipl,parCat,parNum,typSurf,valSurf

  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  isec1(8)=valSurf     ! level
  if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! T
    isec1(6)=130         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.105)) then ! U
    isec1(6)=131         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.105)) then ! V
    isec1(6)=132         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! Q
    isec1(6)=133         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.83).and.(typSurf.eq.105)) then ! clwc
    isec1(6)=246         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.84).and.(typSurf.eq.105)) then ! ciwc
    isec1(6)=247         ! indicatorOfParameter
!ZHG end
! ESO qc(=clwc+ciwc)
  elseif ((parCat.eq.201).and.(parNum.eq.31).and.(typSurf.eq.105)) then ! qc
    isec1(6)=201031      ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.1)) then !SP
    isec1(6)=134         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.32)) then ! W, actually eta dot
    isec1(6)=135         ! indicatorOfParameter
  elseif ((parCat.eq.128).and.(parNum.eq.77)) then ! W, actually eta dot !added bymc to make it consistent with new gridcheck.f90
    isec1(6)=135         ! indicatorOfParameter    !
  elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.101)) then !SLP
    isec1(6)=151         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.103)) then ! 10U
    isec1(6)=165         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.103)) then ! 10V
    isec1(6)=166         ! indicatorOfParameter
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.103)) then ! 2T
    isec1(6)=167         ! indicatorOfParameter
  elseif ((parCat.eq.0).and.(parNum.eq.6).and.(typSurf.eq.103)) then ! 2D
    isec1(6)=168         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SD
    isec1(6)=141         ! indicatorOfParameter
  elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC !added by mc to make it consistent with new gridchek.f90
    isec1(6)=164         ! indicatorOfParameter
 elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP !added by mc to make it consistent with new gridchek.f90
    isec1(6)=142         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
    isec1(6)=143         ! indicatorOfParameter
  elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
    isec1(6)=146         ! indicatorOfParameter
  elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
    isec1(6)=176         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.17) .or. parId .eq. 180) then ! EWSS !added by mc to make it consistent with new gridchek.f90
    isec1(6)=180         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.18) .or. parId .eq. 181) then ! NSSS !added by mc to make it consistent with new gridchek.f90
    isec1(6)=181         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
    isec1(6)=129         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO !added by mc to make it consistent with new gridchek.f90
    isec1(6)=160         ! indicatorOfParameter
  elseif ((discipl.eq.2).and.(parCat.eq.0).and.(parNum.eq.0).and. &
       (typSurf.eq.1)) then ! LSM
    isec1(6)=172         ! indicatorOfParameter
  else
    print*,'***ERROR: undefined GRiB2 message found!',discipl, &
         parCat,parNum,typSurf
  endif
  if(parId .ne. isec1(6) .and. parId .ne. 77) then !added by mc to make it consistent with new gridchek.f90
    write(*,*) 'parId',parId, 'isec1(6)',isec1(6)
!    stop
  endif

  endif

  !get the size and data of the values array
  if (isec1(6).ne.-1) then
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
  if (isec1(6) .eq. 167 .and. (gotGrib.eq.0)) then !added by mc to make it consistent with new gridchek.f90 note that gotGrid must be changed in gotGrib!!
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

  k=isec1(8)
  if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
  if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

  if(isec1(6).eq.129) then
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        oron(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif
  if(isec1(6).eq.172) then
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        lsmn(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif
  if(isec1(6).eq.160) then
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
