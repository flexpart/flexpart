! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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

  do l=1,numbnests
    hflswitch=.false.
    strswitch=.false.
    levdiff2=nlev_ec-nwz+1
    iumax=0
    iwmax=0

    ifile=0
    igrib=0
    iret=0

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

  conversion_factor=1.


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
  call grib_get_int(igrib,'paramId',parId,iret) !added by mc to make it consisitent with new readwind.f90
  call grib_check(iret,gribFunction,gribErrorMsg) !added by mc to make it consisitent with new readwind.f90

  !print*,discipl,parCat,parNum,typSurf,valSurf

  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  isec1(8)=valSurf     ! level
   conversion_factor=1.
  if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! T
    isec1(6)=130         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.105)) then ! U
    isec1(6)=131         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.105)) then ! V
    isec1(6)=132         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! Q
    isec1(6)=133         ! indicatorOfParameter
! ESO Cloud water is in a) fields CLWC and CIWC, *or* b) field QC 
    elseif ((parCat.eq.1).and.(parNum.eq.83).and.(typSurf.eq.105)) then ! clwc
      isec1(6)=246         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.84).and.(typSurf.eq.105)) then ! ciwc
      isec1(6)=247         ! indicatorOfParameter
! ESO qc(=clwc+ciwc):
    elseif ((parCat.eq.201).and.(parNum.eq.31).and.(typSurf.eq.105)) then ! qc
      isec1(6)=201031         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.1)) then !SP
    isec1(6)=134         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.32)) then ! W, actually eta dot !
    isec1(6)=135         ! indicatorOfParameter
  elseif ((parCat.eq.128).and.(parNum.eq.77)) then ! W, actually eta dot !added by mc to make it consisitent with new readwind.f90
    isec1(6)=135         ! indicatorOfParameter    !added by mc to make it consisitent with new readwind.f90
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
    conversion_factor=1000. !added by mc to make it consisitent with new readwind.f90
  elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC !added by mc to make it consisitent with new readwind.f90
    isec1(6)=164         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP !added by mc to make it consisitent with new readwind.f90
    isec1(6)=142         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
    isec1(6)=143         ! indicatorOfParameter
    conversion_factor=1000. !added by mc to make it consisitent with new readwind.f90
  elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
    isec1(6)=146         ! indicatorOfParameter
  elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
    isec1(6)=176         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.38) .or. parId .eq. 180) then ! EWSS !added by mc to make it consisitent with new readwind.f90
    isec1(6)=180         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.37) .or. parId .eq. 181) then ! NSSS !added by mc to make it consisitent with new readwind.f90
    isec1(6)=181         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
    isec1(6)=129         ! indicatorOfParameter
   elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO !added by mc to make it consisitent with new readwind.f90
    isec1(6)=160         ! indicatorOfParameter
  elseif ((discipl.eq.2).and.(parCat.eq.0).and.(parNum.eq.0).and. &
       (typSurf.eq.1)) then ! LSM
    isec1(6)=172         ! indicatorOfParameter
  elseif (parNum.eq.152) then 
      isec1(6)=152         ! avoid warning for lnsp       
  else
    print*,'***WARNING: undefined GRiB2 message found!',discipl, &
         parCat,parNum,typSurf
  endif
  if(parId .ne. isec1(6) .and. parId .ne. 77) then !added by mc to make it consisitent with new readwind.f90
    write(*,*) 'parId',parId, 'isec1(6)',isec1(6)  !
!    stop
  endif

  endif

  !HSO  get the size and data of the values array
  if (isec1(6).ne.-1) then
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

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
 if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then ! !added by mc to make it consisitent with new readwind.f90
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
        if(isec1(6).eq.134) psn(i,j,1,n,l)= &!! SURF. PRESS.
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.135) wwhn(i,j,nlev_ec-k+1,l)= &!! W VELOCITY
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
        if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
        if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

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


      end do
    end do

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
