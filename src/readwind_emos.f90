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

subroutine readwind(indj,n,uuh,vvh,wwh)

  !**********************************************************************
  !                                                                     *
  !             TRAJECTORY MODEL SUBROUTINE READWIND                    *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  !             AUTHOR:      G. WOTAWA                                  *
  !             DATE:        1997-08-05                                 *
  !             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
  !                                                                     *
  !**********************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) in common block
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

  use par_mod
  use com_mod

  implicit none

  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,lunit

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec0(2),isec1(56),isec2(22+nxmax+nymax),isec3(2)
  integer :: isec4(64),inbuff(jpack),ilen,ierr,iword
  !integer iswap
  real :: zsec2(60+2*nuvzmax),zsec3(2),zsec4(jpunp)
  real :: xaux,yaux,xaux0,yaux0
  real,parameter :: eps=1.e-4
  real :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

  character(len=1) :: yoper = 'D'
  logical :: hflswitch,strswitch

  hflswitch=.false.
  strswitch=.false.
  levdiff2=nlev_ec-nwz+1
  iumax=0
  iwmax=0

  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
5   call pbopen(lunit,path(3)(1:length(3))//wfname(indj),'r',ierr)
  if(ierr.lt.0) goto 999

  ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call pbgrib(lunit,inbuff,jpack,ilen,ierr)
  if(ierr.eq.-1) goto 50    ! EOF DETECTED
  if(ierr.lt.-1) goto 888   ! ERROR DETECTED

  ierr=1

  ! Check whether we are on a little endian or on a big endian computer
  !********************************************************************

  !if (inbuff(1).eq.1112101447) then         ! little endian, swap bytes
  !  iswap=1+ilen/4
  !  call swap32(inbuff,iswap)
  !else if (inbuff(1).ne.1196575042) then    ! big endian
  !  stop 'subroutine gridcheck: corrupt GRIB data'
  !endif


  call gribex(isec0,isec1,isec2,zsec2,isec3,zsec3,isec4, &
       zsec4,jpunp,inbuff,jpack,iword,yoper,ierr)
  if (ierr.ne.0) goto 888   ! ERROR DETECTED

  if(ifield.eq.1) then

  ! CHECK GRID SPECIFICATIONS

    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    if(isec2(12)/2-1.ne.nlev_ec) &
         stop 'READWIND: VERTICAL DISCRETIZATION NOT CONSISTENT'
    xaux=real(isec2(5))/1000.+real(nxshift)*dx
    yaux=real(isec2(7))/1000.
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

  do j=0,nymin1
    do i=0,nxfield-1
      k=isec1(8)
      if(isec1(6).eq.130) tth(i,j,nlev_ec-k+2,n)= &!! TEMPERATURE
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.131) uuh(i,j,nlev_ec-k+2)= &!! U VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.132) vvh(i,j,nlev_ec-k+2)= &!! V VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.133) then                      !! SPEC. HUMIDITY
        qvh(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (qvh(i,j,nlev_ec-k+2,n) .lt. 0.) &
             qvh(i,j,nlev_ec-k+2,n) = 0.
  !        this is necessary because the gridded data may contain
  !        spurious negative values
      endif
      if(isec1(6).eq.134) ps(i,j,1,n)= &!! SURF. PRESS.
           zsec4(nxfield*(ny-j-1)+i+1)

      if(isec1(6).eq.135) wwh(i,j,nlev_ec-k+1)= &!! W VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.141) sd(i,j,1,n)= &!! SNOW DEPTH
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.151) msl(i,j,1,n)= &!! SEA LEVEL PRESS.
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.164) tcc(i,j,1,n)= &!! CLOUD COVER
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.165) u10(i,j,1,n)= &!! 10 M U VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.166) v10(i,j,1,n)= &!! 10 M V VELOCITY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.167) tt2(i,j,1,n)= &!! 2 M TEMPERATURE
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.168) td2(i,j,1,n)= &!! 2 M DEW POINT
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.142) then                      !! LARGE SCALE PREC.
        lsprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (lsprec(i,j,1,n).lt.0.) lsprec(i,j,1,n)=0.
      endif
      if(isec1(6).eq.143) then                      !! CONVECTIVE PREC.
        convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
      endif
      if(isec1(6).eq.146) sshf(i,j,1,n)= &!! SENS. HEAT FLUX
           zsec4(nxfield*(ny-j-1)+i+1)
      if((isec1(6).eq.146).and.(zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) &
           hflswitch=.true.    ! Heat flux available
      if(isec1(6).eq.176) then                      !! SOLAR RADIATION
        ssr(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
        if (ssr(i,j,1,n).lt.0.) ssr(i,j,1,n)=0.
      endif
      if(isec1(6).eq.180) ewss(i,j)= &!! EW SURFACE STRESS
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.181) nsss(i,j)= &!! NS SURFACE STRESS
           zsec4(nxfield*(ny-j-1)+i+1)
      if(((isec1(6).eq.180).or.(isec1(6).eq.181)).and. &
           (zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) strswitch=.true.    ! stress available
      if(isec1(6).eq.129) oro(i,j)= &!! ECMWF OROGRAPHY
           zsec4(nxfield*(ny-j-1)+i+1)/ga
      if(isec1(6).eq.160) excessoro(i,j)= &!! STANDARD DEVIATION OF OROGRAPHY
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.172) lsm(i,j)= &!! ECMWF LAND SEA MASK
           zsec4(nxfield*(ny-j-1)+i+1)
      if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
      if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)
    end do
  end do

  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
50   call pbclose(lunit,ierr)     !! FINNISHED READING / CLOSING GRIB FILE

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
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  stop 'Execution terminated'
999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED !!!                    #### '
  stop 'Execution terminated'

end subroutine readwind
