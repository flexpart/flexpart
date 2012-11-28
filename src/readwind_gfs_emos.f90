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
  !             CAHENGE: 16/11/2005, Caroline Forster, GFS data         *
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
  integer :: ii,indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,lunit

  ! NCEP

  integer :: numpt,numpu,numpv,numpw,numprh
  real :: help, temp, ew
  real :: elev
  real :: ulev1(0:nxmax-1,0:nymax-1),vlev1(0:nxmax-1,0:nymax-1)
  real :: tlev1(0:nxmax-1,0:nymax-1)
  real :: qvh2(0:nxmax-1,0:nymax-1)

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec0(2),isec1(56),isec2(22+nxmax+nymax),isec3(2)
  integer :: isec4(64),inbuff(jpack),ilen,iswap,ierr,iword
  real :: zsec2(60+2*nuvzmax),zsec3(2),zsec4(jpunp)
  real :: xaux,yaux,xaux0,yaux0
  real,parameter :: eps=1.e-4
  real :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
  real :: plev1,hlev1,ff10m,fflev1

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
  if (ierr.ne.0) goto 10   ! ERROR DETECTED

  if(ifield.eq.1) then

  ! CHECK GRID SPECIFICATIONS

    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    xaux=real(isec2(5))/1000.+real(nxshift)*dx
    if(xaux.eq.0.) xaux=-179.0     ! NCEP DATA
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
      if((isec1(6).eq.011).and.(isec1(7).eq.100)) then
  ! TEMPERATURE
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpt=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          tth(179+i,j,numpt,n)=help
        else
          tth(i-181,j,numpt,n)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.100)) then
  ! U VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpu=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          uuh(179+i,j,numpu)=help
        else
          uuh(i-181,j,numpu)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.100)) then
  ! V VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpv=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          vvh(179+i,j,numpv)=help
        else
          vvh(i-181,j,numpv)=help
        endif
      endif
      if((isec1(6).eq.052).and.(isec1(7).eq.100)) then
  ! RELATIVE HUMIDITY -> CONVERT TO SPECIFIC HUMIDITY LATER
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numprh=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          qvh(179+i,j,numprh,n)=help
        else
          qvh(i-181,j,numprh,n)=help
        endif
      endif
      if((isec1(6).eq.001).and.(isec1(7).eq.001)) then
  ! SURFACE PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          ps(179+i,j,1,n)=help
        else
          ps(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.039).and.(isec1(7).eq.100)) then
  ! W VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpw=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          wwh(179+i,j,numpw)=help
        else
          wwh(i-181,j,numpw)=help
        endif
      endif
      if((isec1(6).eq.066).and.(isec1(7).eq.001)) then
  ! SNOW DEPTH
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          sd(179+i,j,1,n)=help
        else
          sd(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.002).and.(isec1(7).eq.102)) then
  ! MEAN SEA LEVEL PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          msl(179+i,j,1,n)=help
        else
          msl(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.071).and.(isec1(7).eq.244)) then
  ! TOTAL CLOUD COVER
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          tcc(179+i,j,1,n)=help
        else
          tcc(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.10)) then
  ! 10 M U VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          u10(179+i,j,1,n)=help
        else
          u10(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.10)) then
  ! 10 M V VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          v10(179+i,j,1,n)=help
        else
          v10(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.011).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          tt2(179+i,j,1,n)=help
        else
          tt2(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.017).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M DEW POINT TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          td2(179+i,j,1,n)=help
        else
          td2(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.062).and.(isec1(7).eq.001)) then
  ! LARGE SCALE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          lsprec(179+i,j,1,n)=help
        else
          lsprec(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.063).and.(isec1(7).eq.001)) then
  ! CONVECTIVE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          convprec(179+i,j,1,n)=help
        else
          convprec(i-181,j,1,n)=help
        endif
      endif
  ! SENS. HEAT FLUX
      sshf(i,j,1,n)=0.0     ! not available from gfs.tccz.pgrbfxx files
      hflswitch=.false.    ! Heat flux not available
  ! SOLAR RADIATIVE FLUXES
      ssr(i,j,1,n)=0.0      ! not available from gfs.tccz.pgrbfxx files
  ! EW SURFACE STRESS
      ewss(i,j)=0.0         ! not available from gfs.tccz.pgrbfxx files
  ! NS SURFACE STRESS
      nsss(i,j)=0.0         ! not available from gfs.tccz.pgrbfxx files
      strswitch=.false.    ! stress not available
      if((isec1(6).eq.007).and.(isec1(7).eq.001)) then
  ! TOPOGRAPHY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          oro(179+i,j)=help
          excessoro(179+i,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        else
          oro(i-181,j)=help
          excessoro(i-181,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        endif
      endif
      if((isec1(6).eq.081).and.(isec1(7).eq.001)) then
  ! LAND SEA MASK
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          lsm(179+i,j)=help
        else
          lsm(i-181,j)=help
        endif
      endif
      if((isec1(6).eq.221).and.(isec1(7).eq.001)) then
  ! MIXING HEIGHT
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          hmix(179+i,j,1,n)=help
        else
          hmix(i-181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.052).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M RELATIVE HUMIDITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          qvh2(179+i,j)=help
        else
          qvh2(i-181,j)=help
        endif
      endif
      if((isec1(6).eq.011).and.(isec1(7).eq.107)) then
  ! TEMPERATURE LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          tlev1(179+i,j)=help
        else
          tlev1(i-181,j)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.107)) then
  ! U VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          ulev1(179+i,j)=help
        else
          ulev1(i-181,j)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.107)) then
  ! V VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.180) then
          vlev1(179+i,j)=help
        else
          vlev1(i-181,j)=help
        endif
      endif

    end do
  end do

  if((isec1(6).eq.33).and.(isec1(7).eq.100)) then
  ! NCEP ISOBARIC LEVELS
    iumax=iumax+1
  endif


  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
50   call pbclose(lunit,ierr)     !! FINNISHED READING / CLOSING GRIB FILE

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

end subroutine readwind
