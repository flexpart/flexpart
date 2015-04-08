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
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tthn and qvhn (on eta coordinates) in common block
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,lunit,l

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec0(2),isec1(56),isec2(22+nxmaxn+nymaxn),isec3(2)
  integer :: isec4(64),inbuff(jpack),ilen,ierr,iword
  !integer iswap
  real :: zsec2(60+2*nuvzmax),zsec3(2),zsec4(jpunp)
  real :: xaux,yaux,xaux0,yaux0
  real :: ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

  character(len=1) :: yoper = 'D'
  logical :: hflswitch,strswitch


  do l=1,numbnests
    hflswitch=.false.
    strswitch=.false.
    levdiff2=nlev_ec-nwz+1
    iumax=0
    iwmax=0

  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
5   call pbopen(lunit,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//wfnamen(l,indj),'r',ierr)
    if(ierr.lt.0) goto 999

    ifield=0
10    ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
    call pbgrib(lunit,inbuff,jpack,ilen,ierr)
    if(ierr.eq.-1) goto 50    ! EOF DETECTED
    if(ierr.lt.-1) goto 888   ! ERROR DETECTED

    ierr=1

  ! Check whether we are on a little endian or on a big endian computer
  !********************************************************************

  !  if (inbuff(1).eq.1112101447) then         ! little endian, swap bytes
  !    iswap=1+ilen/4
  !    call swap32(inbuff,iswap)
  !  else if (inbuff(1).ne.1196575042) then    ! big endian
  !    stop 'subroutine gridcheck: corrupt GRIB data'
  !  endif

    call gribex(isec0,isec1,isec2,zsec2,isec3,zsec3,isec4, &
         zsec4,jpunp,inbuff,jpack,iword,yoper,ierr)
    if (ierr.ne.0) goto 888   ! ERROR DETECTED

    if(ifield.eq.1) then

  ! CHECK GRID SPECIFICATIONS

      if(isec2(2).ne.nxn(l)) stop &
           'READWIND: NX NOT CONSISTENT FOR A NESTING LEVEL'
      if(isec2(3).ne.nyn(l)) stop &
           'READWIND: NY NOT CONSISTENT FOR A NESTING LEVEL'
      if(isec2(12)/2-1.ne.nlev_ec) stop 'READWIND: VERTICAL DISCRET&
           &IZATION NOT CONSISTENT FOR A NESTING LEVEL'
      xaux=real(isec2(5))/1000.
      yaux=real(isec2(7))/1000.
      xaux0=xlon0n(l)
      yaux0=ylat0n(l)
      if(xaux.lt.0.) xaux=xaux+360.
      if(yaux.lt.0.) yaux=yaux+360.
      if(xaux0.lt.0.) xaux0=xaux0+360.
      if(yaux0.lt.0.) yaux0=yaux0+360.
      if(xaux.ne.xaux0) &
           stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT FOR A NES&
           &TING LEVEL'
      if(yaux.ne.yaux0) &
           stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT FOR A NEST&
           &ING LEVEL'
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
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
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
          convprecn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
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
      end do
    end do

    goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
50   call pbclose(lunit,ierr)     !! FINISHED READING / CLOSING GRIB FILE

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
