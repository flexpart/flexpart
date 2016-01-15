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

subroutine wetdepo(itime,ltsample,loutnext)
  !                  i      i        i
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! cc [0-1]           total cloud cover                                       *
  ! convp [mm/h]       convective precipitation rate                           *
  ! grfraction [0-1]   fraction of grid, for which precipitation occurs        *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! ldeltat [s]        interval since radioactive decay was computed           *
  ! lfr, cfr           area fraction covered by precipitation for large scale  *
  !                    and convective precipitation (dependent on prec. rate)  *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! lsp [mm/h]         large scale precipitation rate                          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! prec [mm/h]        precipitation rate in subgrid, where precipitation occurs*
  ! wetdeposit         mass that is wet deposited                              *
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat,i,j,ix,jy
  integer :: ngrid,itage,nage,hz,il,interp_time, n, clouds_v
  integer :: ks, kp
!  integer :: n1,n2, icbot,ictop, indcloud !TEST
  real :: S_i, act_temp, cl, cle ! in cloud scavenging
  real :: clouds_h ! cloud height for the specific grid point
  real :: xtn,ytn,lsp,convp,cc,grfraction(3),prec(3),wetscav,totprec
  real :: wetdeposit(maxspec),restmass
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  !save lfr,cfr

  real, parameter :: lfr(5) = (/ 0.5,0.65,0.8,0.9,0.95/)
  real, parameter :: cfr(5) = (/ 0.4,0.55,0.7,0.8,0.9 /)

!ZHG aerosol below-cloud scavenging removal polynomial constants for rain and snow
  real, parameter :: bclr(6) = (/274.35758, 332839.59273, 226656.57259, 58005.91340, 6588.38582, 0.244984/) !rain (Laakso et al 2003)
  real, parameter :: bcls(6) = (/22.7, 0.0, 0.0, 1321.0, 381.0, 0.0/) !now (Kyro et al 2009)
  real :: frac_act, liq_frac, dquer_m

  integer :: blc_count, inc_count
  real    :: Si_dummy, wetscav_dummy



  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

  ! Loop over all particles
  !************************

  blc_count=0
  inc_count=0

  do jpart=1,numpart

    if (itra1(jpart).eq.-999999999) goto 20
    if(ldirect.eq.1)then
      if (itra1(jpart).gt.itime) goto 20
    else
      if (itra1(jpart).lt.itime) goto 20
    endif

  ! Determine age class of the particle
    itage=abs(itra1(jpart)-itramem(jpart))
    do nage=1,nageclass
      if (itage.lt.lage(nage)) goto 33
    end do
33  continue


  ! Determine which nesting level to be used
  !*****************************************

    ngrid=0
    do j=numbnests,1,-1
      if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and. &
      (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
        ngrid=j
        goto 23
      endif
    end do
23  continue


  ! Determine nested grid coordinates
  !**********************************

    if (ngrid.gt.0) then
      xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
      ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
    else
      ix=int(xtra1(jpart))
      jy=int(ytra1(jpart))
    endif


  ! Interpolate large scale precipitation, convective precipitation and
  ! total cloud cover
  ! Note that interpolated time refers to itime-0.5*ltsample [PS]
  !********************************************************************
    interp_time=nint(itime-0.5*ltsample)
    
    if (ngrid.eq.0) then
      call interpol_rain(lsprec,convprec,tcc,nxmax,nymax, &
      1,nx,ny,memind,real(xtra1(jpart)),real(ytra1(jpart)),1, &
      memtime(1),memtime(2),interp_time,lsp,convp,cc)
    else
      call interpol_rain_nests(lsprecn,convprecn,tccn, &
      nxmaxn,nymaxn,1,maxnests,ngrid,nxn,nyn,memind,xtn,ytn,1, &
      memtime(1),memtime(2),interp_time,lsp,convp,cc)
    endif

!  If total precipitation is less than 0.01 mm/h - no scavenging occurs
    if ((lsp.lt.0.01).and.(convp.lt.0.01)) goto 20

  ! get the level were the actual particle is in
    do il=2,nz
      if (height(il).gt.ztra1(jpart)) then
        hz=il-1
        goto 26
      endif
    end do
26  continue

    n=memind(2)
    if (abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) &
         n=memind(1)

    if (ngrid.eq.0) then
      clouds_v=clouds(ix,jy,hz,n)
      clouds_h=cloudsh(ix,jy,n)
    else
      ! new removal not implemented for nests yet 
      clouds_v=cloudsn(ix,jy,hz,n,ngrid)
      clouds_h=cloudsnh(ix,jy,n,ngrid)
    endif

  ! if there is no precipitation or the particle is above the clouds no
  ! scavenging is done

    if (clouds_v.le.1) goto 20
  
  ! 1) Parameterization of the the area fraction of the grid cell where the
  !    precipitation occurs: the absolute limit is the total cloud cover, but
  !    for low precipitation rates, an even smaller fraction of the grid cell
  !    is used. Large scale precipitation occurs over larger areas than
  !    convective precipitation.
  !**************************************************************************

    if (lsp.gt.20.) then
      i=5
    else if (lsp.gt.8.) then
      i=4
    else if (lsp.gt.3.) then
      i=3
    else if (lsp.gt.1.) then
      i=2
    else
      i=1
    endif

    if (convp.gt.20.) then
      j=5
    else if (convp.gt.8.) then
      j=4
    else if (convp.gt.3.) then
      j=3
    else if (convp.gt.1.) then
      j=2
    else
      j=1
    endif


  !ZHG oct 2014 : Calculated for 1) both 2) lsp 3) convp 
  ! Tentatively differentiate the grfraction for lsp and convp for treating differently the two forms
  ! for now they are treated the same
    grfraction(1)=max(0.05,cc*(lsp*lfr(i)+convp*cfr(j))/(lsp+convp))
    grfraction(2)=max(0.05,cc*(lfr(i)))
    grfraction(3)=max(0.05,cc*(cfr(j)))


  ! 2) Computation of precipitation rate in sub-grid cell
  !******************************************************
    prec(1)=(lsp+convp)/grfraction(1)
    prec(2)=(lsp)/grfraction(2)
    prec(3)=(convp)/grfraction(3)


  ! 3) Computation of scavenging coefficients for all species
  !    Computation of wet deposition
  !**********************************************************

    do ks=1,nspec      ! loop over species
      wetdeposit(ks)=0. 
      wetscav=0.   

      !ZHG test if it nested?
      if (ngrid.gt.0) then
        act_temp=ttn(ix,jy,hz,n,ngrid)
      else
        act_temp=tt(ix,jy,hz,n)
      endif
     

  !***********************
  ! BELOW CLOUD SCAVENGING
  !***********************  
      if (clouds_v.ge.4) then !below cloud

        if (weta(ks).gt.0. .or. wetb(ks).gt.0.) then !if positive below-cloud parameters given in SPECIES file (either A or B)
          blc_count=blc_count+1

!GAS
          if (dquer(ks) .le. 0.) then  !is gas
            wetscav=weta(ks)*prec(1)**wetb(ks)
            
!AEROSOL
          else !is particle
!NIK 17.02.2015
! For the calculation here particle size needs to be in meter and not um as dquer is changed to in readreleases
! for particles larger than 10 um use the largest size defined in the parameterizations (10um)
            dquer_m=min(10.,dquer(ks))/1000000. !conversion from um to m
            if (act_temp .ge. 273 .and. weta(ks).gt.0.)  then !Rain 
              ! ZHG 2014 : Particle RAIN scavenging coefficient based on Laakso et al 2003, 
              ! the below-cloud scavenging (rain efficienty) 
              ! parameter A (=weta) from SPECIES file
              wetscav= weta(ks)*10**(bclr(1)+ (bclr(2)*(log10(dquer_m))**(-4))+(bclr(3)*(log10(dquer_m))**(-3))+ (bclr(4)* &
                   (log10(dquer_m))**(-2))+ (bclr(5)*(log10(dquer_m))**(-1))+ bclr(6)* (prec(1))**(0.5))

            elseif (act_temp .lt. 273 .and. wetb(ks).gt.0.)  then ! Snow
              ! ZHG 2014 : Particle SNOW scavenging coefficient based on Kyro et al 2009, 
              ! the below-cloud scavenging (Snow efficiency) 
              ! parameter B (=wetb) from SPECIES file
              wetscav= wetb(ks)*10**(bcls(1)+ (bcls(2)*(log10(dquer_m))**(-4))+(bcls(3)*(log10(dquer_m))**(-3))+ (bcls(4)* &
                   (log10(dquer_m))**(-2))+ (bcls(5)*(log10(dquer_m))**(-1))+ bcls(6)* (prec(1))**(0.5))

            endif

!             write(*,*) 'bl-cloud, act_temp=',act_temp, ',prec=',prec(1),',wetscav=', wetscav, ', jpart=',jpart

          endif !gas or particle
        endif ! positive below-cloud scavenging parameters given in Species file
      endif !end BELOW


  !********************
  ! IN CLOUD SCAVENGING
      !******************************************************
      if (clouds_v.lt.4) then ! In-cloud
! NIK 13 may 2015: only do incloud if positive in-cloud scavenging parameters are given in species file
        if (weta_in(ks).gt.0. .or. wetb_in(ks).gt.0.) then 
          inc_count=inc_count+1
! if negative coefficients (turned off) set to zero for use in equation
          if (weta_in(ks).lt.0.) weta_in(ks)=0.
          if (wetb_in(ks).lt.0.) wetb_in(ks)=0.

          !ZHG 2015 Cloud liquid & ice water (CLWC+CIWC) from ECMWF
          if (readclouds) then                  !icloud_stats(ix,jy,4,n) has units kg/m2
!            cl =icloud_stats(ix,jy,4,n)*(grfraction(1)/cc)
! ESO: stop using icloud_stats
            cl =clw4(ix,jy,n)*(grfraction(1)/cc)
          else                                  !parameterize cloudwater m2/m3
            !ZHG updated parameterization of cloud water to better reproduce the values coming from ECMWF
            cl=1.6E-6*prec(1)**0.36
          endif

            !ZHG: Calculate the partition between liquid and water phase water. 
            if (act_temp .le. 253) then
              liq_frac=0
            else if (act_temp .ge. 273) then
              liq_frac=1
            else
              liq_frac =((act_temp-273)/(273-253))**2
            endif
           ! ZHG: Calculate the aerosol partition based on cloud phase and Ai and Bi
            frac_act = liq_frac*weta_in(ks) +(1-liq_frac)*wetb_in(ks)
 
  !ZHG Use the activated fraction and the liqid water to calculate the washout

  ! AEROSOL
          !**************************************************
          if (dquer(ks).gt. 0.) then ! is particle

            S_i= frac_act/cl

          !*********************
  ! GAS
          else ! is gas
               
            cle=(1-cl)/(henry(ks)*(r_air/3500.)*act_temp)+cl
            !REPLACE to switch old/ new scheme 
            ! S_i=frac_act/cle
            S_i=1/cle
          endif ! gas or particle

  ! scavenging coefficient based on Hertel et al 1995 - using the S_i for either gas or aerosol
           !OLD 
          if (readclouds) then
            wetscav=S_i*(prec(1)/3.6E6)
          else
            wetscav=S_i*(prec(1)/3.6E6)/clouds_h
          endif

!ZHG 2015 TEST
!          Si_dummy=frac_act/2E-7*prec(1)**0.36
!           wetscav_dummy=Si_dummy*(prec(1)/3.6E6)/clouds_h
!           if (clouds_v.lt.4) then
!           talltest=talltest+1
!if(talltest .eq. 1) OPEN(UNIT=199, FILE=utfil,FORM='FORMATTED',STATUS = 'UNKNOWN')
!if(talltest .lt. 100001)  write(199,*) prec(1)/3.6E6, cl, clouds_h*2E-7*prec(1)**0.36,clouds_v,ytra1(jpart)-90
!if(talltest .lt. 100001)  write(199,*) wetscav, wetscav_dummy,prec(1),ytra1(jpart)-90,clouds_v,cl
!if(talltest .eq. 100001) CLOSE(199)
!if(talltest .eq. 100001) STOP
!
!write(*,*)  'PREC kg/m2s CLOUD kg/m2', (prec(1)/3.6E6), cl !, '2E-7*prec(1)**0.36',  2E-7*prec(1)**0.36,'2E-7*prec(1)**0.36*clouds_h',2E-7*prec(1)**0.36*clouds_h
!write(*,*)  'PREC kg/m2s LSP+convp kg/m2', prec(1), convp+lsp
!write(*,*)  wetscav, wetscav_dummy
!write(*,*) cc, grfraction(1), cc/grfraction(1)

!write(*,*)  'Lmbda_old', (prec(1)/3.6E6)/(clouds_h*2E-7*prec(1)**0.36)


!write(*,*) '**************************************************'
!write(*,*)  'clouds_h', clouds_h, 'clouds_v',clouds_v,'abs(ltsample)', abs(ltsample)
!write(*,*)  'readclouds', readclouds, 'wetscav',wetscav, 'wetscav_dummy', wetscav_dummy
!write(*,*)  'S_i', S_i , 'Si_dummy', Si_dummy, 'prec(1)', prec(1) 


!           write(*,*) 'PRECIPITATION ,cl  ECMWF , cl PARAMETIZED, clouds_v, lat' &
!                      ,prec(1)/3.6E6, cl, clouds_h*2E-7*prec(1)**0.36,clouds_v,ytra1(jpart)-90

!endif 

        endif ! positive in-cloud scavenging parameters given in Species file
      endif !incloud
!END ZHG TEST
     
  !**************************************************
  ! CALCULATE DEPOSITION 
  !**************************************************
  !     if (wetscav.le.0) write (*,*) 'neg, or 0 wetscav!'
  !     +          ,wetscav,cle,cl,act_temp,prec,clouds_h,clouds_v

      if (wetscav.gt.0.) then
        wetdeposit(ks)=xmass1(jpart,ks)* &
             (1.-exp(-wetscav*abs(ltsample)))*grfraction(1)  ! wet deposition
!write(*,*) 'MASS DEPOSITED: PREC, WETSCAV, WETSCAVP', prec(1), wetdeposit(ks), xmass1(jpart,ks)* & 
!             (1.-exp(-wetscav_dummy*abs(ltsample)))*grfraction(1), clouds_v


      else ! if no scavenging
        wetdeposit(ks)=0.
      endif

      restmass = xmass1(jpart,ks)-wetdeposit(ks)
      if (ioutputforeachrelease.eq.1) then
        kp=npoint(jpart)
      else
        kp=1
      endif
      if (restmass .gt. smallnum) then
        xmass1(jpart,ks)=restmass
  !   depostatistic
  !   wetdepo_sum(ks,kp)=wetdepo_sum(ks,kp)+wetdeposit(ks)
  !   depostatistic
      else
        xmass1(jpart,ks)=0.
      endif
  !   Correct deposited mass to the last time step when radioactive decay of
  !   gridded deposited mass was calculated
      if (decay(ks).gt.0.) then
        wetdeposit(ks)=wetdeposit(ks)*exp(abs(ldeltat)*decay(ks))
      endif


    end do !all species

  ! Sabine Eckhardt, June 2008 create deposition runs only for forward runs
  ! Add the wet deposition to accumulated amount on output grid and nested output grid
  !*****************************************************************************

    if (ldirect.eq.1) then
      call wetdepokernel(nclass(jpart),wetdeposit,real(xtra1(jpart)), &
           real(ytra1(jpart)),nage,kp)
      if (nested_output.eq.1) call wetdepokernel_nest(nclass(jpart), &
           wetdeposit,real(xtra1(jpart)),real(ytra1(jpart)),nage,kp)
    endif

20  continue
  end do ! all particles

  ! count the total number of below-cloud and in-cloud occurences:
  tot_blc_count=tot_blc_count+blc_count
  tot_inc_count=tot_inc_count+inc_count

end subroutine wetdepo
