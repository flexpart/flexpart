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
  integer :: ngrid,itage,nage,hz,il,interp_time, n
  integer(kind=1) :: clouds_v
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
  logical :: readclouds_this_nest


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
    readclouds_this_nest=.false.

    if (ngrid.gt.0) then
      xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
      ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
      if (readclouds_nest(ngrid)) readclouds_this_nest=.true.
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
!        goto 26
        exit
      endif
    end do
!26  continue

    n=memind(2)
    if (abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) &
         n=memind(1)

    if (ngrid.eq.0) then
      clouds_v=clouds(ix,jy,hz,n)
      clouds_h=cloudsh(ix,jy,n)
    else
      clouds_v=cloudsn(ix,jy,hz,n,ngrid)
      clouds_h=cloudshn(ix,jy,n,ngrid)
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

      if (ngrid.gt.0) then
        act_temp=ttn(ix,jy,hz,n,ngrid)
      else
        act_temp=tt(ix,jy,hz,n)
      endif


!***********************
! BELOW CLOUD SCAVENGING
!***********************  
      if (clouds_v.ge.4) then !below cloud

! For gas: if positive below-cloud parameters (A or B), and dquer<=0
!******************************************************************
        if ((dquer(ks).le.0.).and.(weta_gas(ks).gt.0..or.wetb_gas(ks).gt.0.)) then
          !        if (weta(ks).gt.0. .or. wetb(ks).gt.0.) then 
          blc_count=blc_count+1
          wetscav=weta_gas(ks)*prec(1)**wetb_gas(ks)

! For aerosols: if positive below-cloud parameters (Crain/Csnow or B), and dquer>0
!*********************************************************************************
        else if ((dquer(ks).gt.0.).and.(crain_aero(ks).gt.0..or.csnow_aero(ks).gt.0.)) then
          blc_count=blc_count+1

!NIK 17.02.2015
! For the calculation here particle size needs to be in meter and not um as dquer is
! changed in readreleases
! For particles larger than 10 um use the largest size defined in the parameterizations (10um)
          dquer_m=min(10.,dquer(ks))/1000000. !conversion from um to m

! Rain:
          if (act_temp .ge. 273. .and. crain_aero(ks).gt.0.)  then

! ZHG 2014 : Particle RAIN scavenging coefficient based on Laakso et al 2003, 
! the below-cloud scavenging (rain efficienty) parameter Crain (=crain_aero) from SPECIES file
            wetscav=crain_aero(ks)*10**(bclr(1)+(bclr(2)*(log10(dquer_m))**(-4))+ &
                 & (bclr(3)*(log10(dquer_m))**(-3))+ (bclr(4)*(log10(dquer_m))**(-2))+&
                 &(bclr(5)*(log10(dquer_m))**(-1))+bclr(6)* (prec(1))**(0.5))

! Snow:
          elseif (act_temp .lt. 273. .and. csnow_aero(ks).gt.0.)  then 
! ZHG 2014 : Particle SNOW scavenging coefficient based on Kyro et al 2009, 
! the below-cloud scavenging (Snow efficiency) parameter Csnow (=csnow_aero) from SPECIES file
            wetscav=csnow_aero(ks)*10**(bcls(1)+(bcls(2)*(log10(dquer_m))**(-4))+&
                 &(bcls(3)*(log10(dquer_m))**(-3))+ (bcls(4)*(log10(dquer_m))**(-2))+&
                 &(bcls(5)*(log10(dquer_m))**(-1))+ bcls(6)* (prec(1))**(0.5))

          endif
          
!             write(*,*) 'bl-cloud, act_temp=',act_temp, ',prec=',prec(1),',wetscav=', wetscav, ', jpart=',jpart

        endif ! gas or particle
!      endif ! positive below-cloud scavenging parameters given in Species file
      endif !end BELOW

!********************
! IN CLOUD SCAVENGING
!********************
      if (clouds_v.lt.4) then ! In-cloud
! NIK 13 may 2015: only do incloud if positive in-cloud scavenging parameters are
! given in species file, or if gas and positive Henry's constant
        if ((ccn_aero(ks).gt.0. .or. in_aero(ks).gt.0.).or.(henry(ks).gt.0.and.dquer(ks).le.0)) then 
          inc_count=inc_count+1
! if negative coefficients (turned off) set to zero for use in equation
          if (ccn_aero(ks).lt.0.) ccn_aero(ks)=0.
          if (in_aero(ks).lt.0.) in_aero(ks)=0.

!ZHG 2015 Cloud liquid & ice water (CLWC+CIWC) from ECMWF
! nested fields
          if (ngrid.gt.0.and.readclouds_this_nest) then
            cl=ctwcn(ix,jy,n,ngrid)*(grfraction(1)/cc)
          else if (ngrid.eq.0.and.readclouds) then
            cl=ctwc(ix,jy,n)*(grfraction(1)/cc)
          else                                  !parameterize cloudwater m2/m3
!ZHG updated parameterization of cloud water to better reproduce the values coming from ECMWF
            cl=1.6E-6*prec(1)**0.36
          endif

!ZHG: Calculate the partition between liquid and water phase water. 
          if (act_temp .le. 253.) then
            liq_frac=0
          else if (act_temp .ge. 273.) then
            liq_frac=1
          else
            liq_frac =((act_temp-273.)/(273.-253.))**2.
          end if
! ZHG: Calculate the aerosol partition based on cloud phase and Ai and Bi
          frac_act = liq_frac*ccn_aero(ks) +(1-liq_frac)*in_aero(ks)

!ZHG Use the activated fraction and the liqid water to calculate the washout

! AEROSOL
!**************************************************
          if (dquer(ks).gt.0.) then ! is particle

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
          if ((readclouds.and.ngrid.eq.0).or.(readclouds_this_nest.and.ngrid.gt.0)) then
            wetscav=incloud_ratio*S_i*(prec(1)/3.6E6)
          else
            wetscav=incloud_ratio*S_i*(prec(1)/3.6E6)/clouds_h
          endif


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
