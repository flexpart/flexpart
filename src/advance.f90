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

subroutine advance(itime,nrelpoint,ldt,up,vp,wp, &
       usigold,vsigold,wsigold,nstop,xt,yt,zt,prob,icbt)
  !                     i    i  i/oi/oi/o
  !  i/o     i/o     i/o     o  i/oi/oi/o i/o  i/o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of turbulent particle trajectories utilizing a                *
  !  zero-acceleration scheme, which is corrected by a numerically more        *
  !  accurate Petterssen scheme whenever possible.                             *
  !                                                                            *
  !  Particle positions are read in, incremented, and returned to the calling  *
  !  program.                                                                  *
  !                                                                            *
  !  In different regions of the atmosphere (PBL vs. free troposphere),        *
  !  different parameters are needed for advection, parameterizing turbulent   *
  !  velocities, etc. For efficiency, different interpolation routines have    *
  !  been written for these different cases, with the disadvantage that there  *
  !  exist several routines doing almost the same. They all share the          *
  !  included file 'interpol_mod'. The following                               *
  !  interpolation routines are used:                                          *
  !                                                                            *
  !  interpol_all(_nests)     interpolates everything (called inside the PBL)  *
  !  interpol_misslev(_nests) if a particle moves vertically in the PBL,       *
  !                           additional parameters are interpolated if it     *
  !                           crosses a model level                            *
  !  interpol_wind(_nests)    interpolates the wind and determines the         *
  !                           standard deviation of the wind (called outside   *
  !                           PBL) also interpolates potential vorticity       *
  !  interpol_wind_short(_nests) only interpolates the wind (needed for the    *
  !                           Petterssen scheme)                               *
  !  interpol_vdep(_nests)    interpolates deposition velocities               *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Changes:                                                                  *
  !                                                                            *
  !  8 April 2000: Deep convection parameterization                            *
  !                                                                            *
  !  May 2002: Petterssen scheme introduced                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! icbt               1 if particle not transferred to forbidden state,       *
  !                    else -1                                                 *
  ! dawsave            accumulated displacement in along-wind direction        *
  ! dcwsave            accumulated displacement in cross-wind direction        *
  ! dxsave             accumulated displacement in longitude                   *
  ! dysave             accumulated displacement in latitude                    *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          time at which this subroutine is entered                *
  ! itimec [s]         actual time, which is incremented in this subroutine    *
  ! href [m]           height for which dry deposition velocity is calculated  *
  ! ladvance [s]       Total integration time period                           *
  ! ldirect            1 forward, -1 backward                                  *
  ! ldt [s]            Time step for the next integration                      *
  ! lsynctime [s]      Synchronisation interval of FLEXPART                    *
  ! ngrid              index which grid is to be used                          *
  ! nrand              index for a variable to be picked from rannumb          *
  ! nstop              if > 1 particle has left domain and must be stopped     *
  ! prob               probability of absorption due to dry deposition         *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! rhoa               air density                                             *
  ! rhograd            vertical gradient of the air density                    *
  ! up,vp,wp           random velocities due to turbulence (along wind, cross  *
  !                    wind, vertical wind                                     *
  ! usig,vsig,wsig     mesoscale wind fluctuations                             *
  ! usigold,vsigold,wsigold  like usig, etc., but for the last time step       *
  ! vdepo              Deposition velocities for all species                   *
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod
  use cmapf_mod
  use random_mod, only: ran3

  implicit none

  real(kind=dp) :: xt,yt
  real :: zt,xts,yts,weight
  integer :: itime,itimec,nstop,ldt,i,j,k,nrand,loop,memindnext,mind
  integer :: ngr,nix,njy,ks,nsp,nrelpoint
  real :: dz,dz1,dz2,xlon,ylat,xpol,ypol,gridsize
  real :: ru,rv,rw,dt,ux,vy,cosfact,xtn,ytn,tropop
  real :: prob(maxspec),up,vp,wp,dxsave,dysave,dawsave
  real :: dcwsave
  real :: usigold,vsigold,wsigold,r,rs
  real :: uold,vold,wold,vdepo(maxspec)
  !real uprof(nzmax),vprof(nzmax),wprof(nzmax)
  !real usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
  !real rhoprof(nzmax),rhogradprof(nzmax)
  real :: rhoa,rhograd,delz,dtf,rhoaux,dtftlw,uxscale,wpscale
  integer(kind=2) :: icbt
  real,parameter :: eps=nxmax/3.e5,eps2=1.e-9
  real :: ptot_lhh,Q_lhh,phi_lhh,ath,bth !modified by mc 
  real :: old_wp_buf,dcas,dcas1,del_test !added by mc
  integer :: i_well,jj,flagrein !test well mixed: modified by mc


  !!! CHANGE: TEST OF THE WELL-MIXED CRITERION
  !  integer,parameter :: iclass=10
  !  real(kind=dp) :: zacc,tacc,t(iclass),th(0:iclass),hsave
  !  logical dump
  !  save zacc,tacc,t,th,hsave,dump
  !!! CHANGE

  integer :: idummy = -7
  real    :: settling = 0.


  !!! CHANGE: TEST OF THE WELL-MIXED CRITERION
  !if (idummy.eq.-7) then
  !open(550,file='WELLMIXEDTEST')
  !do 17 i=0,iclass
  !7      th(i)=real(i)/real(iclass)
  !endif
  !!! CHANGE


  nstop=0
  do i=1,nmixz
    indzindicator(i)=.true.
  end do


  if (DRYDEP) then    ! reset probability for deposition
    do ks=1,nspec
      depoindicator(ks)=.true.
      prob(ks)=0.
    end do
  endif

  dxsave=0.           ! reset position displacements
  dysave=0.           ! due to mean wind
  dawsave=0.          ! and turbulent wind
  dcwsave=0.

  itimec=itime

  nrand=int(ran3(idummy)*real(maxrand-1))+1


  ! Determine whether lat/long grid or polarstereographic projection
  ! is to be used
  ! Furthermore, determine which nesting level to be used
  !*****************************************************************

  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngrid=j
        goto 23
      endif
    end do
23   continue
  endif


  !***************************
  ! Interpolate necessary data
  !***************************

  if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
    memindnext=1
  else
    memindnext=2
  endif

  ! Determine nested grid coordinates
  !**********************************

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)
  else
    ix=int(xt)
    jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
  endif
  ixp=ix+1
  jyp=jy+1


  ! Compute maximum mixing height around particle position
  !*******************************************************

  h=0.
  if (ngrid.le.0) then
    do k=1,2
      mind=memind(k) ! eso: compatibility with 3-field version
      do j=jy,jyp
        do i=ix,ixp
          if (hmix(i,j,1,mind).gt.h) h=hmix(i,j,1,mind)
        end do
      end do
    end do
    tropop=tropopause(nix,njy,1,1)
  else
    do k=1,2
      mind=memind(k)
      do j=jy,jyp
        do i=ix,ixp
          if (hmixn(i,j,1,mind,ngrid).gt.h) h=hmixn(i,j,1,mind,ngrid)
        end do
      end do
    end do
    tropop=tropopausen(nix,njy,1,1,ngrid)
  endif

  zeta=zt/h



  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************

  if (zeta.le.1.) then

  ! BEGIN TIME LOOP
  !================

    loop=0
100   loop=loop+1
      if (method.eq.1) then
        ldt=min(ldt,abs(lsynctime-itimec+itime))
        itimec=itimec+ldt*ldirect
      else
        ldt=abs(lsynctime)
        itimec=itime+lsynctime
      endif
      dt=real(ldt)

      zeta=zt/h


      if (loop.eq.1) then
        if (ngrid.le.0) then
          xts=real(xt)
          yts=real(yt)
          call interpol_all(itime,xts,yts,zt)
        else
          call interpol_all_nests(itime,xtn,ytn,zt)
        endif

      else


  ! Determine the level below the current position for u,v,rho
  !***********************************************************

        do i=2,nz
          if (height(i).gt.zt) then
            indz=i-1
            indzp=i
            goto 6
          endif
        end do
6       continue

  ! If one of the levels necessary is not yet available,
  ! calculate it
  !*****************************************************

        do i=indz,indzp
          if (indzindicator(i)) then
            if (ngrid.le.0) then
              call interpol_misslev(i)
            else
              call interpol_misslev_nests(i)
            endif
          endif
        end do
      endif


  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************

      dz=1./(height(indzp)-height(indz))
      dz1=(zt-height(indz))*dz
      dz2=(height(indzp)-zt)*dz

      u=dz1*uprof(indzp)+dz2*uprof(indz)
      v=dz1*vprof(indzp)+dz2*vprof(indz)
      w=dz1*wprof(indzp)+dz2*wprof(indz)
      rhoa=dz1*rhoprof(indzp)+dz2*rhoprof(indz)
      rhograd=dz1*rhogradprof(indzp)+dz2*rhogradprof(indz)


  ! Compute the turbulent disturbances
  ! Determine the sigmas and the timescales
  !****************************************

      if (turbswitch) then
        call hanna(zt)
      else
        call hanna1(zt)
      endif


  !*****************************************
  ! Determine the new diffusivity velocities
  !*****************************************

  ! Horizontal components
  !**********************

      if (nrand+1.gt.maxrand) nrand=1
      if (dt/tlu.lt..5) then
        up=(1.-dt/tlu)*up+rannumb(nrand)*sigu*sqrt(2.*dt/tlu)
      else
        ru=exp(-dt/tlu)
        up=ru*up+rannumb(nrand)*sigu*sqrt(1.-ru**2)
      endif
      if (dt/tlv.lt..5) then
        vp=(1.-dt/tlv)*vp+rannumb(nrand+1)*sigv*sqrt(2.*dt/tlv)
      else
        rv=exp(-dt/tlv)
        vp=rv*vp+rannumb(nrand+1)*sigv*sqrt(1.-rv**2)
      endif
      nrand=nrand+2


      if (nrand+ifine.gt.maxrand) nrand=1
      rhoaux=rhograd/rhoa
      dtf=dt*fine

      dtftlw=dtf/tlw

  ! Loop over ifine short time steps for vertical component
  !********************************************************

      do i=1,ifine

  ! Determine the drift velocity and density correction velocity
  !*************************************************************

        if (turbswitch) then
          if (dtftlw.lt..5) then
  !*************************************************************
  !************** CBL options added by mc see routine cblf90 ***
            if (cblflag.eq.1) then  !modified by mc
              if (-h/ol.gt.5) then  !modified by mc
              !if (ol.lt.0.) then   !modified by mc  
              !if (ol.gt.0.) then   !modified by mc : for test
                  !print  *,zt,wp,ath,bth,tlw,dtf,'prima'
                  flagrein=0
                  nrand=nrand+1
                  old_wp_buf=wp
                  call cbl(wp,zt,ust,wst,h,rhoa,rhograd,sigw,dsigwdz,tlw,ptot_lhh,Q_lhh,phi_lhh,ath,bth,ol,flagrein) !inside the routine for inverse time
                  wp=(wp+ath*dtf+bth*rannumb(nrand)*sqrt(dtf))*real(icbt) 
                  ! wp=(wp+ath*dtf+bth*gasdev2(mydum)*sqrt(dtf))*real(icbt) 
                  delz=wp*dtf
                  if (flagrein.eq.1) then
                      call re_initialize_particle(zt,ust,wst,h,sigw,old_wp_buf,nrand,ol)
                      wp=old_wp_buf
                      delz=wp*dtf
                      nan_count=nan_count+1
                  end if
                  !print  *,zt,wp,ath,bth,tlw,dtf,rannumb(nrand+i),icbt
                  !pause                  
              else 
                  nrand=nrand+1
                  old_wp_buf=wp
                  ath=-wp/tlw+sigw*dsigwdz+wp*wp/sigw*dsigwdz+sigw*sigw/rhoa*rhograd  !1-note for inverse time should be -wp/tlw*ldirect+... calculated for wp=-wp
                                                                                      !2-but since ldirect =-1 for inverse time and this must be calculated for (-wp) and
                                                                                      !3-the gaussian pdf is symmetric (i.e. pdf(w)=pdf(-w) ldirect can be discarded
                  bth=sigw*rannumb(nrand)*sqrt(2.*dtftlw)
                  wp=(wp+ath*dtf+bth)*real(icbt)  
                  delz=wp*dtf
                  del_test=(1.-wp)/wp !catch infinity value
                  if (isnan(wp).or.isnan(del_test)) then 
                      nrand=nrand+1                      
                      wp=sigw*rannumb(nrand)
                      delz=wp*dtf
                      nan_count2=nan_count2+1
                      !print *,'NaN coutner equal to:', nan_count,'reduce ifine if this number became a non-negligible fraction of the particle number'
                  end if  
              end if
  !******************** END CBL option *******************************            
  !*******************************************************************            
            else
                 wp=((1.-dtftlw)*wp+rannumb(nrand+i)*sqrt(2.*dtftlw) &
                 +dtf*(dsigwdz+rhoaux*sigw))*real(icbt) 
                 delz=wp*sigw*dtf
            end if
          else
            rw=exp(-dtftlw)
            wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2) &
                 +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*real(icbt)
            delz=wp*sigw*dtf
          endif
          
        else
          rw=exp(-dtftlw)
          wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw &
               +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*real(icbt)
          delz=wp*dtf
        endif

  !****************************************************
  ! Compute turbulent vertical displacement of particle
  !****************************************************

        if (abs(delz).gt.h) delz=mod(delz,h)

  ! Determine if particle transfers to a "forbidden state" below the ground
  ! or above the mixing height
  !************************************************************************

        if (delz.lt.-zt) then         ! reflection at ground
          icbt=-1
          zt=-zt-delz
        else if (delz.gt.(h-zt)) then ! reflection at h
          icbt=-1
          zt=-zt-delz+2.*h
        else                         ! no reflection
          icbt=1
          zt=zt+delz
        endif

        if (i.ne.ifine) then
          zeta=zt/h
          call hanna_short(zt)
        endif

      end do
      if (cblflag.ne.1) nrand=nrand+i

  ! Determine time step for next integration
  !*****************************************

      if (turbswitch) then
        ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5), &
             0.5/abs(dsigwdz))*ctl)
      else
        ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5))*ctl)
      endif
      ldt=max(ldt,mintime)


  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases, or if particle
  ! represents more than one species
  !*************************************************************************

      if (mdomainfill.eq.0) then
! ESO 05.2015  Changed this to fix MQUASILAG option, where nrelpoint is
!              particle number and thus xmass array goes out of bounds
!        do nsp=1,nspec
!           if (xmass(nrelpoint,nsp).gt.eps2) goto 887
!         end do
! 887     nsp=min(nsp,nspec)
        if (nspec.eq.1.and.density(1).gt.0.) then
          call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
        end if
        w=w+settling
      endif

  ! Horizontal displacements during time step dt are small real values compared
  ! to the position; adding the two, would result in large numerical errors.
  ! Thus, displacements are accumulated during lsynctime and are added to the
  ! position at the end
  !****************************************************************************

      dxsave=dxsave+u*dt
      dysave=dysave+v*dt
      dawsave=dawsave+up*dt
      dcwsave=dcwsave+vp*dt
      zt=zt+w*dt*real(ldirect)

      ! HSO/AL: Particle managed to go over highest level -> interpolation error in goto 700
      !          alias interpol_wind (division by zero)
      if (zt.ge.height(nz)) zt=height(nz)-100.*eps

      if (zt.gt.h) then
        if (itimec.eq.itime+lsynctime) goto 99
        goto 700    ! complete the current interval above PBL
      endif


  !!! CHANGE: TEST OF THE WELL-MIXED CRITERION
  !!! These lines may be switched on to test the well-mixed criterion
  !if (zt.le.h) then
  !  zacc=zacc+zt/h*dt
  !  hsave=hsave+h*dt
  !  tacc=tacc+dt
  !  do 67 i=1,iclass
  !    if ((zt/h.gt.th(i-1)).and.(zt/h.le.th(i)))
  !    +    t(i)=t(i)+dt
  !7        continue
  !endif
  !if ((mod(itime,10800).eq.0).and.dump) then
  ! dump=.false.
  ! write(550,'(i6,12f10.3)') itime,hsave/tacc,zacc/tacc,
  !    + (t(i)/tacc*real(iclass),i=1,iclass)
  !  zacc=0.
  !  tacc=0.
  !  do 68 i=1,iclass
  !8        t(i)=0.
  !  hsave=0.
  !endif
  !if (mod(itime,10800).ne.0) dump=.true.
  !!! CHANGE
      
  ! Determine probability of deposition
  !************************************

      if ((DRYDEP).and.(zt.lt.2.*href)) then
        do ks=1,nspec
          if (DRYDEPSPEC(ks)) then
            if (depoindicator(ks)) then
              if (ngrid.le.0) then
                call interpol_vdep(ks,vdepo(ks))
              else
                call interpol_vdep_nests(ks,vdepo(ks))
              endif
            endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
               prob(ks)=1.+(prob(ks)-1.)* &
                    exp(-vdepo(ks)*abs(dt)/(2.*href))
          endif
        end do
      endif

      if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection

      if (itimec.eq.(itime+lsynctime)) then
        usig=0.5*(usigprof(indzp)+usigprof(indz))
        vsig=0.5*(vsigprof(indzp)+vsigprof(indz))
        wsig=0.5*(wsigprof(indzp)+wsigprof(indz))
        goto 99  ! finished
      endif
      goto 100

  ! END TIME LOOP
  !==============


  endif



  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************


  ! Interpolate the wind
  !*********************

700   continue
  if (ngrid.le.0) then
    xts=real(xt)
    yts=real(yt)
    call interpol_wind(itime,xts,yts,zt)
  else
    call interpol_wind_nests(itime,xtn,ytn,zt)
  endif


  ! Compute everything for above the PBL

  ! Assume constant, uncorrelated, turbulent perturbations
  ! In the stratosphere, use a small vertical diffusivity d_strat,
  ! in the troposphere, use a larger horizontal diffusivity d_trop.
  ! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
  !******************************************************************

  ldt=abs(lsynctime-itimec+itime)
  dt=real(ldt)

  if (zt.lt.tropop) then  ! in the troposphere
    uxscale=sqrt(2.*d_trop/dt)
    if (nrand+1.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    nrand=nrand+2
    wp=0.
  else if (zt.lt.tropop+1000.) then     ! just above the tropopause: make transition
    weight=(zt-tropop)/1000.
    uxscale=sqrt(2.*d_trop/dt*(1.-weight))
    if (nrand+2.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    wpscale=sqrt(2.*d_strat/dt*weight)
    wp=rannumb(nrand+2)*wpscale+d_strat/1000.
    nrand=nrand+3
  else                 ! in the stratosphere
    if (nrand.gt.maxrand) nrand=1
    ux=0.
    vy=0.
    wpscale=sqrt(2.*d_strat/dt)
    wp=rannumb(nrand)*wpscale
    nrand=nrand+1
  endif


  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases
  !*************************************************************************



    if (mdomainfill.eq.0) then
! ESO 05.2015  Changed this to fix MQUASILAG option, where nrelpoint is
!              particle number and thus xmass array goes out of bounds

!      do nsp=1,nspec
!         if (xmass(nrelpoint,nsp).gt.eps2) goto 888
!       end do
! 888   nsp=min(nsp,nspec)
!        if (density(nsp).gt.0.) then
      if (nspec.eq.1.and.density(1).gt.0.) then
        call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
      end if
      w=w+settling
    endif

  ! Calculate position at time step itime+lsynctime
  !************************************************

  dxsave=dxsave+(u+ux)*dt
  dysave=dysave+(v+vy)*dt
  zt=zt+(w+wp)*dt*real(ldirect)
  if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection

99   continue



  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************


  ! Mesoscale wind velocity fluctuations are obtained by scaling
  ! with the standard deviation of the grid-scale winds surrounding
  ! the particle location, multiplied by a factor turbmesoscale.
  ! The autocorrelation time constant is taken as half the
  ! time interval between wind fields
  !****************************************************************

  r=exp(-2.*real(abs(lsynctime))/real(lwindinterv))
  rs=sqrt(1.-r**2)
  if (nrand+2.gt.maxrand) nrand=1
  usigold=r*usigold+rs*rannumb(nrand)*usig*turbmesoscale
  vsigold=r*vsigold+rs*rannumb(nrand+1)*vsig*turbmesoscale
  wsigold=r*wsigold+rs*rannumb(nrand+2)*wsig*turbmesoscale

  dxsave=dxsave+usigold*real(lsynctime)
  dysave=dysave+vsigold*real(lsynctime)

  zt=zt+wsigold*real(lsynctime)
  if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion

  !*************************************************************
  ! Transform along and cross wind components to xy coordinates,
  ! add them to u and v, transform u,v to grid units/second
  ! and calculate new position
  !*************************************************************

  call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
  dxsave=dxsave+ux   ! comment by mc: comment this line to stop the particles horizontally for test reasons 
  dysave=dysave+vy
  if (ngrid.ge.0) then
    cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
    xt=xt+real(dxsave*cosfact*real(ldirect),kind=dp)
    yt=yt+real(dysave*dyconst*real(ldirect),kind=dp)
  else if (ngrid.eq.-1) then      ! around north pole
    xlon=xlon0+xt*dx                                !comment by mc: compute old particle position
    ylat=ylat0+yt*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)   !convert old particle position in polar stereographic
    gridsize=1000.*cgszll(northpolemap,ylat,xlon)   !calculate size in m of grid element in polar stereographic coordinate
    dxsave=dxsave/gridsize                          !increment from meter to grdi unit
    dysave=dysave/gridsize
    xpol=xpol+dxsave*real(ldirect)                  !position in grid unit polar stereographic
    ypol=ypol+dysave*real(ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)  !convert to lat long coordinate
    xt=(xlon-xlon0)/dx                             !convert to grid units in lat long coordinate, comment by mc
    yt=(ylat-ylat0)/dy
  else if (ngrid.eq.-2) then    ! around south pole
    xlon=xlon0+xt*dx
    ylat=ylat0+yt*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat,xlon)
    dxsave=dxsave/gridsize
    dysave=dysave/gridsize
    xpol=xpol+dxsave*real(ldirect)
    ypol=ypol+dysave*real(ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    xt=(xlon-xlon0)/dx
    yt=(ylat-ylat0)/dy
  endif


  ! If global data are available, use cyclic boundary condition
  !************************************************************

  if (xglobal) then
    if (xt.ge.real(nxmin1)) xt=xt-real(nxmin1)
    if (xt.lt.0.) xt=xt+real(nxmin1)
    if (xt.le.eps) xt=eps
    if (abs(xt-real(nxmin1)).le.eps) xt=real(nxmin1)-eps
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************

  if ( yt.lt.0. ) then
    xt=mod(xt+180.,360.)
    yt=-yt
  else if ( yt.gt.real(nymin1) ) then
    xt=mod(xt+180.,360.)
    yt=2*real(nymin1)-yt
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************

  if ((xt.lt.0.).or.(xt.ge.real(nxmin1)).or.(yt.lt.0.).or. &
       (yt.gt.real(nymin1))) then
    nstop=3
    return
  endif

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************

  if (zt.ge.height(nz)) zt=height(nz)-100.*eps


  !************************************************************************
  ! Now we could finish, as this was done in FLEXPART versions up to 4.0.
  ! However, truncation errors of the advection can be significantly
  ! reduced by doing one iteration of the Petterssen scheme, if this is
  ! possible.
  ! Note that this is applied only to the grid-scale winds, not to
  ! the turbulent winds.
  !************************************************************************

  ! The Petterssen scheme can only applied with long time steps (only then u
  ! is the "old" wind as required by the scheme); otherwise do nothing
  !*************************************************************************

  if (ldt.ne.abs(lsynctime)) return

  ! The Petterssen scheme can only be applied if the ending time of the time step
  ! (itime+ldt*ldirect) is still between the two wind fields held in memory;
  ! otherwise do nothing
  !******************************************************************************

  if (abs(itime+ldt*ldirect).gt.abs(memtime(2))) return

  ! Apply it also only if starting and ending point of current time step are on
  ! the same grid; otherwise do nothing
  !*****************************************************************************
  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngr=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngr=-2
  else
    ngr=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngr=j
        goto 43
      endif
    end do
43   continue
  endif

  if (ngr.ne.ngrid) return

  ! Determine nested grid coordinates
  !**********************************

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
  else
    ix=int(xt)
    jy=int(yt)
  endif
  ixp=ix+1
  jyp=jy+1


  ! Memorize the old wind
  !**********************

  uold=u
  vold=v
  wold=w

  ! Interpolate wind at new position and time
  !******************************************

  if (ngrid.le.0) then
    xts=real(xt)
    yts=real(yt)
    call interpol_wind_short(itime+ldt*ldirect,xts,yts,zt)
  else
    call interpol_wind_short_nests(itime+ldt*ldirect,xtn,ytn,zt)
  endif

  if (mdomainfill.eq.0) then
! ESO 05.2015  Changed this to fix MQUASILAG option, where nrelpoint is
!              particle number and thus xmass array goes out of bounds
!    do nsp=1,nspec
!       if (xmass(nrelpoint,nsp).gt.eps2) goto 889
!     end do
! 889   nsp=min(nsp,nspec)
!      if (density(nsp).gt.0.) then
    if (nspec.eq.1.and.density(1).gt.0.) then
      call get_settling(itime+ldt,real(xt),real(yt),zt,nsp,settling)  !bugfix
    end if
    w=w+settling
  endif


  ! Determine the difference vector between new and old wind
  ! (use half of it to correct position according to Petterssen)
  !*************************************************************

  u=(u-uold)/2.
  v=(v-vold)/2.
  w=(w-wold)/2.


  ! Finally, correct the old position
  !**********************************

  zt=zt+w*real(ldt*ldirect)
  if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection
  if (ngrid.ge.0) then
    cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
    xt=xt+real(u*cosfact*real(ldt*ldirect),kind=dp)
    yt=yt+real(v*dyconst*real(ldt*ldirect),kind=dp)
  else if (ngrid.eq.-1) then      ! around north pole
    xlon=xlon0+xt*dx
    ylat=ylat0+yt*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(northpolemap,ylat,xlon)
    u=u/gridsize
    v=v/gridsize
    xpol=xpol+u*real(ldt*ldirect)
    ypol=ypol+v*real(ldt*ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
    xt=(xlon-xlon0)/dx
    yt=(ylat-ylat0)/dy
  else if (ngrid.eq.-2) then    ! around south pole
    xlon=xlon0+xt*dx
    ylat=ylat0+yt*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat,xlon)
    u=u/gridsize
    v=v/gridsize
    xpol=xpol+u*real(ldt*ldirect)
    ypol=ypol+v*real(ldt*ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    xt=(xlon-xlon0)/dx
    yt=(ylat-ylat0)/dy
  endif

  ! If global data are available, use cyclic boundary condition
  !************************************************************

  if (xglobal) then
    if (xt.ge.real(nxmin1)) xt=xt-real(nxmin1)
    if (xt.lt.0.) xt=xt+real(nxmin1)
    if (xt.le.eps) xt=eps
    if (abs(xt-real(nxmin1)).le.eps) xt=real(nxmin1)-eps
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************

  if ( yt.lt.0. ) then
    xt=mod(xt+180.,360.)
    yt=-yt
  else if ( yt.gt.real(nymin1) ) then
    xt=mod(xt+180.,360.)
    yt=2*real(nymin1)-yt
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************

  if ((xt.lt.0.).or.(xt.ge.real(nxmin1)).or.(yt.lt.0.).or. &
       (yt.gt.real(nymin1))) then
    nstop=3
    return
  endif

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************

  if (zt.ge.height(nz)) zt=height(nz)-100.*eps


end subroutine advance

