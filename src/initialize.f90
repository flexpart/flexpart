subroutine initialize(itime,ldt,up,vp,wp, &
       usigold,vsigold,wsigold,xt,yt,zt,icbt)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of trajectories utilizing a zero-acceleration scheme. The time*
  !  step is determined by the Courant-Friedrichs-Lewy (CFL) criterion. This   *
  !  means that the time step must be so small that the displacement within    *
  !  this time step is smaller than 1 grid distance. Additionally, a temporal  *
  !  CFL criterion is introduced: the time step must be smaller than the time  *
  !  interval of the wind fields used for interpolation.                       *
  !  For random walk simulations, these are the only time step criteria.       *
  !  For the other options, the time step is also limited by the Lagrangian    *
  !  time scale.                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Literature:                                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          current temporal position                               *
  ! ldt [s]            Suggested time step for next integration                *
  ! ladvance [s]       Total integration time period                           *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! up,vp,wp           random velocities due to turbulence                     *
  ! usig,vsig,wsig     uncertainties of wind velocities due to interpolation   *
  ! usigold,vsigold,wsigold  like usig, etc., but for the last time step       *
  ! xt,yt,zt           Next time step's spatial position of trajectory         *
  !                                                                            *
  !                                                                            *
  ! Constants:                                                                 *
  ! cfl                factor, by which the time step has to be smaller than   *
  !                    the spatial CFL-criterion                               *
  ! cflt               factor, by which the time step has to be smaller than   *
  !                    the temporal CFL-criterion                              *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod
  use random_mod, only: ran3

  implicit none

  integer :: itime
  integer :: ldt,nrand
  integer(kind=2) :: icbt
  real :: zt,dz,dz1,dz2,up,vp,wp,usigold,vsigold,wsigold
  real(kind=dp) :: xt,yt
  save idummy

  integer :: idummy = -7

  icbt=1           ! initialize particle to "no reflection"

  nrand=int(ran3(idummy)*real(maxrand-1))+1


  !******************************
  ! 2. Interpolate necessary data
  !******************************

  ! Compute maximum mixing height around particle position
  !*******************************************************

  ix=int(xt)
  jy=int(yt)
  ixp=ix+1
  jyp=jy+1

  h=max(hmix(ix ,jy,1,memind(1)), &
       hmix(ixp,jy ,1,memind(1)), &
       hmix(ix ,jyp,1,memind(1)), &
       hmix(ixp,jyp,1,memind(1)), &
       hmix(ix ,jy ,1,memind(2)), &
       hmix(ixp,jy ,1,memind(2)), &
       hmix(ix ,jyp,1,memind(2)), &
       hmix(ixp,jyp,1,memind(2)))

  zeta=zt/h


  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************

  if (zeta.le.1.) then

    call interpol_all(itime,real(xt),real(yt),zt)


  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************

    dz1=zt-height(indz)
    dz2=height(indzp)-zt
    dz=1./(dz1+dz2)

    u=(dz1*uprof(indzp)+dz2*uprof(indz))*dz
    v=(dz1*vprof(indzp)+dz2*vprof(indz))*dz
    w=(dz1*wprof(indzp)+dz2*wprof(indz))*dz


  ! Compute the turbulent disturbances

  ! Determine the sigmas and the timescales
  !****************************************

    if (turbswitch) then
      call hanna(zt)
    else
      call hanna1(zt)
    endif


  ! Determine the new diffusivity velocities
  !*****************************************

    if (nrand+2.gt.maxrand) nrand=1
    up=rannumb(nrand)*sigu
    vp=rannumb(nrand+1)*sigv
    wp=rannumb(nrand+2)
    if (.not.turbswitch) then     ! modified by mc
      wp=wp*sigw
    else if (cblflag.eq.1) then   ! modified by mc
      if(-h/ol.gt.5) then
!if (ol.lt.0.) then
!if (ol.gt.0.) then !by mc : only for test correct is lt.0
        call initialize_cbl_vel(idummy,zt,ust,wst,h,sigw,wp,ol)
      else
        wp=wp*sigw
      end if
    end if


  ! Determine time step for next integration
  !*****************************************

    if (turbswitch) then
      ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5), &
           0.5/abs(dsigwdz),600.)*ctl)
    else
      ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5),600.)*ctl)
    endif
    ldt=max(ldt,mintime)


    usig=(usigprof(indzp)+usigprof(indz))/2.
    vsig=(vsigprof(indzp)+vsigprof(indz))/2.
    wsig=(wsigprof(indzp)+wsigprof(indz))/2.

  else



  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************


  ! Interpolate the wind
  !*********************

    call interpol_wind(itime,real(xt),real(yt),zt)


  ! Compute everything for above the PBL

  ! Assume constant turbulent perturbations
  !****************************************

    ldt=abs(lsynctime)

    if (nrand+1.gt.maxrand) nrand=1
    up=rannumb(nrand)*0.3
    vp=rannumb(nrand+1)*0.3
    nrand=nrand+2
    wp=0.
    sigw=0.

  endif

  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************


  ! It is assumed that the average interpolation error is 1/2 sigma
  ! of the surrounding points, autocorrelation time constant is
  ! 1/2 of time interval between wind fields
  !****************************************************************

  if (nrand+2.gt.maxrand) nrand=1
  usigold=rannumb(nrand)*usig
  vsigold=rannumb(nrand+1)*vsig
  wsigold=rannumb(nrand+2)*wsig

end subroutine initialize
