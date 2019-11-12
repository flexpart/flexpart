subroutine hanna_short(z)
  !                       i
  !*****************************************************************************
  !                                                                            *
  !   Computation of \sigma_i and \tau_L based on the scheme of Hanna (1982)   *
  !                                                                            *
  !   Author: A. Stohl                                                         *
  !                                                                            *
  !   4 December 1997                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dsigwdz [1/s]     vertical gradient of sigw                                *
  ! ol [m]            Obukhov length                                           *
  ! sigu, sigv, sigw  standard deviations of turbulent velocity fluctuations   *
  ! tlu [s]           Lagrangian time scale for the along wind component.      *
  ! tlv [s]           Lagrangian time scale for the cross wind component.      *
  ! tlw [s]           Lagrangian time scale for the vertical wind component.   *
  ! ust, ustar [m/s]  friction velocity                                        *
  ! wst, wstar [m/s]  convective velocity scale                                *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use hanna_mod

  implicit none

  real :: z



  !**********************
  ! 1. Neutral conditions
  !**********************

  if (h/abs(ol).lt.1.) then
    ust=max(1.e-4,ust)
    sigw=1.3*exp(-2.e-4*z/ust)
    dsigwdz=-2.e-4*sigw
    sigw=sigw*ust+1.e-2
    tlw=0.5*z/sigw/(1.+1.5e-3*z/ust)


  !***********************
  ! 2. Unstable conditions
  !***********************

  else if (ol.lt.0.) then


  ! Determine sigmas
  !*****************

    sigw=sqrt(1.2*wst**2*(1.-.9*zeta)*zeta**0.66666+ &
         (1.8-1.4*zeta)*ust**2)+1.e-2
    dsigwdz=0.5/sigw/h*(-1.4*ust**2+wst**2* &
         (0.8*max(zeta,1.e-3)**(-.33333)-1.8*zeta**0.66666))


  ! Determine average Lagrangian time scale
  !****************************************

    if (z.lt.abs(ol)) then
      tlw=0.1*z/(sigw*(0.55-0.38*abs(z/ol)))
    else if (zeta.lt.0.1) then
      tlw=0.59*z/sigw
    else
      tlw=0.15*h/sigw*(1.-exp(-5*zeta))
    endif


  !*********************
  ! 3. Stable conditions
  !*********************

  else
    sigw=1.e-2+1.3*ust*(1.-zeta)
    dsigwdz=-1.3*ust/h
    tlw=0.1*h/sigw*zeta**0.8
  endif


  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)
  if (dsigwdz.eq.0.) dsigwdz=1.e-10

end subroutine hanna_short
