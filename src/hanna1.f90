! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine hanna1(z)
  !                  i
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

  real :: z,s1,s2



  !**********************
  ! 1. Neutral conditions
  !**********************

  if (h/abs(ol).lt.1.) then

    ust=max(1.e-4,ust)
    sigu=2.0*ust*exp(-3.e-4*z/ust)
    sigu=max(sigu,1.e-5)
    sigv=1.3*ust*exp(-2.e-4*z/ust)
    sigv=max(sigv,1.e-5)
    sigw=sigv
    dsigw2dz=-6.76e-4*ust*exp(-4.e-4*z/ust)
    tlu=0.5*z/sigw/(1.+1.5e-3*z/ust)
    tlv=tlu
    tlw=tlu


  !***********************
  ! 2. Unstable conditions
  !***********************

  else if (ol.lt.0.) then


  ! Determine sigmas
  !*****************

    sigu=ust*(12-0.5*h/ol)**0.33333
    sigu=max(sigu,1.e-6)
    sigv=sigu

    if (zeta.lt.0.03) then
      sigw=0.96*wst*(3*zeta-ol/h)**0.33333
      dsigw2dz=1.8432*wst*wst/h*(3*zeta-ol/h)**(-0.33333)
    else if (zeta.lt.0.4) then
      s1=0.96*(3*zeta-ol/h)**0.33333
      s2=0.763*zeta**0.175
      if (s1.lt.s2) then
        sigw=wst*s1
        dsigw2dz=1.8432*wst*wst/h*(3*zeta-ol/h)**(-0.33333)
      else
        sigw=wst*s2
        dsigw2dz=0.203759*wst*wst/h*zeta**(-0.65)
      endif
    else if (zeta.lt.0.96) then
      sigw=0.722*wst*(1-zeta)**0.207
      dsigw2dz=-.215812*wst*wst/h*(1-zeta)**(-0.586)
    else if (zeta.lt.1.00) then
      sigw=0.37*wst
      dsigw2dz=0.
    endif
    sigw=max(sigw,1.e-6)


  ! Determine average Lagrangian time scale
  !****************************************

    tlu=0.15*h/sigu
    tlv=tlu
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
    sigu=2.*ust*(1.-zeta)
    sigv=1.3*ust*(1.-zeta)
    sigu=max(sigu,1.e-6)
    sigv=max(sigv,1.e-6)
    sigw=sigv
    dsigw2dz=3.38*ust*ust*(zeta-1.)/h
    tlu=0.15*h/sigu*(sqrt(zeta))
    tlv=0.467*tlu
    tlw=0.1*h/sigw*zeta**0.8
  endif




  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)


end subroutine hanna1
