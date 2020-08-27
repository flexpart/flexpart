! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!##################################################################
!##################################################################
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

function f_qvsat( p, t )

  !PURPOSE:
  !
  !Calculate the saturation specific humidity using enhanced Teten's
  !formula.
  !
  !AUTHOR: Yuhe Liu
  !01/08/1998
  !
  !MODIFICATION HISTORY:
  !
  !INPUT :
  !  p        Pressure (Pascal)
  !  t        Temperature (K)
  !OUTPUT:
  !  f_qvsat  Saturation water vapor specific humidity (kg/kg).
  !
  !Variable Declarations.
  !

  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: f_qvsat   ! Saturation water vapor specific humidity (kg/kg)
  real :: f_esl,f_esi,fespt

  real,parameter ::  rd     = 287.0 ! Gas constant for dry air  (m**2/(s**2*K))
  real,parameter ::  rv     = 461.0 ! Gas constant for water vapor  (m**2/(s**2*K)).
  real,parameter ::  rddrv  = rd/rv


  ! Change by A. Stohl to save computation time:
  ! IF ( t.ge.273.15 ) THEN     ! for water
  if ( t.ge.253.15 ) then      ! modification Petra Seibert
                               ! (supercooled water may be present)
    fespt=f_esl(p,t)
  else
    fespt=f_esi(p,t)
  endif

!!$  f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)      !old
  if (p-(1.0-rddrv)*fespt == 0.) then                     !bugfix
     f_qvsat = 1.
  else
     f_qvsat = rddrv * fespt / (p-(1.0-rddrv)*fespt)
  end if

  return
end function f_qvsat


function f_esl( p, t )

  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: f_esl     ! Saturation water vapor pressure over liquid water

  real :: f

  !#######################################################################
  !
  !Saturation specific humidity parameters used in enhanced Teten's
  !formula. (See A. Buck, JAM 1981)
  !
  !#######################################################################

  real,parameter ::  satfwa = 1.0007
  real,parameter ::  satfwb = 3.46e-8  ! for p in Pa

  real,parameter ::  satewa = 611.21   ! es in Pa
  real,parameter ::  satewb = 17.502
  real,parameter ::  satewc = 32.18

  real,parameter ::  satfia = 1.0003
  real,parameter ::  satfib = 4.18e-8  ! for p in Pa

  real,parameter ::  sateia = 611.15   ! es in Pa
  real,parameter ::  sateib = 22.452
  real,parameter ::  sateic = 0.6

  f = satfwa + satfwb * p
  f_esl = f * satewa * exp( satewb*(t-273.15)/(t-satewc) )

  return
end function f_esl

function f_esi( p, t )

  implicit none

  real :: p         ! Pressure (Pascal)
  real :: t         ! Temperature (K)
  real :: f_esi     ! Saturation water vapor pressure over ice (Pa)

  real :: f

  !#######################################################################
  !
  !Saturation specific humidity parameters used in enhanced Teten's
  !formula. (See A. Buck, JAM 1981)
  !
  !#######################################################################
  !
  real,parameter ::  satfwa = 1.0007
  real,parameter ::  satfwb = 3.46e-8  ! for p in Pa

  real,parameter ::  satewa = 611.21   ! es in Pa
  real,parameter ::  satewb = 17.502
  real,parameter ::  satewc = 32.18

  real,parameter ::  satfia = 1.0003
  real,parameter ::  satfib = 4.18e-8  ! for p in Pa

  real,parameter ::  sateia = 611.15   ! es in Pa
  real,parameter ::  sateib = 22.452
  real,parameter ::  sateic = 0.6

  f = satfia + satfib * p
  f_esi = f * sateia * exp( sateib*(t-273.15)/(t-sateic) )

  return
end function f_esi
