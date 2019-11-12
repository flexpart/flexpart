real function photo_O1D(sza)

  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    Nov 2014                                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !    INPUT:                                                                  *
  !    sza        solar zenith angle (degrees)                                 *
  !                                                                            *
  !    OUTPUT:                                                                 *
  !    photo_O1D  J(O1D) photoylsis rate                                       *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: iz,ik
  real :: sza
  real :: z1,z2,zg,f1,f2,dummy
  real :: photo_NO2
  integer, parameter :: nzenith=11
  real, parameter :: pi=3.1415927
  real, dimension(nzenith) :: zangle,fact_photo

  ! zangle: zenith angles for which fact_photo is tabulated
  ! fact_photo: conversion of photolysis rate of NO2 to photolysis 
  !     rate of O3 into O1D as a function of solar zenith angle

  zangle=(/0.,10.,20.,30.,40.,50.,60.,70.,78.,86.,90.0001/)
  fact_photo=(/0.4616E-02,0.4478E-02,0.4131E-02,0.3583E-02,0.2867E-02,&
    &0.2081E-02,0.1235E-02,0.5392E-03,0.2200E-03,0.1302E-03,0.0902E-03/)

  if (sza.lt.90.) then
    do iz=1,nzenith-1
      if(sza.ge.zangle(iz)) ik=iz
    end do
    z1=1./cos(zangle(ik)*pi/180.)
    z2=1./cos(zangle(ik+1)*pi/180.)
    zg=1./cos(sza*pi/180.)
    dummy=(zg-z1)/(z2-z1)
    f1=alog(fact_photo(ik))
    f2=alog(fact_photo(ik+1))
    photo_NO2=1.45e-2*exp(-0.4/cos(sza*pi/180.))
    photo_O1D=photo_NO2*exp(f1+(f2-f1)*dummy)
  else
    photo_O1D=0.
  endif

  return

end function photo_O1D


