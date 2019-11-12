!-----------------------------------------------------------------------
function distance(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH.
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN DEGREES
  !rlon1    - REAL LONGITUDE OF POINT 1 IN DEGREES
  !rlat2    - REAL LATITUDE OF POINT 2 IN DEGREES
  !rlon2    - REAL LONGITUDE OF POINT 2 IN DEGREES
  !
  !   OUTPUT ARGUMENT LIST:
  !distance - REAL GREAT CIRCLE DISTANCE IN KILOMETERS
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real          :: rlat1,rlon1,rlat2,rlon2,distance
  real(kind=dp) :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp, dpr=180.0_dp/pi
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ((abs(rlat1-rlat2).lt.0.03).and. &
       (abs(rlon1-rlon2).lt.0.03)) then
    distance=0.
  else
    clat1=cos(real(rlat1,kind=dp)/dpr)
    slat1=sin(real(rlat1,kind=dp)/dpr)
    clat2=cos(real(rlat2,kind=dp)/dpr)
    slat2=sin(real(rlat2,kind=dp)/dpr)
    cdlon=cos(real((rlon1-rlon2),kind=dp)/dpr)
    crd=slat1*slat2+clat1*clat2*cdlon
    distance=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance
