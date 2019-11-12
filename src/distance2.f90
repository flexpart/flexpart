!-----------------------------------------------------------------------
function distance2(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH. COORDINATES ARE GIVEN IN RADIANS!
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN RADIANS
  !rlon1    - REAL LONGITUDE OF POINT 1 IN RADIANS
  !rlat2    - REAL LATITUDE OF POINT 2 IN RADIANS
  !rlon2    - REAL LONGITUDE OF POINT 2 IN RADIANS
  !
  !   OUTPUT ARGUMENT LIST:
  !distance2   - REAL GREAT CIRCLE DISTANCE IN KM
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real                    :: rlat1,rlon1,rlat2,rlon2,distance2
  real(kind=dp)           :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp

  if ((abs(rlat1-rlat2).lt.0.0003).and. &
       (abs(rlon1-rlon2).lt.0.0003)) then
    distance2=0.0_dp
  else

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    clat1=cos(real(rlat1,kind=dp))
    slat1=sin(real(rlat1,kind=dp))
    clat2=cos(real(rlat2,kind=dp))
    slat2=sin(real(rlat2,kind=dp))
    cdlon=cos(real(rlon1-rlon2,kind=dp))
    crd=slat1*slat2+clat1*clat2*cdlon
    distance2=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance2
