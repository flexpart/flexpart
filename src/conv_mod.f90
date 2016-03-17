!*******************************************************************************
!   Include file for convection
!   This file contains a global common block used by convect
!   and other subroutines
!   Author: P. Ferstl
!
!   Feb 2001
!
!*******************************************************************************

module conv_mod

  use par_mod, only: nconvlevmax, na, nxmax, nymax, nxmaxn, nymaxn, maxnests

  implicit none

  !integer,parameter :: nconvlevmax = nuvzmax-1, &
  !                     na = nconvlevmax+1
  !these parameters are defined in par_mod now!

  real :: pconv(nconvlevmax),phconv(na),dpr(nconvlevmax)
  real :: pconv_hpa(nconvlevmax),phconv_hpa(na)

  real :: ft(nconvlevmax), fq(nconvlevmax)
  real :: fmass(nconvlevmax,nconvlevmax),sub(nconvlevmax)
  real :: fmassfrac(nconvlevmax,nconvlevmax)
  real :: cbaseflux(0:nxmax-1,0:nymax-1)
  real :: cbasefluxn(0:nxmaxn-1,0:nymaxn-1,maxnests)
  real :: tconv(na),qconv(na),qsconv(na)
  real :: psconv,tt2conv,td2conv

  integer :: nconvlev,nconvtop

end module conv_mod
