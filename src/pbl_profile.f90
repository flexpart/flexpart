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

subroutine pbl_profile(ps,td2m,zml1,t2m,tml1,u10m,uml1,stress,hf)

  !********************************************************************
  !                                                                   *
  !                    G. WOTAWA, 1995-07-07                          *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! DESCRIPTION: CALCULATION OF FRICTION VELOCITY AND SURFACE SENS-   *
  !              IBLE HEAT FLUX USING THE PROFILE METHOD (BERKOVICZ   *
  !              AND PRAHM, 1982)                                     *
  !                                                                   *
  ! Output now is surface stress instead of ustar                     *
  !                                                                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! INPUT:                                                            *
  !                                                                   *
  !                                                                   *
  ! ps      surface pressure(Pa)                                      *
  ! td2m    two metre dew point(K)                                    *
  ! zml1    heigth of first model level (m)                           *
  ! t2m     two metre temperature (K)                                 *
  ! tml1    temperature first model level (K)                         *
  ! u10m    ten metre wind speed (ms-1)                               *
  ! uml1    wind speed first model level (ms-1)                       *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! OUTPUT:                                                           *
  !                                                                   *
  ! stress  surface stress (i.e., friction velocity (ms-1) squared    *
  !                         multiplied with air density)              *
  ! hf      surface sensible heat flux (Wm-2)                         *
  !                                                                   *
  !********************************************************************
  ! ustar   friction velocity (ms-1)                                  *
  ! maxiter maximum number of iterations                              *
  !********************************************************************

  use par_mod

  implicit none

  integer :: iter
  real :: ps,td2m,rhoa,zml1,t2m,tml1,u10m,uml1,ustar,hf
  real :: al,alold,aldiff,tmean,crit
  real :: deltau,deltat,thetastar,psim,psih,e,ew,tv,stress
  integer,parameter :: maxiter=10
  real,parameter    :: r1=0.74

  e=ew(td2m)               ! vapor pressure
  tv=t2m*(1.+0.378*e/ps)   ! virtual temperature
  rhoa=ps/(r_air*tv)       ! air density

  deltau=uml1-u10m         !! Wind Speed difference between
                           !! Model level 1 and 10 m

  if(deltau.le.0.001) then    !! Monin-Obukhov Theory not
    al=9999.               !! applicable --> Set dummy values
    ustar=0.01
    stress=ustar*ustar*rhoa
    hf=0.0
    return
  endif
  deltat=tml1-t2m+0.0098*(zml1-2.)  !! Potential temperature difference
                                    !! between model level 1 and 10 m

  if(abs(deltat).le.0.03) then    !! Neutral conditions
    hf=0.0
    al=9999.
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    stress=ustar*ustar*rhoa
    return
  endif

  tmean=0.5*(t2m+tml1)
  crit=(0.0219*tmean*(zml1-2.0)*deltau**2)/ &
       (deltat*(zml1-10.0)**2)
  if((deltat.gt.0).and.(crit.le.1.)) then
                                    !! Successive approximation will
    al=50.                          !! not converge
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1/2.)-psih(zml1,al)+psih(2.,al))
    hf=rhoa*cpa*ustar*thetastar
    stress=ustar*ustar*rhoa
    return
  endif

  al=9999.                 ! Start iteration assuming neutral conditions
  do iter=1,maxiter
    alold=al
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1/2.)-psih(zml1,al)+psih(2.,al))
    al=(tmean*ustar**2)/(ga*vonkarman*thetastar)
    aldiff=abs((al-alold)/alold)
    if(aldiff.lt.0.01) goto 30  !! Successive approximation successful
  end do
30   hf=rhoa*cpa*ustar*thetastar
  if(al.gt.9999.) al=9999.
  if(al.lt.-9999.) al=-9999.

  stress=ustar*ustar*rhoa

end subroutine pbl_profile
