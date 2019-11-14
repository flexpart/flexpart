! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine getvdep(n,ix,jy,ust,temp,pa,L,gr,rh,rr,snow,vdepo)
  !                   i i  i   i   i   i  i i  i  i    i o
  !*****************************************************************************
  !                                                                            *
  !  This routine calculates the dry deposition velocities.                    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 December 1996                                                       *
  !     Sabine Eckhardt, Jan 07                                                *
  !     if the latitude is negative: add half a year to the julian day         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! gr [W/m2]         global radiation                                         *
  ! L [m]             Obukhov length                                           *
  ! nyl               kinematic viscosity                                      *
  ! pa [Pa]           surface air pressure                                     *
  ! ra [s/m]          aerodynamic resistance                                   *
  ! raquer [s/m]      average aerodynamic resistance                           *
  ! rh [0-1]          relative humidity                                        *
  ! rhoa              density of the air                                       *
  ! rr [mm/h]         precipitation rate                                       *
  ! temp [K]          2m temperature                                           *
  ! tc [C]            2m temperature                                           *
  ! ust [m/s]         friction velocity                                        *
  ! snow [m of water equivalent] snow depth                                    *
  ! xlanduse          fractions of numclasS landuses for each model grid point *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: yyyymmdd,hhmmss,yyyy,mmdd,n,lseason,i,j,ix,jy
  real :: vdepo(maxspec),vd,rb(maxspec),rc(maxspec),raquer,ylat
  real :: raerod,ra,ust,temp,tc,pa,L,gr,rh,rr,myl,nyl,rhoa,diffh2o,snow
  real :: slanduse(numclass)
  real,parameter :: eps=1.e-5
  real(kind=dp) :: jul

  ! Calculate month and determine the seasonal category
  !****************************************************

  jul=bdate+real(wftime(n),kind=dp)/86400._dp

  ylat=jy*dy+ylat0
  if (ylat.lt.0) then
      jul=jul+365/2
  endif


  call caldate(jul,yyyymmdd,hhmmss)
  yyyy=yyyymmdd/10000
  mmdd=yyyymmdd-10000*yyyy

  if ((ylat.gt.-20).and.(ylat.lt.20)) then
     mmdd=600 ! summer
  endif

  if ((mmdd.ge.1201).or.(mmdd.le.301)) then
    lseason=4
  else if ((mmdd.ge.1101).or.(mmdd.le.331)) then
    lseason=3
  else if ((mmdd.ge.401).and.(mmdd.le.515)) then
    lseason=5
  else if ((mmdd.ge.516).and.(mmdd.le.915)) then
    lseason=1
  else
    lseason=2
  endif

  ! Calculate diffusivity of water vapor
  !************************************
  diffh2o=2.11e-5*(temp/273.15)**1.94*(101325/pa)

  ! Conversion of temperature from K to C
  !**************************************

  tc=temp-273.15

  ! Calculate dynamic viscosity
  !****************************

  if (tc.lt.0) then
    myl=(1.718+0.0049*tc-1.2e-05*tc**2)*1.e-05
  else
    myl=(1.718+0.0049*tc)*1.e-05
  endif

  ! Calculate kinematic viscosity
  !******************************

  rhoa=pa/(287.*temp)
  nyl=myl/rhoa


  ! 0. Set all deposition velocities zero
  !**************************************

  do i=1,nspec
    vdepo(i)=0.
  end do


  ! 1. Compute surface layer resistances rb
  !****************************************

  call getrb(nspec,ust,nyl,diffh2o,reldiff,rb)

  ! change for snow
  do j=1,numclass
    if (snow.gt.0.001) then ! 10 mm
       if (j.eq.12) then
         slanduse(j)=1.
       else
         slanduse(j)=0.
       endif
    else
       slanduse(j)=xlanduse(ix,jy,j)
    endif
  end do

  raquer=0.
  do j=1,numclass            ! loop over all landuse classes

    if (slanduse(j).gt.eps)  then

  ! 2. Calculate aerodynamic resistance ra
  !***************************************

      ra=raerod(L,ust,z0(j))
      raquer=raquer+ra*slanduse(j)

  ! 3. Calculate surface resistance for gases
  !******************************************

      call getrc(nspec,lseason,j,tc,gr,rh,rr,rc)

  ! 4. Calculate deposition velocities for gases and ...
  ! 5. ... sum deposition velocities for all landuse classes
  !*********************************************************

      do i=1,nspec
        if (reldiff(i).gt.0.) then
          if ((ra+rb(i)+rc(i)).gt.0.) then
            vd=1./(ra+rb(i)+rc(i))
          else
            vd=9.999
          endif
          vdepo(i)=vdepo(i)+vd*slanduse(j)
        endif
      end do
    endif
  end do


  ! 6. Calculate deposition velocities for particles
  !*************************************************

  call partdep(nspec,density,fract,schmi,vset,raquer,ust,nyl,vdepo)

  !if (debug_mode) then
  !  print*,'getvdep:188: vdepo=', vdepo
    !stop
  !endif

  ! 7. If no detailed parameterization available, take constant deposition
  !    velocity if that is available
  !***********************************************************************

  do i=1,nspec
    if ((reldiff(i).lt.0.).and.(density(i).lt.0.).and. &
         (dryvel(i).gt.0.)) then
      vdepo(i)=dryvel(i)
    endif
  end do


end subroutine getvdep
