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

subroutine caldate(juldate,yyyymmdd,hhmiss)
  !                      i       o       o
  !*****************************************************************************
  !                                                                            *
  !     Calculates the Gregorian date from the Julian date                     *
  !                                                                            *
  !     AUTHOR: Andreas Stohl (21 January 1994), adapted from Numerical Recipes*
  !                                                                            *
  !     Variables:                                                             *
  !     dd             Day                                                     *
  !     hh             Hour                                                    *
  !     hhmiss         Hour, Minute, Second                                    *
  !     ja,jb,jc,jd,je help variables                                          *
  !     jalpha         help variable                                           *
  !     juldate        Julian Date                                             *
  !     julday         help variable                                           *
  !     mi             Minute                                                  *
  !     mm             Month                                                   *
  !     ss             Seconds                                                 *
  !     yyyy           Year                                                    *
  !     yyyymmdd       Year, Month, Day                                        *
  !                                                                            *
  !     Constants:                                                             *
  !     igreg          help constant                                           *
  !                                                                            *
  !*****************************************************************************

  use par_mod, only: dp

  implicit none

  integer           :: yyyymmdd,yyyy,mm,dd,hhmiss,hh,mi,ss
  integer           :: julday,ja,jb,jc,jd,je,jalpha
  real(kind=dp)     :: juldate
  integer,parameter :: igreg=2299161

  julday=int(juldate)
  if(julday.ge.igreg)then
    jalpha=int(((julday-1867216)-0.25)/36524.25)
    ja=julday+1+jalpha-int(0.25*jalpha)
  else
    ja=julday
  endif
  jb=ja+1524
  jc=int(6680.+((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)
  dd=jb-jd-int(30.6001*je)
  mm=je-1
  if (mm.gt.12) mm=mm-12
  yyyy=jc-4715
  if (mm.gt.2) yyyy=yyyy-1
  if (yyyy.le.0) yyyy=yyyy-1

  yyyymmdd=10000*yyyy+100*mm+dd
  hh=int(24._dp*(juldate-real(julday,kind=dp)))
  mi=int(1440._dp*(juldate-real(julday,kind=dp))-60._dp*real(hh,kind=dp))
  ss=nint(86400._dp*(juldate-real(julday,kind=dp))-3600._dp*real(hh,kind=dp)- &
       60._dp*real(mi,kind=dp))
  if (ss.eq.60) then  ! 60 seconds = 1 minute
    ss=0
    mi=mi+1
  endif
  if (mi.eq.60) then
    mi=0
    hh=hh+1
  endif
  hhmiss=10000*hh+100*mi+ss

end subroutine caldate
