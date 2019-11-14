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

real function zenithangle(ylat,xlon,jul)

  !*********************************************************************
  !                                                                    *
  !                      Author: G. WOTAWA                             *
  !                      Date: 1993-11-17                              *
  !                      Project: POP-M                                *
  !                      Last update:                                  *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !     DESCRIPTION: This function returns the sinus of solar          *
  !                  elevation as a function of geographic longitude,  *
  !                  latitude and GMT-Time.                            *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !     INPUT:                                                         *
  !                                                                    *
  !            ylat          geographical latitude  [DEG]              *
  !            xlon          geographical longitude [DEG]              *
  !            jjjj          Year                                      *
  !            mm            Month                                     *
  !            dd            Day                                       *
  !            hh            Hour                                      *
  !            minute        Minute                                    *
  !                                                                    *
  !*********************************************************************

  use par_mod, only: dp

  implicit none

  integer :: jjjj,mm,id,iu,minute,yyyymmdd,hhmmss
  integer :: ndaynum
  real :: sinsol,solelev,ylat,xlon
  real :: rnum,rylat,ttime,dekl,rdekl,eq
  real,parameter :: pi=3.1415927
  real(kind=dp)  :: jul

  call caldate(jul,yyyymmdd,hhmmss)
  jjjj=yyyymmdd/10000
  mm=yyyymmdd/100-jjjj*100
  id=yyyymmdd-jjjj*10000-mm*100
  iu=hhmmss/10000
  minute=hhmmss/100-100*iu

  ndaynum=31*(mm-1)+id
  if(mm.gt.2) ndaynum=ndaynum-int(0.4*mm+2.3)
  if((mm.gt.2).and.(jjjj/4*4.eq.jjjj)) ndaynum=ndaynum+1

  rnum=2.*pi*ndaynum/365.
  rylat=pi*ylat/180.
  ttime=real(iu)+real(minute)/60.

  dekl=0.396+3.631*sin(rnum)+0.038*sin(2.*rnum)+0.077*sin(3.*rnum)- &
       22.97*cos(rnum)-0.389*cos(2.*rnum)-0.158*cos(3.*rnum)
  rdekl=pi*dekl/180.

  eq=(0.003-7.343*sin(rnum)-9.47*sin(2.*rnum)- &
       0.329*sin(3.*rnum)-0.196*sin(4.*rnum)+ &
       0.552*cos(rnum)-3.020*cos(2.*rnum)- &
       0.076*cos(3.*rnum)-0.125*cos(4.*rnum))/60.

  sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)* &
       cos((ttime-12.+xlon/15.+eq)*pi/12.)
  ! Calculate the maximum solar elevation on that day
  !sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)*
  !    &       cos((eq)*pi/12.)
  solelev=asin(sinsol)*180./pi
  zenithangle=90.-solelev

  return
end function zenithangle
