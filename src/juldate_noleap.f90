!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann, Ignacio Pisso                                      *
!                                                                     *
! This file is part of FLEXPART-NorESM                                *
!                                                                     *
! FLEXPART-NorESM is free software: you can redistribute it           *
! and/or modify                                                       *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART-NorESM is distributed in the hope that it will be useful,  *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART-NorESM.                                         *
!  If not, see <http://www.gnu.org/licenses/>.                        * 
!**********************************************************************

      DOUBLE PRECISION FUNCTION juldate(YYYYMMDD,HHMISS)
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*                                                                             *
!*     Calculates days from a starting date  using a no leap calendar          *
!*     NOTE  this is not really a julian date but we kept the name of          *
!*     the original FLEXPART routine, see commennts below if                   *
!*      one needs it to be a real julian date routine                          *
!*     NOTE: read comments below if one wants to include leap years            *
!*     NOTE: read comments below if one wants to a real juldate  with          *
!*           gregorian calndar
!*                                                                             *
!*     AUTHOR: Massimo Cassiani 2016                                           *
!*                                                                             *
!*     Variables:                                                              *
!*     DD             Day                                                      *
!*     HH             Hour                                                     *
!*     HHMISS         Hour, minute + second                                    *
!*     JA,JM,JY       help variables                                           *
!*     juldate        "julian date" from arbitrary starting date yearzer       *
!*     JULDAY         help variable                                            *
!*     MI             Minute                                                   *
!*     MM             Month                                                    *
!*     SS             Second                                                   *
!*     YYYY           Year                                                     *
!*     YYYYMMDDHH     Date and Time                                            *
!*                                                                             *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      use par_mod
      IMPLICIT NONE
     
      
      INTEGER YYYYMMDD,YYYY,MM,DD,HH,MI,SS,HHMISS
      INTEGER JULDAY,JY,JM,JA,IGREG,YYYY1,MM1
      INTEGER AP,MP,DPM,YAP
 !     DOUBLE PRECISION juldate
!      PARAMETER (IGREG=15+31*(10+12*1582))
      
           
      YYYY=YYYYMMDD/10000
      !YYYY1=YYYYMMDD/10000
      MM=(YYYYMMDD-10000*YYYY)/100
      !MM1=(YYYYMMDD-10000*YYYY)/100
      DD=YYYYMMDD-10000*YYYY-100*MM
      HH=HHMISS/10000
      MI=(HHMISS-10000*HH)/100
      SS=HHMISS-10000*HH-100*MI
      IF (YYYY.EQ.0) PAUSE &
     'There is Year ZERO? then modify JULDATE_NOLEAP routine!!!'
      AP= (14- MM)/12    ! will result in a 1 for January (month 1) and February (month 2).  The result is 0 for the other 10 months.
      MP=MM +12*AP-3     ! results in a 10 for January, 11 for February, 0 for March, 1 for April, ..., and a 9 for December
      DPM=(153*MP+2)/5   ! will give the number of days in the previous months starting from month zero
      YAP= YYYY - yearorigin - AP  !adds yearorigin to the year so that we will start counting years  from the year (yearorigin=-YYYY)
                                   !note that year 1 corresponds to 1 CE, year 0 corresponds to 1 BCE, year 1 corresponds to 2 BCE, and so onThere is no year between 1 CE and 1 BCE
      JULDAY=DD+DPM+365*YAP        !+ YAP/4-YAP/100+YAP/400 !add this part if you want to count leap year:                                    
                                   !-32405 !!NOTE, subtract this further number and put origin to 4800 BCE if you want true julian date starting at 4713 BCE (24 nov 4714BC). The -32045 comes from (4800-4712)*365+22 (leap years) - (6+31)(nov+dec) - (31+29) (jan+feb_leap)
      JULDAY=JULDAY+59-365        !-59 for correcting from starting 1 of march see MP. -365 because we consider that there is no zero year.
                                  !NOte put this correction to zero i.e. JULDAY=JULDAY and look comments just above if you want true julian date from gregorian calendar with no zero.
            
      
      juldate=DBLE((JULDAY))+DBLE(DBLE(HH)/24.)+ &
      DBLE(DBLE(MI)/1440.)+DBLE(DBLE(SS)/86400.)

      RETURN
      END
