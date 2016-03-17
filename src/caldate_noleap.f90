!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann                                                     *
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
 

 SUBROUTINE CALDATE(JULDATE,YYYYMMDD,HHMISS)
!c                     i       o       o
!
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!*                                                                             *
!*     Calculates a pseudo-Gregorian date (from year zero) using  the          *
!*     juldate varaiable                                                       *
!*     juldate actually contains days from arbitrary origin date (yearorigin)  *
!*     using no leap calendar                                                  *
!*     Author:                                                                 *
!*     M. Cassiani  2016                                                       *
!*                                                                             *
!*                                                                             *
!                                                                              *
!*                                                                             *
!*                                                                             *
!      Variables:                                                              *
!*     DD             Day                                                      *
!*     HH             Hour                                                     *
!*     HHMISS         Hour, Minute, Second                                     *
!*     JA,JB,JC,JD,JE help variables                                           *
!*     JALPHA         help variable                                            *
!*     JULDATE        Julian Date                                              *
!*     JULDAY         help variable                                            *
!*     MI             Minute                                                   *
!*     MM             Month                                                    *
!*     SS             Seconds                                                  *
!*     YYYY           Year                                                     *
!*     YYYYMMDD       Year, Month, Day                                         *
!*                                                                             *
!*     Constants:                                                              *
!*     IGREG          help constant                                            *
!*                                                                             *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      use par_mod
      IMPLICIT NONE

      INTEGER YYYYMMDD,YYYY,MM,DD,HHMISS,HH,MI,SS
      INTEGER JULDAY,JA,JB,JC,JD,JE,IGREG,JALPHA,DAYS,ddpm
      DOUBLE PRECISION JULDATE
      INTEGER daysannozero
      
            
      YYYY=INT(JULDATE)/365   
      
      
      JULDAY=INT(JULDATE)      
      DAYS=JULDAY-YYYY*365
      if (DAYS.gt.59) then
          MM=int((((DAYS-0.1-59.)*5.)-2.)/153.)
          ddpm=((153*MM)+2)/5
          DD=DAYS-ddpm-59
      else 
          MM=int((((DAYS-0.1+306.)*5.)-2.)/153.)
          ddpm=((153*MM)+2)/5
          DD=DAYS+306-ddpm
      end if
      
      if (MM.le.9) then
      MM=MM+3
      else
      MM=MM-9
      end if
      
      YYYY=YYYY+yearorigin
      if (YYYY.ge.0.and.DAYS.ne.0) YYYY=YYYY+1
      
      
      YYYYMMDD=10000*YYYY+100*MM+DD
      
      HH=INT(24.*(JULDATE-FLOAT(JULDAY)))
      MI=INT(1440.*(JULDATE-FLOAT(JULDAY))-60.*FLOAT(HH))
      SS=NINT(86400.*(JULDATE-FLOAT(JULDAY))-3600.*FLOAT(HH)) &
     -60.*FLOAT(MI)
      IF (SS.EQ.60) THEN  ! 60 seconds = 1 minute
        SS=0
        MI=MI+1
      ENDIF
      IF (MI.EQ.60) THEN
        MI=0
        HH=HH+1
      ENDIF
      HHMISS=10000*HH+100*MI+SS

      END
