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


Function dewpoint(TK,Q,P)
!*****************************************************************************
!                                                                            *
!     This routine calculates dew point in kelvin                             *
!     files from NorESM                                                      *
!                                                                            *
!     Author:                                                                *
!     M. Cassiani  2016                                                      *
!                                                                            *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! TK temperatur in kelvin                                                    *
! Q specific humidity                                                        *
! P Surface pressure                                                                 *
!*****************************************************************************


    implicit none
    real(kind=4) :: TK,Q,P,es,t,e,dewpoint,a,b,c,dewpoint2
    
    t=TK-273.15 !from kelvin to celsius
    if (t.gt.0) then !from Campbell and Norman introduction to environmental BioPhysics 1998 p. 41-43
                     !considering t<0 above ice t>0 above water
       a=611 !Pa
       b=17.502 !
       c=240.97 !
    else 
       a=611 !Pa
       b=21.87 !
       c=265.5 !  
    end if
    !es=a*exp(b*t/(t+c)) ! from Campbell and Norman, , springer verlag.
    !es=exp(53.67957-6473.7/(t+273.16)-4.8451*log(t+273.16))  !saturation pressure
    e=Q*P/0.622 !approximate partial pressure of water wapor in Pa, 0.622 is the ratio of molar mass of water and air
                !also (1-0.622)e was considered negligible compared to P.
    !u=e/es !relative humidity
    dewpoint=(c*log(e/a)/(b-log(e/a)))+273.15  !from Campbell and Norman, see above
    !dewpoint2=5.42*10**3/(log(2.53*10.**11*0.622/(Q*P))) !-273 !Rogers and Yau's "A Short Course in Cloud Physics" 
    
end
