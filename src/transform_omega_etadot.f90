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
     

    subroutine transform_omega_etadot(uuh,vvh,wwh,n)
                                      !i   i  i/o i     

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!*                                                                             *
!*  calculate eta_dot*dp_eta = wwh-dp_dt-u*dp_dx-v*dp_dy on costant eta levels *
!*  put eta_dot*dp_eta in wwh                                                  *
!*  here we use ps to obatin p,                                                *
!*  other variables used and in common block are:                              *
!*  ps(ix,jy,1,n),akz(kz),bkz(kz),ps_tplus1_and_min1(ix,jy,indj-1 or indj+1)   *
!*  if methodw is not  1 routne is not used                                    * 
!*                                                                             *
!*                                                                             *
!*     Author:                                                                 *
!*     M. Cassiani  2016                                                       *
!*                                                                             *
!*                                                                             *
!*                                                                             *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    
    use par_mod
    use com_mod  
    
    implicit none
    
    integer :: ix,jy,kz,n,ixp,ixm,deltat
    real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
    real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
    real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
    real(kind=4) :: pressure(0:nxmax-1,0:nymax-1,nwzmax)
    real(kind=4) :: dpressuredt(0:nxmax-1,0:nymax-1,nwzmax)
    real(kind=4) :: etadotdpdeta(0:nxmax-1,0:nymax-1,nwzmax) !not necessary but introduced for clarity eventually we coul rewrite wwh
    real(kind=4) :: oneoverdxinm(0:nymax-1)
    real(kind=4) :: oneoverdyinm,dpdx,dpdy
    
    
    
     do jy=1,ny-2
       oneoverdxinm(jy)=dxconst/cos((real(jy)*dy+ylat0)*pi180)
     end do
        oneoverdyinm=dyconst 
         
     
 
      do jy=0,nymin1
        do ix=0, nxfield-1           
          do kz=2,nuvz  !here nuvz is 27 (top of the domain) we skyp level one (extra kayer at the ground)
            pressure(ix,jy,kz)=akz(kz)+bkz(kz)*ps(ix,jy,1,n) 
            dpressuredt(ix,jy,kz)=bkz(kz)*(ps_tplus1_and_min1(ix,jy,2)-ps_tplus1_and_min1(ix,jy,1)) !ps_tplus1_and_min1 is already divided for the time interval
          end do
        end do
      end do
      
      
      do jy=1,nymin1-1 !skip the poles
        do ix=0, nxfield-1  
          ixp=ix+1
          ixm=ix-1
          if (ixm.eq.-1)  ixm=nxfield-1 !special treatment for 0 and nxfield-1,
          if (ixp.eq.nxfield) ixp=0    !special treatment for 0 and nxfield-1,
          do kz=2,nuvz  !here nuvz is 27 (top of the domain) we skip level one (extra layer at the ground)
            dpdx=(pressure(ixp,jy,kz)-pressure(ixm,jy,kz))*oneoverdxinm(jy)*0.5 !ix increase eastward  i.e. in the direction positive u (positive towards the east)
            dpdy=(pressure(ix,jy+1,kz)-pressure(ix,jy-1,kz))*oneoverdyinm*0.5  !jy incrrease northward i.e. in the direction positive v (positivwe towards the north)
            etadotdpdeta(ix,jy,kz)=wwh(ix,jy,kz)-dpressuredt(ix,jy,kz)-uuh(ix,jy,kz)*dpdx-vvh(ix,jy,kz)*dpdy
            ixp=ixp
          end do
        end do
      end do
      
      do jy=1,nymin1-1 !skip the poles
        do ix=0, nxfield-1  !special treatment for 0 and nxfield-1,          
          do kz=2,nuvz  !here nuvz is 27 (top of the domain) we skip level one (extra layer at the ground)
            wwh(ix,jy,kz)=etadotdpdeta(ix,jy,kz)          
          end do
        end do
      end do
       
      wwh(:,0,:)=0.
      wwh(:,nymin1,:)=0.
       
      do ix=0, nxfield-1  !         
        do kz=2,nuvz  !here nuvz is 27 (top of the domain) we skip level one (extra layer at the ground)
          wwh(ix,0,kz)=wwh(ix,0,kz)+etadotdpdeta(ix,1,kz)/nxfield  !set to averaged value of next circle 
          wwh(ix,nymin1,kz)= wwh(ix,nymin1,kz)+etadotdpdeta(ix,nymin1-1,kz)/nxfield  !set to averaged value of previous circle 
        end do
      end do
      
       
       
      return
      end
