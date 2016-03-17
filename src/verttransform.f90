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


subroutine verttransform(n,uuh,vvh,wwh,pvh,itime,indj)
  !                      i   i  i   i   i   i      i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !*****************************************************************************
  !     Modified by M. Cassiani 2106 to trasnform etadotdpdeta using           *
  !     NorESM levels                                                          *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) from common block
  !*****************************************************************************
  ! Sabine Eckhardt, March 2007
  ! added the variable cloud for use with scavenging - descr. in com_mod
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! clouds(0:nxmax,0:nymax,0:nzmax,2) cloud field for wet deposition           *
  ! uu(0:nxmax,0:nymax,nzmax,2)     wind components in x-direction [m/s]       *
  ! vv(0:nxmax,0:nymax,nzmax,2)     wind components in y-direction [m/s]       *
  ! ww(0:nxmax,0:nymax,nzmax,2)     wind components in z-direction [deltaeta/s]*
  ! tt(0:nxmax,0:nymax,nzmax,2)     temperature [K]                            *
  ! pv(0:nxmax,0:nymax,nzmax,2)     potential voriticity (pvu)                 *
  ! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll

  implicit none

  integer, save :: indice(3)=0,helpint=0,counter=0,izp,iz1,tempo(0:3)=0
  integer :: indj
  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
  integer :: rain_cloud_above,kz_inv
  real :: f_qvsat,pressure
  real :: rh,lsp,convp
  real :: uvzlev(nuvzmax),rhoh(nuvzmax),pinmconv(nzmax)
  real :: ew,pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: wzlev(nwzmax),uvwzlev(0:nxmax-1,0:nymax-1,nzmax)
  
  real,parameter :: const=r_air/ga

  logical :: init = .true.
 !below added by mc for test reason
  real latitudine, longitudine, quotap, quotaz
  character :: adate*8,atime*6
  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: skyp_check_w=1 !,timeintervalins=86400  !this flag test w
  character*7 :: stringwftime ! added for testing by mc
  
!C------------------ TESt FILE TO BE DELETED IN FINAL VERSION ------------

      if (skyp_check_w.eq.0) then
        write(stringwftime,'(I7.7)')wftime(indj)      
        open(79,file='..\windtrans\w_etadotdpdeta_'//stringwftime//'.dat') !        
      end if
      
!c-------------------------------------------------------------------------
! Determine current calendar date, needed for the file name
!**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss
  ! Open output file and write the output
  !**************************************
  !                                       
  !  open(137,file=path(2)(1:length(2))//'testveloc_'//adate// &  ! for testing u v interpolation against original data : by mc
  !       atime,form='formatted')
  !  
  !100 format(7e20.10)    
  !********* above added for test reSAON BY MC *****************
    
  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
  ! the vertical resolution in the z system is doubled. As reference point,*
  ! the lower left corner of the grid is used.                             *
  ! Unlike in the eta system, no difference between heights for u,v and    *
  ! heights for w exists.                                                  *
  !*************************************************************************


  !do kz=1,nuvz
  !  write (*,*) 'akz: ',akz(kz),'bkz',bkz(kz)
  !end do

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************

    do jy=0,nymin1
      do ix=0,nxmin1
        if (ps(ix,jy,1,n).gt.100000.) then
          ixm=ix
          jym=jy
          goto 3
        endif
      end do
    end do
3   continue


    tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n))/ &
         ps(ixm,jym,1,n))
    pold=ps(ixm,jym,1,n)
    height(1)=0.

    do kz=2,nuvz
      pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
      tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))


  ! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
  ! upon the transformation to z levels. In order to save computer memory, this is
  ! not done anymore in the standard version. However, this option can still be
  ! switched on by replacing the following lines with those below, that are
  ! currently commented out.
  ! Note that two more changes are necessary in this subroutine below.
  ! One change is also necessary in gridcheck.f, and another one in verttransform_nests.
  !*****************************************************************************

      if (abs(tv-tvold).gt.0.2) then
        height(kz)= &
             height(kz-1)+const*log(pold/pint)* &
             (tv-tvold)/log(tv/tvold)
      else
        height(kz)=height(kz-1)+ &
             const*log(pold/pint)*tv
      endif

  ! Switch on following lines to use doubled vertical resolution
  !*************************************************************
  !    if (abs(tv-tvold).gt.0.2) then
  !      height((kz-1)*2)=
  !    +      height(max((kz-2)*2,1))+const*log(pold/pint)*
  !    +      (tv-tvold)/log(tv/tvold)
  !    else
  !      height((kz-1)*2)=height(max((kz-2)*2,1))+
  !    +      const*log(pold/pint)*tv
  !    endif
  ! End doubled vertical resolution

      tvold=tv
      pold=pint
    end do


  ! Switch on following lines to use doubled vertical resolution
  !*************************************************************
  !  do 7 kz=3,nz-1,2
  !    height(kz)=0.5*(height(kz-1)+height(kz+1))
  !  height(nz)=height(nz-1)+height(nz-1)-height(nz-2)
  ! End doubled vertical resolution


  ! Determine highest levels that can be within PBL
  !************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        goto 9
      endif
    end do
9   continue

  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.

  endif !on init


  ! Loop over the whole grid
  !*************************

  do jy=0,nymin1
    do ix=0,nxmin1
      latitudine= xlon0+(ix)*dx  !added for test reason by mc
      longitudine=ylat0+(jy)*dy   !added for test reason by mc
        
        
      tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
      pold=ps(ix,jy,1,n)
      uvzlev(1)=0.
      wzlev(1)=0.
      rhoh(1)=pold/(r_air*tvold)
      !rhoperu(ix,jy,1)=uuh(ix,jy,1)*rhoh(1) !aaded for testing by mc 15-11-2013
      !rhoperv(ix,jy,1)=vvh(ix,jy,1)*rhoh(1) !aaded for testing by mc 15-11-2013
  ! Compute heights of eta levels
  !******************************

      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
        rhoh(kz)=pint/(r_air*tv)
        !rhoperu(ix,jy,kz)=uuh(ix,jy,kz)*rhoh(kz) !aaded for testing by mc 15-11-2013
        !rhoperv(ix,jy,kz)=vvh(ix,jy,kz)*rhoh(kz) !aaded for testing by mc 15-11-2013
        if (abs(tv-tvold).gt.0.2) then
          uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)* &
               (tv-tvold)/log(tv/tvold)
        else
          uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)*tv
        endif

        tvold=tv
        pold=pint
      end do


      do kz=2,netadot-1  !--NOTE: netadot is used in NORESM in ECMWF it was nwz. netadot is used since omega in Pa/s 
                         !--is co-located with u and etadot levels (akm, bkm) are not co-located and therefore used for computing pressur gradients
        wzlev(kz)=(uvzlev(kz+1)+uvzlev(kz))/2. !high of etadot levels
      end do
      wzlev(netadot)=wzlev(netadot-1)+ &
           uvzlev(nuvz)-uvzlev(nuvz-1)

      uvwzlev(ix,jy,1)=0.0
      do kz=2,nuvz
        uvwzlev(ix,jy,kz)=uvzlev(kz)
      end do
      
  ! doubling of vertical resolution is not implemeted yet in NORESM version --
  ! Switch on following lines to use doubled vertical resolution
  ! Switch off the three lines above.
  !*************************************************************
  !22          uvwzlev(ix,jy,(kz-1)*2)=uvzlev(kz)
  !     do 23 kz=2,nwz
  !23          uvwzlev(ix,jy,(kz-1)*2+1)=wzlev(kz)
  ! End doubled vertical resolution

  ! pinmconv=(h2-h1)/(p2-p1)

      pinmconv(1)=(uvzlev(2)-wzlev(1))/ &  !noRESM. note: akm are etadot levels staggered from u,v,omega levels see e.g. CAM3 user guide page 52
           ((akz(2)+bkz(2)*ps(ix,jy,1,n))- & !this is not actually used
           (akm(1)+bkm(1)*ps(ix,jy,1,n)))
      
      do kz=2,nz-1
        pinmconv(kz)=(wzlev(kz)-wzlev(kz-1))/ & !NORESM
             ((akm(kz)+bkm(kz)*ps(ix,jy,1,n))- &
             (akm(kz-1)+bkm(kz-1)*ps(ix,jy,1,n)))
      end do
      pinmconv(nz)=(wzlev(nz)-wzlev(nz-1))/ & !NORESM
           ((akm(nz)+bkm(nz)*ps(ix,jy,1,n))- &
           (akm(nz-1)+bkm(nz-1)*ps(ix,jy,1,n)))

  ! Levels, where u,v,t and q are given
  !************************************
      ww(ix,jy,1,n)=wwh(ix,jy,1) ! NORESM   note that here teh value of zero is assigned to vertical velocity at the ground   see below
      uu(ix,jy,1,n)=uuh(ix,jy,1)
      vv(ix,jy,1,n)=vvh(ix,jy,1)
      tt(ix,jy,1,n)=tth(ix,jy,1,n)
      qv(ix,jy,1,n)=qvh(ix,jy,1,n)
      pv(ix,jy,1,n)=pvh(ix,jy,1)
      rho(ix,jy,1,n)=rhoh(1)
      uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
      vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
      ww(ix,jy,nz,n)=wwh(ix,jy,nuvz)*pinmconv(nz) !NORESM noet that wwh(..,nz) has not been assigned zero here in readwind (zero assginment done in ECMWF)
      tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
      qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
      pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
      rho(ix,jy,nz,n)=rhoh(nuvz)
      kmin=2
      do iz=2,nz-1
        do kz=kmin,nuvz
          if(height(iz).gt.uvzlev(nuvz)) then
            uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
            vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
            tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
            qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
            pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
            rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
            goto 30
          endif
          if ((height(iz).gt.uvzlev(kz-1)).and. &
               (height(iz).le.uvzlev(kz))) then
            dz1=height(iz)-uvzlev(kz-1)
            dz2=uvzlev(kz)-height(iz)
            dz=dz1+dz2
            uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
            vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
            if  (height(iz).gt.uvzlev(2)) then                       !NORESM
              ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(kz-1)*dz2 &   !NORESM
                +wwh(ix,jy,kz)*pinmconv(kz)*dz1)/dz                  !NORESM
            else if (height(iz).le.uvzlev(2)) then                   !NORESM
              ww(ix,jy,iz,n)=wwh(ix,jy,2)*pinmconv(2)                !NORESM
              !else                                                  !NORESM
              !ww(ix,jy,iz,n)=0.                                     !NORESM ww will be assigned below
            end if
           
          ! write (137,100)1.*n, latitudine, longitudine, height(iz), uu(ix,jy,iz,n), vv(ix,jy,iz,n), ww(ix,jy,iz,n)  !ADDDEd FOR TEST REASON BY MC !TO BE REMOMEVED
           
            tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
                +tth(ix,jy,kz,n)*dz1)/dz
            qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
                +qvh(ix,jy,kz,n)*dz1)/dz
            pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
            rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
            kmin=kz
            goto 30
          endif
        end do
30      continue
      end do

      
    

  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/ &
           (height(2)-height(1))
      do kz=2,nz-1
        drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
             (height(kz+1)-height(kz-1))
      end do
      drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    end do
  end do


  !****************************************************************
  ! Compute slope of eta levels (in NORESM) (eta also in ECMWF) in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,ny-2
    do ix=1,nx-2

      kmin=2
      do iz=2,nz-1

        ui=uu(ix,jy,iz,n)*dxconst/cos((real(jy)*dy+ylat0)*pi180)
        vi=vv(ix,jy,iz,n)*dyconst
        

        do kz=kmin,nz
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
               (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            kmin=kz
            goto 47
          endif
        end do

47      ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1
        
        !drhoudx1=(rhoperu(ixp,jy,kl)-rhoperu(ix1,jy,kl))/2.*dxconst/cos((real(jy)*dy+ylat0)*pi180) !added for test by mc 15-11-2013
        !drhoudx2=(rhoperu(ixp,jy,klp)-rhoperu(ix1,jy,klp))/2.*dxconst/cos((real(jy)*dy+ylat0)*pi180) !added for test by mc 15-11-2013
        !drhoudx(ix,jy,iz)=(drhoudx1*dz2+drhoudx2*dz1)/dz !added for test by mc 15-11-2013
        
        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz
        
        !drhovdy1=(rhoperv(ix,jyp,kl)-rhoperv(ix,jy1,kl))/2.*dyconst !added for test by mc 15-11-2013
        !drhovdy2=(rhoperv(ix,jyp,klp)-rhoperv(ix,jy1,klp))/2.*dyconst !added for test by mc 15-11-2013
        !drhovdy(ix,jy,iz)=(drhovdy1*dz2+drhovdy2*dz1)/dz !added for test by mc 15-11-2013

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz
        
        if (height(iz).gt.uvwzlev(ix,jy,2)) then !  test 11-11-2013 for correction of w
          ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi) !valid also for NORESM
        else if (height(iz).le.uvwzlev(ix,jy,2)) then !!  test 11-11-2013 for correction of w 
          ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+dzdx2*uuh(ix,jy,2)*dxconst/cos((real(jy)*dy+ylat0)*pi180)+dzdy2*vvh(ix,jy,2)*dyconst !  test 11-11-2013 for correction of w
        end if !  test 11-11-2013 for correction of w
        !!else
        !!ww(ix,jy,iz,n)=0. ! to be assigned below for NORESM  
        !end if
        
      
        
      end do

    end do
  end do
    
    
   !-------additional part for NORESM use value of w in terrain following coordinate 
   !-------assign w=0 in terrain following coordinate at terrain level and, to fill the value not already set,
   !------- i.e value for height less then uvwzle(2), interpolate between zero and the first assigned value
   
    do jy=1,ny-2
      do ix=1,nx-2
        do iz=nz,2,-1
          if (height(iz).lt.uvwzlev(ix,jy,2)) then
            dz=uvwzlev(ix,jy,2)-height(1)
            dz1=height(iz)-height(1)
            dz2=dz-dz1
            !print *,ww(ix,jy,iz,n)
            ww(ix,jy,iz,n)=ww(ix,jy,1,n)*dz2/dz+ww(ix,jy,iz,n)*dz1/dz
            !print *,ww(ix,jy,iz,n)
          end if
        end do
      end do
    end do
   

  !C-------------- write some diagnostic by mc  --------------
   if (skyp_check_w.eq.0) then
     do jy=0,nymin1
       do ix=0, nxfield-1 
         do kz=1,26
           write (79,'(3f12.4,4E15.6)')xlon0+dx*ix,ylat0+dy*jy,height(kz),ww(ix,jy,kz,n)           
         end do
       end do
     end do
   end if !    
!*************************************************************************************

  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do jy=int(switchnorthg)-2,nymin1
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do


    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+ &
           vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else if (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

  do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
  end do

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do jy=0,int(switchsouthg)+3
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do

    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+ &
           vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else if (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do
  endif


  !write (*,*) 'initializing clouds, n:',n,nymin1,nxmin1,nz
  !   create a cloud and rainout/washout field, clouds occur where rh>80%
  !   total cloudheight is stored at level 0
  do jy=0,nymin1
    do ix=0,nxmin1
      rain_cloud_above=0
      lsp=lsprec(ix,jy,1,n)
      convp=convprec(ix,jy,1,n)
      cloudsh(ix,jy,n)=0
      do kz_inv=1,nz-1
         kz=nz-kz_inv+1
         pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
         rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
         clouds(ix,jy,kz,n)=0
         if (rh.gt.0.8) then ! in cloud
            if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
               rain_cloud_above=1
               cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+ &
                    height(kz)-height(kz-1)
               if (lsp.ge.convp) then
                  clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
               else
                  clouds(ix,jy,kz,n)=2 ! convp dominated rainout
               endif
            else ! no precipitation
                  clouds(ix,jy,kz,n)=1 ! cloud
            endif
         else ! no cloud
            if (rain_cloud_above.eq.1) then ! scavenging
               if (lsp.ge.convp) then
                  clouds(ix,jy,kz,n)=5 ! lsp dominated washout
               else
                  clouds(ix,jy,kz,n)=4 ! convp dominated washout
               endif
            endif
         endif
      end do
    end do
  end do

  !do 102 kz=1,nuvz
  !write(an,'(i02)') kz+10
  !write(*,*) nuvz,nymin1,nxmin1,'--',an,'--'
  !open(4,file='/nilu_wrk2/sec/cloudtest/cloud'//an,form='formatted')
  !do 101 jy=0,nymin1
  !    write(4,*) (clouds(ix,jy,kz,n),ix=1,nxmin1)
  !101   continue
  ! close(4)
  !102   continue

  ! open(4,file='/nilu_wrk2/sec/cloudtest/height',form='formatted')
  ! do 103 jy=0,nymin1
  !     write (4,*)
  !+       (height(kz),kz=1,nuvz)
  !103   continue
  ! close(4)

  !open(4,file='/nilu_wrk2/sec/cloudtest/p',form='formatted')
  ! do 104 jy=0,nymin1
  !     write (4,*)
  !+       (r_air*tt(ix,jy,1,n)*rho(ix,jy,1,n),ix=1,nxmin1)
  !104   continue
  ! close(4)
end subroutine verttransform
