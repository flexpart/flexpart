! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine richardson(psurf,ust,ttlev,qvlev,ulev,vlev,nuvz, &
       akz,bkz,hf,tt2,td2,h,wst,hmixplus,metdata_format)
  !                        i    i    i     i    i    i    i
  ! i   i  i   i   i  o  o     o
  !****************************************************************************
  !                                                                           *
  !     Calculation of mixing height based on the critical Richardson number. *
  !     Calculation of convective time scale.                                 *
  !     For unstable conditions, one iteration is performed. An excess        *
  !     temperature (dependent on hf and wst) is calculated, added to the     *
  !     temperature at the lowest model level. Then the procedure is repeated.*
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     22 August 1996                                                        *
  !                                                                           *
  !     Literature:                                                           *
  !     Vogelezang DHP and Holtslag AAM (1996): Evaluation and model impacts  *
  !     of alternative boundary-layer height formulations. Boundary-Layer     *
  !     Meteor. 81, 245-269.                                                  *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  !     Update: 1999-02-01 by G. Wotawa                                       *
  !                                                                           *
  !     Two meter level (temperature, humidity) is taken as reference level   *
  !     instead of first model level.                                         *
  !     New input variables tt2, td2 introduced.                              *
  !                                                                           *
  !     CHANGE: 17/11/2005 Caroline Forster NCEP GFS version                  * 
  !                                                                           *
  !     Unified ECMWF and GFS builds                                          *
  !     Marian Harustak, 12.5.2017                                            *
  !       - Merged richardson and richardson_gfs into one routine using       *
  !         if-then for meteo-type dependent code                             *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  ! h                          mixing height [m]                              *
  ! hf                         sensible heat flux                             *
  ! psurf                      surface pressure at point (xt,yt) [Pa]         *
  ! tv                         virtual temperature                            *
  ! wst                        convective velocity scale                      *
  ! metdata_format             format of metdata (ecmwf/gfs)                  *
  !                                                                           *
  ! Constants:                                                                *
  ! ric                        critical Richardson number                     *
  !                                                                           *
  !****************************************************************************

  use par_mod
  use class_gribfile

  implicit none

  integer :: metdata_format
  integer :: i,k,nuvz,iter,llev,loop_start
  real :: tv,tvold,zref,z,zold,pint,pold,theta,thetaref,ri
  real :: akz(nuvz),bkz(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2,ew
  real :: psurf,ust,ttlev(nuvz),qvlev(nuvz),h,excess
  real :: thetaold,zl,ul,vl,thetal,ril,hmixplus,wspeed,bvfsq,bvf
  real :: f_qvsat,rh,rhold,rhl,theta1,theta2,zl1,zl2,thetam
  real,parameter    :: const=r_air/ga, ric=0.25, b=100., bs=8.5
  integer,parameter :: itmax=3

  excess=0.0
  iter=0

  if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
    ! NCEP version: find first model level above ground
    !**************************************************

     llev = 0
     do i=1,nuvz
       if (psurf.lt.akz(i)) llev=i
     end do
     llev = llev+1
    ! sec llev should not be 1!
     if (llev.eq.1) llev = 2
     if (llev.gt.nuvz) llev = nuvz-1
    ! NCEP version
  end if


  ! Compute virtual temperature and virtual potential temperature at
  ! reference level (2 m)
  !*****************************************************************

30   iter=iter+1

  pold=psurf
  tvold=tt2*(1.+0.378*ew(td2)/psurf)
  zold=2.0
  zref=zold
  rhold=ew(td2)/ew(tt2)

  thetaref=tvold*(100000./pold)**(r_air/cpa)+excess
  thetaold=thetaref


  ! Integrate z up to one level above zt
  !*************************************
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    loop_start=2
  else
    loop_start=llev
  end if
  do k=loop_start,nuvz
    pint=akz(k)+bkz(k)*psurf  ! pressure on model layers
    tv=ttlev(k)*(1.+0.608*qvlev(k))

    if (abs(tv-tvold).gt.0.2) then
      z=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
    else
      z=zold+const*log(pold/pint)*tv
    endif

    theta=tv*(100000./pint)**(r_air/cpa)
  ! Petra
    rh = qvlev(k) / f_qvsat( pint, ttlev(k) )


  ! Calculate Richardson number at each level
  !****************************************

    ri=ga/thetaref*(theta-thetaref)*(z-zref)/ &
         max(((ulev(k)-ulev(2))**2+(vlev(k)-vlev(2))**2+b*ust**2),0.1)

  !  addition of second condition: MH should not be placed in an
  !  unstable layer (PS / Feb 2000)
    if (ri.gt.ric .and. thetaold.lt.theta) goto 20

    tvold=tv
    pold=pint
    rhold=rh
    thetaold=theta
    zold=z
  end do
  k=k-1 ! ESO: make sure k <= nuvz (ticket #139)
20   continue

  ! Determine Richardson number between the critical levels
  !********************************************************

  zl1=zold
  theta1=thetaold
  do i=1,20
    zl=zold+real(i)/20.*(z-zold)
    ul=ulev(k-1)+real(i)/20.*(ulev(k)-ulev(k-1))
    vl=vlev(k-1)+real(i)/20.*(vlev(k)-vlev(k-1))
    thetal=thetaold+real(i)/20.*(theta-thetaold)
    rhl=rhold+real(i)/20.*(rh-rhold)
    ril=ga/thetaref*(thetal-thetaref)*(zl-zref)/ &
         max(((ul-ulev(2))**2+(vl-vlev(2))**2+b*ust**2),0.1)
    zl2=zl
    theta2=thetal
    if (ril.gt.ric) goto 25
    zl1=zl
    theta1=thetal
  end do

25   continue
  h=zl
  thetam=0.5*(theta1+theta2)
  wspeed=sqrt(ul**2+vl**2)                    ! Wind speed at z=hmix
  bvfsq=(ga/thetam)*(theta2-theta1)/(zl2-zl1) ! Brunt-Vaisala frequency
                                              ! at z=hmix

  ! Under stable conditions, limit the maximum effect of the subgrid-scale topography
  ! by the maximum lifting possible from the available kinetic energy
  !*****************************************************************************

  if(bvfsq.le.0.) then
    hmixplus=9999.
  else
    bvf=sqrt(bvfsq)
    hmixplus=wspeed/bvf*convke
  endif


  ! Calculate convective velocity scale
  !************************************

  if (hf.lt.0.) then
    wst=(-h*ga/thetaref*hf/cpa)**0.333
    excess=-bs*hf/cpa/wst
    if (iter.lt.itmax) goto 30
  else
    wst=0.
  endif

end subroutine richardson
