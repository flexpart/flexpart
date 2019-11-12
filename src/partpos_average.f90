
subroutine partpos_average(itime,j)


!**********************************************************************
! This subroutine averages particle quantities, to be used for particle
! dump (in partoutput.f90). Averaging is done over output interval.
!**********************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: itime,j,ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: xlon,ylat,x,y,z
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,hm(2),hmixi,pv1(2),pvprof(2),pvi,qv1(2),qvprof(2),qvi
  real :: tt1(2),ttprof(2),tti,rho1(2),rhoprof(2),rhoi
  real :: uu1(2),uuprof(2),uui,vv1(2),vvprof(2),vvi
  real :: tr(2),tri,energy



 ! Some variables needed for temporal interpolation
  !*************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  xlon=xlon0+xtra1(j)*dx
  ylat=ylat0+ytra1(j)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

  ix=xtra1(j)
  jy=ytra1(j)
  ixp=ix+1
  jyp=jy+1
  ddx=xtra1(j)-real(ix)
  ddy=ytra1(j)-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

! eso: Temporary fix for particle exactly at north pole
  if (jyp >= nymax) then
  !  write(*,*) 'WARNING: conccalc.f90 jyp >= nymax'
    jyp=jyp-1
  end if

  ! Topography
  !***********

  topo=p1*oro(ix,jy)+p2*oro(ixp,jy)+p3*oro(ix,jyp)+p4*oro(ixp,jyp)

 ! Potential vorticity, specific humidity, temperature, and density
  !*****************************************************************

  do il=2,nz
    if (height(il).gt.ztra1(j)) then
      indz=il-1
      indzp=il
      goto 6
    endif
  end do
6 continue

  dz1=ztra1(j)-height(indz)
  dz2=height(indzp)-ztra1(j)
  dz=1./(dz1+dz2)


  do ind=indz,indzp
    do m=1,2
      indexh=memind(m)

  ! Potential vorticity
      pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
            +p2*pv(ixp,jy ,ind,indexh) &
            +p3*pv(ix ,jyp,ind,indexh) &
            +p4*pv(ixp,jyp,ind,indexh)
  ! Specific humidity
      qv1(m)=p1*qv(ix ,jy ,ind,indexh) &
            +p2*qv(ixp,jy ,ind,indexh) &
            +p3*qv(ix ,jyp,ind,indexh) &
            +p4*qv(ixp,jyp,ind,indexh)
  ! Temperature
      tt1(m)=p1*tt(ix ,jy ,ind,indexh) &
            +p2*tt(ixp,jy ,ind,indexh) &
            +p3*tt(ix ,jyp,ind,indexh) &
            +p4*tt(ixp,jyp,ind,indexh)
  ! U wind
      uu1(m)=p1*uu(ix ,jy ,ind,indexh) &
            +p2*uu(ixp,jy ,ind,indexh) &
            +p3*uu(ix ,jyp,ind,indexh) &
            +p4*uu(ixp,jyp,ind,indexh)
  ! V wind
      vv1(m)=p1*vv(ix ,jy ,ind,indexh) &
            +p2*vv(ixp,jy ,ind,indexh) &
            +p3*vv(ix ,jyp,ind,indexh) &
            +p4*vv(ixp,jyp,ind,indexh)
  ! Density
      rho1(m)=p1*rho(ix ,jy ,ind,indexh) &
             +p2*rho(ixp,jy ,ind,indexh) &
             +p3*rho(ix ,jyp,ind,indexh) &
             +p4*rho(ixp,jyp,ind,indexh)
    end do
    pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
    qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
    ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
    uuprof(ind-indz+1)=(uu1(1)*dt2+uu1(2)*dt1)*dtt
    vvprof(ind-indz+1)=(vv1(1)*dt2+vv1(2)*dt1)*dtt
    rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
  end do
  pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
  qvi=(dz1*qvprof(2)+dz2*qvprof(1))*dz
  tti=(dz1*ttprof(2)+dz2*ttprof(1))*dz
  uui=(dz1*uuprof(2)+dz2*uuprof(1))*dz
  vvi=(dz1*vvprof(2)+dz2*vvprof(1))*dz
  rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz

  ! Tropopause and PBL height
  !**************************

  do m=1,2
    indexh=memind(m)

  ! Tropopause
    tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
        + p2*tropopause(ixp,jy ,1,indexh) &
        + p3*tropopause(ix ,jyp,1,indexh) &
        + p4*tropopause(ixp,jyp,1,indexh)

  ! PBL height
    hm(m)=p1*hmix(ix ,jy ,1,indexh) &
        + p2*hmix(ixp,jy ,1,indexh) &
        + p3*hmix(ix ,jyp,1,indexh) &
        + p4*hmix(ixp,jyp,1,indexh)
  end do

  hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
  tri=(tr(1)*dt2+tr(2)*dt1)*dtt


  energy=tti*cpa+(ztra1(j)+topo)*9.81+qvi*2501000.+(uui**2+vvi**2)/2.

  ! Add new values to sum and increase counter by one
  !**************************************************

  npart_av(j)=npart_av(j)+1

  ! Calculate Cartesian 3D coordinates suitable for averaging
  !**********************************************************

  xlon=xlon*pi180
  ylat=ylat*pi180
  x = cos(ylat)*sin(xlon)
  y = -1.*cos(ylat)*cos(xlon)
  z = sin(ylat)

  part_av_cartx(j)=part_av_cartx(j)+x
  part_av_carty(j)=part_av_carty(j)+y
  part_av_cartz(j)=part_av_cartz(j)+z
  part_av_z(j)=part_av_z(j)+ztra1(j)
  part_av_topo(j)=part_av_topo(j)+topo
  part_av_pv(j)=part_av_pv(j)+pvi
  part_av_qv(j)=part_av_qv(j)+qvi
  part_av_tt(j)=part_av_tt(j)+tti
  part_av_uu(j)=part_av_uu(j)+uui
  part_av_vv(j)=part_av_vv(j)+vvi
  part_av_rho(j)=part_av_rho(j)+rhoi
  part_av_tro(j)=part_av_tro(j)+tri
  part_av_hmix(j)=part_av_hmix(j)+hmixi
  part_av_energy(j)=part_av_energy(j)+energy


return
end subroutine partpos_average 
