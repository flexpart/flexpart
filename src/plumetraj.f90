subroutine plumetraj(itime)
  !                       i
  !*****************************************************************************
  !                                                                            *
  ! Determines a plume centroid trajectory for each release site, and manages  *
  ! clustering of particle locations. Certain parameters (average PV,          *
  ! tropopause height, etc., are provided along the plume trajectories.        *
  ! At the end, output is written to file 'trajectories.txt'.                  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 January 2002                                                        *
  !                                                                            *
  ! Variables:                                                                 *
  ! fclust          fraction of particles belonging to each cluster            *
  ! hmixcenter      mean mixing height for all particles                       *
  ! ncluster        number of clusters to be used                              *
  ! pvcenter        mean PV for all particles                                  *
  ! pvfract         fraction of particles with PV<2pvu                         *
  ! rms             total horizontal rms distance after clustering             *
  ! rmsdist         total horizontal rms distance before clustering            *
  ! rmsclust        horizontal rms distance for each individual cluster        *
  ! topocenter      mean topography underlying all particles                   *
  ! tropocenter     mean tropopause height at the positions of particles       *
  ! tropofract      fraction of particles within the troposphere               *
  ! zrms            total vertical rms distance after clustering               *
  ! zrmsdist        total vertical rms distance before clustering              *
  ! xclust,yclust,  Cluster centroid positions                                 *
  ! zclust                                                                     *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use mean_mod

  implicit none

  integer :: itime,ix,jy,ixp,jyp,indexh,i,j,k,m,n,il,ind,indz,indzp
  real :: xl(maxpart),yl(maxpart),zl(maxpart)
  real :: xcenter,ycenter,zcenter,dist,distance,rmsdist,zrmsdist

  real :: xclust(ncluster),yclust(ncluster),zclust(ncluster)
  real :: fclust(ncluster),rms,rmsclust(ncluster),zrms

  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,topocenter,hm(2),hmixi,hmixfract,hmixcenter
  real :: pv1(2),pvprof(2),pvi,pvcenter,pvfract,tr(2),tri,tropofract
  real :: tropocenter


  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  ! Loop about all release points
  !******************************

  do j=1,numpoint
    if (abs(ireleasestart(j)-itime).gt.lage(nageclass)) goto 10
    topocenter=0.
    hmixcenter=0.
    hmixfract=0.
    tropocenter=0.
    tropofract=0.
    pvfract=0.
    pvcenter=0.
    rmsdist=0.
    zrmsdist=0.

    n=0
    do i=1,numpart
      if (itra1(i).ne.itime) goto 20
      if (npoint(i).ne.j) goto 20
      n=n+1
      xl(n)=xlon0+xtra1(i)*dx
      yl(n)=ylat0+ytra1(i)*dy
      zl(n)=ztra1(i)


  ! Interpolate PBL height, PV, and tropopause height to each
  ! particle position in order to determine fraction of particles
  ! within the PBL, above tropopause height, and average PV.
  ! Interpolate topography, too, and convert to altitude asl
  !**************************************************************

      ix=int(xtra1(i))
      jy=int(ytra1(i))
      ixp=ix+1
      jyp=jy+1
      ddx=xtra1(i)-real(ix)
      ddy=ytra1(i)-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

  ! Topography
  !***********

      topo=p1*oro(ix ,jy) &
           + p2*oro(ixp,jy) &
           + p3*oro(ix ,jyp) &
           + p4*oro(ixp,jyp)
      topocenter=topocenter+topo

  ! Potential vorticity
  !********************

      do il=2,nz
        if (height(il).gt.zl(n)) then
          indz=il-1
          indzp=il
          goto 6
        endif
      end do
6     continue

      dz1=zl(n)-height(indz)
      dz2=height(indzp)-zl(n)
      dz=1./(dz1+dz2)


      do ind=indz,indzp
        do m=1,2
          indexh=memind(m)
          pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
               +p2*pv(ixp,jy ,ind,indexh) &
               +p3*pv(ix ,jyp,ind,indexh) &
               +p4*pv(ixp,jyp,ind,indexh)
        end do
        pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
      end do
      pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
      pvcenter=pvcenter+pvi
      if (yl(n).gt.0.) then
        if (pvi.lt.2.) pvfract=pvfract+1.
      else
        if (pvi.gt.-2.) pvfract=pvfract+1.
      endif


  ! Tropopause and PBL height
  !**************************

      do m=1,2
        indexh=memind(m)

        tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
             + p2*tropopause(ixp,jy ,1,indexh) &
             + p3*tropopause(ix ,jyp,1,indexh) &
             + p4*tropopause(ixp,jyp,1,indexh)

        hm(m)=p1*hmix(ix ,jy ,1,indexh) &
             + p2*hmix(ixp,jy ,1,indexh) &
             + p3*hmix(ix ,jyp,1,indexh) &
             + p4*hmix(ixp,jyp,1,indexh)
      end do

      hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
      tri=(tr(1)*dt2+tr(2)*dt1)*dtt
      if (zl(n).lt.tri) tropofract=tropofract+1.
      tropocenter=tropocenter+tri+topo
      if (zl(n).lt.hmixi) hmixfract=hmixfract+1.
      zl(n)=zl(n)+topo        ! convert to height asl
      hmixcenter=hmixcenter+hmixi


20    continue
    end do


  ! Make statistics for all plumes with n>0 particles
  !**************************************************

    if (n.gt.0) then
      topocenter=topocenter/real(n)
      hmixcenter=hmixcenter/real(n)
      pvcenter=pvcenter/real(n)
      tropocenter=tropocenter/real(n)
      hmixfract=100.*hmixfract/real(n)
      pvfract=100.*pvfract/real(n)
      tropofract=100.*tropofract/real(n)

  ! Cluster the particle positions
  !*******************************

      call clustering(xl,yl,zl,n,xclust,yclust,zclust,fclust,rms, &
           rmsclust,zrms)


  ! Determine center of mass position on earth and average height
  !**************************************************************

      call centerofmass(xl,yl,n,xcenter,ycenter)
      call mean(zl,zcenter,zrmsdist,n)

  ! Root mean square distance from center of mass
  !**********************************************

      do k=1,n
        dist=distance(yl(k),xl(k),ycenter,xcenter)
        rmsdist=rmsdist+dist*dist
      end do
      if (rmsdist.gt.0.) rmsdist=sqrt(rmsdist/real(n))
      rmsdist=max(rmsdist,0.)

  ! Write out results in trajectory data file
  !******************************************

      write(unitouttraj,'(i5,i8,2f9.4,4f8.1,f8.2,4f8.1,3f6.1,&
           &5(2f8.3,f7.0,f6.1,f8.1))')&
           &j,itime-(ireleasestart(j)+ireleaseend(j))/2, &
           xcenter,ycenter,zcenter,topocenter,hmixcenter,tropocenter, &
           pvcenter,rmsdist,rms,zrmsdist,zrms,hmixfract,pvfract, &
           tropofract, &
           (xclust(k),yclust(k),zclust(k),fclust(k),rmsclust(k), &
           k=1,ncluster)
    endif


10  continue
  end do


end subroutine plumetraj
