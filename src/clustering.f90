subroutine clustering(xl,yl,zl,n,xclust,yclust,zclust,fclust,rms, &
       rmsclust,zrms)
  !                      i  i  i  i   o      o      o      o     o
  !   o      o
  !*****************************************************************************
  !                                                                            *
  !   This routine clusters the particle position into ncluster custers.       *
  !   Input are the longitudes (xl) and latitudes (yl) of the individual       *
  !   points, output are the cluster mean positions (xclust,yclust).           *
  !   Vertical positions are not directly used for the clustering.             *
  !                                                                            *
  !   For clustering, the procedure described in Dorling et al. (1992) is used.*
  !                                                                            *
  !   Dorling, S.R., Davies, T.D. and Pierce, C.E. (1992):                     *
  !   Cluster analysis: a technique for estimating the synoptic meteorological *
  !   controls on air and precipitation chemistry - method and applications.   *
  !   Atmospheric Environment 26A, 2575-2581.                                  *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 February 2002                                                        *
  !                                                                            *
  ! Variables:                                                                 *
  ! fclust          fraction of particles belonging to each cluster            *
  ! ncluster        number of clusters to be used                              *
  ! rms             total horizontal rms distance after clustering             *
  ! rmsclust        horizontal rms distance for each individual cluster        *
  ! zrms            total vertical rms distance after clustering               *
  ! xclust,yclust,  Cluster centroid positions                                 *
  ! zclust                                                                     *
  ! xl,yl,zl        particle positions                                         *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  integer :: n,i,j,l,nclust(maxpart),numb(ncluster),ncl
  real :: xl(n),yl(n),zl(n),xclust(ncluster),yclust(ncluster),x,y,z
  real :: zclust(ncluster),distance2,distances,distancemin,rms,rmsold
  real :: xav(ncluster),yav(ncluster),zav(ncluster),fclust(ncluster)
  real :: rmsclust(ncluster)
  real :: zdist,zrms



  if (n.lt.ncluster) return
  rmsold=-5.

  ! Convert longitude and latitude from degrees to radians
  !*******************************************************

  do i=1,n
    nclust(i)=i
    xl(i)=xl(i)*pi180
    yl(i)=yl(i)*pi180
  end do


  ! Generate a seed for each cluster
  !*********************************

  do j=1,ncluster
    zclust(j)=0.
    xclust(j)=xl(j*n/ncluster)
    yclust(j)=yl(j*n/ncluster)
  end do


  ! Iterative loop to compute the cluster means
  !********************************************

  do l=1,100

  ! Assign each particle to a cluster: criterion minimum distance to the
  ! cluster mean position
  !*********************************************************************


    do i=1,n
      distancemin=10.**10.
      do j=1,ncluster
        distances=distance2(yl(i),xl(i),yclust(j),xclust(j))
        if (distances.lt.distancemin) then
          distancemin=distances
          ncl=j
        endif
      end do
      nclust(i)=ncl
    end do


  ! Recalculate the cluster centroid position: convert to 3D Cartesian coordinates,
  ! calculate mean position, and re-project this point onto the Earth's surface
  !*****************************************************************************

    do j=1,ncluster
      xav(j)=0.
      yav(j)=0.
      zav(j)=0.
      rmsclust(j)=0.
      numb(j)=0
    end do
    rms=0.

    do i=1,n
      numb(nclust(i))=numb(nclust(i))+1
      distances=distance2(yl(i),xl(i), &
           yclust(nclust(i)),xclust(nclust(i)))

  ! rms is the total rms of all particles
  ! rmsclust is the rms for a particular cluster
  !*********************************************

      rms=rms+distances*distances
      rmsclust(nclust(i))=rmsclust(nclust(i))+distances*distances

  ! Calculate Cartesian 3D coordinates from longitude and latitude
  !***************************************************************

      x = cos(yl(i))*sin(xl(i))
      y = -1.*cos(yl(i))*cos(xl(i))
      z = sin(yl(i))
      xav(nclust(i))=xav(nclust(i))+x
      yav(nclust(i))=yav(nclust(i))+y
      zav(nclust(i))=zav(nclust(i))+z
    end do

    rms=sqrt(rms/real(n))


  ! Find the mean location in Cartesian coordinates
  !************************************************

    do j=1,ncluster
      if (numb(j).gt.0) then
        rmsclust(j)=sqrt(rmsclust(j)/real(numb(j)))
        xav(j)=xav(j)/real(numb(j))
        yav(j)=yav(j)/real(numb(j))
        zav(j)=zav(j)/real(numb(j))

  ! Project the point back onto Earth's surface
  !********************************************

        xclust(j)=atan2(xav(j),-1.*yav(j))
        yclust(j)=atan2(zav(j),sqrt(xav(j)*xav(j)+yav(j)*yav(j)))
      endif
    end do


  ! Leave the loop if the RMS distance decreases only slightly between 2 iterations
  !*****************************************************************************

    if ((l.gt.1).and.(abs(rms-rmsold)/rmsold.lt.0.005)) goto 99
    rmsold=rms

  end do

99   continue

  ! Convert longitude and latitude from radians to degrees
  !*******************************************************

  do i=1,n
    xl(i)=xl(i)/pi180
    yl(i)=yl(i)/pi180
    zclust(nclust(i))=zclust(nclust(i))+zl(i)
  end do

  do j=1,ncluster
    xclust(j)=xclust(j)/pi180
    yclust(j)=yclust(j)/pi180
    if (numb(j).gt.0) zclust(j)=zclust(j)/real(numb(j))
    fclust(j)=100.*real(numb(j))/real(n)
  end do

  ! Determine total vertical RMS deviation
  !***************************************

  zrms=0.
  do i=1,n
    zdist=zl(i)-zclust(nclust(i))
    zrms=zrms+zdist*zdist
  end do
  if (zrms.gt.0.) zrms=sqrt(zrms/real(n))

end subroutine clustering
