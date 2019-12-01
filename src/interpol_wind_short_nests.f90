subroutine interpol_wind_short_nests(itime,xt,yt,zt)
  !                                       i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod

  implicit none

  integer :: itime
  real :: xt,yt,zt

  ! Auxiliary variables needed for interpolation
  real :: dz1,dz2,dz
  real :: u1(2),v1(2),w1(2),uh(2),vh(2),wh(2)
  integer :: i,m,n,indexh,indzh


  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

  ! Calculate variables for time interpolation
  !*******************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  ! Determine the level below the current position for u,v
  !*******************************************************

  do i=2,nz
    if (height(i).gt.zt) then
      indz=i-1
      goto 6
    endif
  end do
6   continue


  ! Vertical distance to the level below and above current position
  !****************************************************************

  dz=1./(height(indz+1)-height(indz))
  dz1=(zt-height(indz))*dz
  dz2=(height(indz+1)-zt)*dz


  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Loop over 2 time steps and 2 levels
  !************************************

  do m=1,2
    indexh=memind(m)
    do n=1,2
      indzh=indz+n-1

      u1(n)=p1*uun(ix ,jy ,indzh,indexh,ngrid) &
           +p2*uun(ixp,jy ,indzh,indexh,ngrid) &
           +p3*uun(ix ,jyp,indzh,indexh,ngrid) &
           +p4*uun(ixp,jyp,indzh,indexh,ngrid)
      v1(n)=p1*vvn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*vvn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*vvn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*vvn(ixp,jyp,indzh,indexh,ngrid)
      w1(n)=p1*wwn(ix ,jy ,indzh,indexh,ngrid) &
           +p2*wwn(ixp,jy ,indzh,indexh,ngrid) &
           +p3*wwn(ix ,jyp,indzh,indexh,ngrid) &
           +p4*wwn(ixp,jyp,indzh,indexh,ngrid)

    end do


  !**********************************
  ! 2.) Linear vertical interpolation
  !**********************************

    uh(m)=dz2*u1(1)+dz1*u1(2)
    vh(m)=dz2*v1(1)+dz1*v1(2)
    wh(m)=dz2*w1(1)+dz1*w1(2)
  end do


  !************************************
  ! 3.) Temporal interpolation (linear)
  !************************************

  u=(uh(1)*dt2+uh(2)*dt1)*dtt
  v=(vh(1)*dt2+vh(2)*dt1)*dtt
  w=(wh(1)*dt2+wh(2)*dt1)*dtt

end subroutine interpol_wind_short_nests
