subroutine centerofmass(xl,yl,n,xcenter,ycenter)
  !                        i  i  i    o       o
  !*****************************************************************************
  !                                                                            *
  !   This routine calculates the center of mass of n points on the Earth.     *
  !   Input are the longitudes (xl) and latitudes (yl) of the individual       *
  !   points, output is the longitude and latitude of the centre of mass.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 January 2002                                                        *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  integer :: n,l
  real :: xl(n),yl(n),xll,yll,xav,yav,zav,x,y,z,xcenter,ycenter


  xav=0.
  yav=0.
  zav=0.

  do l=1,n

  ! Convert longitude and latitude from degrees to radians
  !*******************************************************

    xll=xl(l)*pi180
    yll=yl(l)*pi180

  ! Calculate 3D coordinates from longitude and latitude
  !*****************************************************

    x = cos(yll)*sin(xll)
    y = -1.*cos(yll)*cos(xll)
    z = sin(yll)


  ! Find the mean location in Cartesian coordinates
  !************************************************

    xav=xav+x
    yav=yav+y
    zav=zav+z
  end do

  xav=xav/real(n)
  yav=yav/real(n)
  zav=zav/real(n)


  ! Project the point back onto Earth's surface
  !********************************************

  xcenter=atan2(xav,-1.*yav)
  ycenter=atan2(zav,sqrt(xav*xav+yav*yav))

  ! Convert back to degrees
  !************************

  xcenter=xcenter/pi180
  ycenter=ycenter/pi180

end subroutine centerofmass
