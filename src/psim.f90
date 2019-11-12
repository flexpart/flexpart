real function psim(z,al)

  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION: CALCULATION OF THE STABILITY CORRECTION FUNCTION FOR   *
  !              MOMENTUM AS FUNCTION OF HEIGHT Z AND OBUKHOV SCALE     *
  !              HEIGHT L                                               *
  !                                                                     *
  !**********************************************************************

  use par_mod

  implicit none

  real :: z,al,zeta,x,a1,a2

  zeta=z/al
  if(zeta.le.0.) then
  ! UNSTABLE CASE
    x=(1.-15.*zeta)**0.25
    a1=((1.+x)/2.)**2
    a2=(1.+x**2)/2.
    psim=log(a1*a2)-2.*atan(x)+pi/2.
  else
  ! STABLE CASE
    psim=-4.7*zeta
  endif

end function psim
