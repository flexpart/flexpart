! Function calculates dynamic viscosity of air (kg/m/s) as function of
! temperature (K) using Sutherland's formula

real function viscosity(t)

  implicit none

  real :: t
  real,parameter :: c=120.,t_0=291.15,eta_0=1.827e-5

  viscosity=eta_0*(t_0+c)/(t+c)*(t/t_0)**1.5

  return

end function viscosity
