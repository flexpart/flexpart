! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

real function ew(x)

  !****************************************************************
  !SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
  !NACH DER GOFF-GRATCH-FORMEL.
  !****************************************************************

  implicit none

  real :: x, y, a, c, d

  ew=0.
  if(x.le.0.) stop 'sorry: t not in [k]'
  y=373.16/x
  a=-7.90298*(y-1.)
  a=a+(5.02808*0.43429*alog(y))
  c=(1.-(1./y))*11.344
  c=-1.+(10.**c)
  c=-1.3816*c/(10.**7)
  d=(1.-y)*3.49149
  d=-1.+(10.**d)
  d=8.1328*d/(10.**3)
  y=a+c+d
  ew=101324.6*(10.**y)       ! Saettigungsdampfdruck in Pa

end function ew
