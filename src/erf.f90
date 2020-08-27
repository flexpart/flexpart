! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

! To be used, if the non-standard Fortran function erf does not exist on
! your machine
!
!aus:  Numerical Recipes (FORTRAN) / Chapter 6.
!
!6.1  FUNCTION GAMMLN
!6.2  FUNCTION GAMMP   <6.2:GSER/6.2:GCF/6.1:GAMMLN>
!6.2  FUNCTION GAMMQ   <6.2:GSER/6.2:GCF/6.1:GAMMLN>
!6.2  SUBROUTINE GSER    <6.1:GAMMLN>
!6.2  SUBROUTINE GCF     <6.1:GAMMLN>
!6.2  FUNCTION ERF     <6.2:GAMMP/6.2:GSER/6.2:GCF/6.1:GAMMLN>
!6.2  FUNCTION ERFC    <6.2.:GAMMP/6.2:GAMMQ/6.2:GSER/
!                       6.2:GCF/6.1:GAMMLN>
!6.2  FUNCTION ERFCC

function gammln(xx)

  use par_mod, only: dp

  implicit none

  integer :: j
  real :: x,tmp,ser,xx,gammln
  real :: cof(6) = (/ &
       76.18009173_dp, -86.50532033_dp, 24.01409822_dp, &
       -1.231739516_dp, .120858003e-2_dp, -.536382e-5_dp    /)
  real :: stp = 2.50662827465_dp
  real :: half = 0.5_dp, one = 1.0_dp, fpf = 5.5_dp

  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  end do
  gammln=tmp+log(stp*ser)
end function gammln

function gammp(a,x)

  implicit none

  real :: a, x, gln, gamser, gammp, gammcf

  if(x .lt. 0. .or. a .le. 0.) then
     print*, 'gammp'
     stop
  end if
  if(x.lt.a+1.)then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.-gammcf
  endif
end function gammp

function gammq(a,x)

  implicit none

  real :: a, x, gln, gamser, gammq, gammcf

  if(x.lt.0..or.a.le.0.) then
     print*, 'gammq'
     stop
  end if
  if(x.lt.a+1.)then
    call gser(gamser,a,x,gln)
    gammq=1.-gamser
  else
    call gcf(gammcf,a,x,gln)
    gammq=gammcf
  endif
end function gammq

subroutine gser(gamser,a,x,gln)

  implicit none

  integer :: n
  real :: gamser, a, x, gln, ap, summ, del
  real, external :: gammln

  integer,parameter :: itmax=100
  real,parameter    :: eps=3.e-7

  gln=gammln(a)
  if(x.le.0.)then
    if(x.lt.0.) then
       print*, 'gser'
       stop
    end if
    gamser=0.
    return
  endif
  ap=a
  summ=1./a
  del=summ
  do n=1,itmax
    ap=ap+1.
    del=del*x/ap
    summ=summ+del
    if(abs(del).lt.abs(summ)*eps)go to 1
  end do
  print*, 'gser: a too large, itmax too small'
  stop
1   gamser=summ*exp(-x+a*log(x)-gln)
end subroutine gser

subroutine gcf(gammcf,a,x,gln)

  implicit none

  integer :: n
  real :: gammcf, a, x, gln, gold, a0, a1, b0, b1, fac, an, anf, ana, g
  real, external :: gammln

  integer,parameter :: itmax=100
  real,parameter    :: eps=3.e-7

  gln=gammln(a)
  gold=0.
  a0=1.
  a1=x
  b0=0.
  b1=1.
  fac=1.
  do n=1,itmax
    an=real(n)
    ana=an-a
    a0=(a1+a0*ana)*fac
    b0=(b1+b0*ana)*fac
    anf=an*fac
    a1=x*a0+anf*a1
    b1=x*b0+anf*b1
    if(a1.ne.0.)then
      fac=1./a1
      g=b1*fac
      if(abs((g-gold)/g).lt.eps)go to 1
      gold=g
    endif
  end do
  print*, 'gcf: a too large, itmax too small'
  stop
1   gammcf=exp(-x+a*alog(x)-gln)*g
end subroutine gcf

function erf(x)

  implicit none

  real :: x, erf
  real, external :: gammp

  if(x.lt.0.)then
    erf=-gammp(.5,x**2)
  else
    erf=gammp(.5,x**2)
  endif
end function erf

function erfc(x)

  implicit none

  real :: x, erfc
  real, external :: gammp, gammq

  if(x.lt.0.)then
    erfc=1.+gammp(.5,x**2)
  else
    erfc=gammq(.5,x**2)
  endif
end function erfc

function erfcc(x)

  implicit none

  real :: x, z, t, erfcc

  z=abs(x)
  t=1./(1.+0.5*z)
  erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
       t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
       t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  if (x.lt.0.) erfcc=2.-erfcc
end function erfcc
