! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!  Taken from Press et al., Numerical Recipes

module random_mod
  
  implicit none

contains
  
  function ran1(idum)

    implicit none

    integer :: idum
    real    :: ran1
    integer,parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
    integer,parameter :: ntab=32, ndiv=1+(im-1)/ntab
    real,parameter    :: am=1./im, eps=1.2e-7, rnmx=1.-eps
    integer :: j, k
    integer :: iv(ntab) = (/ (0,j=1,ntab) /)
    integer :: iy=0

    if (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
      do j=ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        if (j.le.ntab) iv(j)=idum
      enddo
      iy=iv(1)
    endif
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    if (idum.lt.0) idum=idum+im
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=idum
    ran1=min(am*iy,rnmx)
  end function ran1


  function gasdev(idum)

    implicit none

    integer :: idum
    real    :: gasdev, fac, r, v1, v2
    integer :: iset = 0
    real    :: gset = 0.
    
    if (iset.eq.0) then
1     v1=2.*ran3(idum)-1.
      v2=2.*ran3(idum)-1.
      r=v1**2+v2**2
      if(r.ge.1.0 .or. r.eq.0.0) go to 1
      fac=sqrt(-2.*log(r)/r)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
    else
      gasdev=gset
      iset=0
    endif
  end function gasdev


  subroutine gasdev1(idum,random1,random2)

    implicit none

    integer :: idum
    real :: random1, random2, fac, v1, v2, r

1   v1=2.*ran3(idum)-1.
    v2=2.*ran3(idum)-1.
    r=v1**2+v2**2
    if(r.ge.1.0 .or. r.eq.0.0) go to 1
    fac=sqrt(-2.*log(r)/r)
    random1=v1*fac
    random2=v2*fac
! Limit the random numbers to lie within the interval -3 and +3
!**************************************************************
    if (random1.lt.-3.) random1=-3.
    if (random2.lt.-3.) random2=-3.
    if (random1.gt.3.) random1=3.
    if (random2.gt.3.) random2=3.
  end subroutine gasdev1


  function ran3(idum)

    implicit none

    integer :: idum
    real :: ran3

    integer,parameter :: mbig=1000000000, mseed=161803398, mz=0
    real,parameter    :: fac=1./mbig
    integer :: i,ii,inext,inextp,k
    integer :: mj,mk,ma(55)

    save inext,inextp,ma
    integer :: iff = 0

    if(idum.lt.0.or.iff.eq.0)then
      iff=1
      mj=mseed-iabs(idum)
      mj=mod(mj,mbig)
      ma(55)=mj
      mk=1
      do i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.mz)mk=mk+mbig
        mj=ma(ii)
      end do
      do k=1,4
        do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if(ma(i).lt.mz)ma(i)=ma(i)+mbig
        end do
      end do
      inext=0
      inextp=31
      idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz)mj=mj+mbig
    ma(inext)=mj
    ran3=mj*fac
  end function ran3
!  (C) Copr. 1986-92 Numerical Recipes Software US.

end module random_mod
