! From numerical recipes
! Change by A. Stohl: Use of integer instead of real values

subroutine sort2(n,arr,brr)

  implicit none

  integer :: n
  integer :: arr(n),brr(n)
  integer,parameter :: m=7,nstack=50
  integer :: i,ir,j,jstack,k,l,istack(nstack)
  integer :: a,b,temp
  jstack=0
  l=1
  ir=n
1   if(ir-l.lt.m)then
    do j=l+1,ir
      a=arr(j)
      b=brr(j)
      do i=j-1,1,-1
        if(arr(i).le.a)goto 2
        arr(i+1)=arr(i)
        brr(i+1)=brr(i)
      end do
      i=0
2     arr(i+1)=a
      brr(i+1)=b
    end do
    if(jstack.eq.0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+ir)/2
    temp=arr(k)
    arr(k)=arr(l+1)
    arr(l+1)=temp
    temp=brr(k)
    brr(k)=brr(l+1)
    brr(l+1)=temp
    if(arr(l+1).gt.arr(ir))then
      temp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=temp
      temp=brr(l+1)
      brr(l+1)=brr(ir)
      brr(ir)=temp
    endif
    if(arr(l).gt.arr(ir))then
      temp=arr(l)
      arr(l)=arr(ir)
      arr(ir)=temp
      temp=brr(l)
      brr(l)=brr(ir)
      brr(ir)=temp
    endif
    if(arr(l+1).gt.arr(l))then
      temp=arr(l+1)
      arr(l+1)=arr(l)
      arr(l)=temp
      temp=brr(l+1)
      brr(l+1)=brr(l)
      brr(l)=temp
    endif
    i=l+1
    j=ir
    a=arr(l)
    b=brr(l)
3   continue
      i=i+1
    if(arr(i).lt.a)goto 3
4   continue
      j=j-1
    if(arr(j).gt.a)goto 4
    if(j.lt.i)goto 5
    temp=arr(i)
    arr(i)=arr(j)
    arr(j)=temp
    temp=brr(i)
    brr(i)=brr(j)
    brr(j)=temp
    goto 3
5   arr(l)=arr(j)
    arr(j)=a
    brr(l)=brr(j)
    brr(j)=b
    jstack=jstack+2
    if(jstack.gt.nstack) then
       print*, 'nstack too small in sort2'
       stop
    end if
    if(ir-i+1.ge.j-l)then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    endif
  endif
  goto 1
end subroutine sort2
!  (C) Copr. 1986-92 Numerical Recipes Software us.
