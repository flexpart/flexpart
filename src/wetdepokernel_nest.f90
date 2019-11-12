subroutine wetdepokernel_nest(nunc,deposit,x,y,nage,kp)
  !                           i    i       i i i    i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     nested deposition fields using a uniform kernel with bandwidths        *
  !     dxoutn and dyoutn.                                                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !      2 September 2004: Adaptation from wetdepokernel.                      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use par_mod
  use com_mod

  implicit none

  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,kp,nunc,nage

  xl=(x*dx+xoutshiftn)/dxoutn
  yl=(y*dy+youtshiftn)/dyoutn

  ! old:
  ! ix=int(xl) 
  ! jy=int(yl)

  ! ESO: for xl,yl in range <-.5,-1> we get ix,jy=0 and thus negative
  ! wx,wy as the int function rounds upwards for negative numbers.
  ! Either use the floor function, or (perhaps more correctly?) use "(xl.gt.-0.5)" 
  ! in place of "(ix.ge.0)" and similar for the upper boundary.

  ! new:
  ix=floor(xl)
  jy=floor(yl)

  ddx=xl-real(ix)                   ! distance to left cell border
  ddy=yl-real(jy)                   ! distance to lower cell border


  if (ddx.gt.0.5) then
    ixp=ix+1
    wx=1.5-ddx
  else
    ixp=ix-1
    wx=0.5+ddx
  endif

  if (ddy.gt.0.5) then
    jyp=jy+1
    wy=1.5-ddy
  else
    jyp=jy-1
    wy=0.5+ddy
  endif

! Determine mass fractions for four grid points
!**********************************************

  do ks=1,nspec
    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=wx*wy
      wetgriduncn(ix,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=(1.-wx)*(1.-wy)
      wetgriduncn(ixp,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=(1.-wx)*wy
      wetgriduncn(ixp,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=wx*(1.-wy)
      wetgriduncn(ix,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

  end do
end subroutine wetdepokernel_nest
