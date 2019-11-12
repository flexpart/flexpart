subroutine wetdepokernel(nunc,deposit,x,y,nage,kp)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     deposition fields using a uniform kernel with bandwidths dxout and dyout.*
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
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
  ! Changes:
  ! eso 10/2016: Added option to disregard kernel 
  ! 
  !*****************************************************************************

  use unc_mod
  use par_mod
  use com_mod

  implicit none

  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,nunc,nage,ks,kp

  xl=(x*dx+xoutshift)/dxout
  yl=(y*dy+youtshift)/dyout
  ix=int(xl)
  jy=int(yl)
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

  ! If no kernel is used, direct attribution to grid cell
  !******************************************************

  if (.not.lusekerneloutput) then
    do ks=1,nspec
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
           (jy.le.numygrid-1)) then
        wetgridunc(ix,jy,ks,kp,nunc,nage)= &
             wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
      end if
    end do
  else ! use kernel 
    
  ! Determine mass fractions for four grid points
  !**********************************************

  do ks=1,nspec

    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=wx*wy
      wetgridunc(ix,jy,ks,kp,nunc,nage)= &
           wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=(1.-wx)*(1.-wy)
      wetgridunc(ixp,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=(1.-wx)*wy
      wetgridunc(ixp,jy,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=wx*(1.-wy)
      wetgridunc(ix,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

  end do
  end if

end subroutine wetdepokernel
