! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine drydepokernel(nunc,deposit,x,y,nage,kp)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition to the grid using a uniform kernel with  *
  !     bandwidths dx and dy.                                                  *
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

  real(dep_prec), dimension(maxspec) :: deposit
  real :: x,y,ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,nunc,nage,kp


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
      if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then
        if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
             (jy.le.numygrid-1)) then
          drygridunc(ix,jy,ks,kp,nunc,nage)= &
               drygridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
        end if
      end if
    end do
  else ! use kernel 


  ! Determine mass fractions for four grid points
  !**********************************************
  do ks=1,nspec

   if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then

      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
        (jy.le.numygrid-1)) then
        w=wx*wy
        drygridunc(ix,jy,ks,kp,nunc,nage)= &
           drygridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
     endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
    w=(1.-wx)*(1.-wy)
      drygridunc(ixp,jyp,ks,kp,nunc,nage)= &
           drygridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=(1.-wx)*wy
      drygridunc(ixp,jy,ks,kp,nunc,nage)= &
           drygridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=wx*(1.-wy)
      drygridunc(ix,jyp,ks,kp,nunc,nage)= &
           drygridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
    endif

    endif ! deposit>0
  end do
end if

end subroutine drydepokernel
