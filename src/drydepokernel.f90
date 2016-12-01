!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

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

  use unc_mod
  use par_mod
  use com_mod

  implicit none

  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
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


  ! Determine mass fractions for four grid points
  !**********************************************
    do ks=1,nspec

    if ((abs(deposit(ks)).gt.0).and.DRYDEPSPEC(ks)) then

    if (.not.usekernel) then
       drygridunc(ix,jy,ks,kp,nunc,nage)= &
           drygridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
    else
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

    endif ! kernel
    endif ! deposit>0

    end do

end subroutine drydepokernel
