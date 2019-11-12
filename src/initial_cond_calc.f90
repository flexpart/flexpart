subroutine initial_cond_calc(itime,i)
  !                               i   i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the sensitivity to initial conditions for BW runs       *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     15 January 2010                                                        *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: itime,i,ix,jy,ixp,jyp,kz,ks
  integer :: il,ind,indz,indzp,nrelpointer
  real :: rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: ddx,ddy
  real :: rhoprof(2),rhoi,xl,yl,wx,wy,w
  integer :: mind2
  ! mind2        eso: pointer to 2nd windfield in memory


  ! For forward simulations, make a loop over the number of species;
  ! for backward simulations, make an additional loop over the release points
  !**************************************************************************


  if (itra1(i).ne.itime) return

  ! Depending on output option, calculate air density or set it to 1
  ! linit_cond: 1=mass unit, 2=mass mixing ratio unit
  !*****************************************************************


  if (linit_cond.eq.1) then     ! mass unit

    ix=int(xtra1(i))
    jy=int(ytra1(i))
    ixp=ix+1
    jyp=jy+1
    ddx=xtra1(i)-real(ix)
    ddy=ytra1(i)-real(jy)
    rddx=1.-ddx
    rddy=1.-ddy
    p1=rddx*rddy
    p2=ddx*rddy
    p3=rddx*ddy
    p4=ddx*ddy

    do il=2,nz
      if (height(il).gt.ztra1(i)) then
        indz=il-1
        indzp=il
        goto 6
      endif
    end do
6   continue

    dz1=ztra1(i)-height(indz)
    dz2=height(indzp)-ztra1(i)
    dz=1./(dz1+dz2)

  ! Take density from 2nd wind field in memory (accurate enough, no time interpolation needed)
  !*****************************************************************************
    mind2=memind(2)

    do ind=indz,indzp
      rhoprof(ind-indz+1)=p1*rho(ix ,jy ,ind,mind2) &
           +p2*rho(ixp,jy ,ind,mind2) &
           +p3*rho(ix ,jyp,ind,mind2) &
           +p4*rho(ixp,jyp,ind,mind2)
    end do
    rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz
  elseif (linit_cond.eq.2) then    ! mass mixing ratio unit
    rhoi=1.
  endif


  !****************************************************************************
  ! 1. Evaluate grid concentrations using a uniform kernel of bandwidths dx, dy
  !****************************************************************************


  ! For backward simulations, look from which release point the particle comes from
  ! For domain-filling trajectory option, npoint contains a consecutive particle
  ! number, not the release point information. Therefore, nrelpointer is set to 1
  ! for the domain-filling option.
  !*****************************************************************************

  if ((ioutputforeachrelease.eq.0).or.(mdomainfill.eq.1)) then
    nrelpointer=1
  else
    nrelpointer=npoint(i)
  endif

  do kz=1,numzgrid                ! determine height of cell
    if (outheight(kz).gt.ztra1(i)) goto 21
  end do
21   continue
  if (kz.le.numzgrid) then           ! inside output domain


    xl=(xtra1(i)*dx+xoutshift)/dxout
    yl=(ytra1(i)*dy+youtshift)/dyout
    ix=int(xl)
    if (xl.lt.0.) ix=ix-1
    jy=int(yl)
    if (yl.lt.0.) jy=jy-1


  ! If a particle is close to the domain boundary, do not use the kernel either
  !****************************************************************************

    if ((xl.lt.0.5).or.(yl.lt.0.5).or. &
         (xl.gt.real(numxgrid-1)-0.5).or. &
         (yl.gt.real(numygrid-1)-0.5)) then             ! no kernel, direct attribution to grid cell
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
           (jy.le.numygrid-1)) then
        do ks=1,nspec
          init_cond(ix,jy,kz,ks,nrelpointer)= &
               init_cond(ix,jy,kz,ks,nrelpointer)+ &
               xmass1(i,ks)/rhoi
        end do
      endif

    else                                 ! attribution via uniform kernel

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

      if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          w=wx*wy
          do ks=1,nspec
            init_cond(ix,jy,kz,ks,nrelpointer)= &
                 init_cond(ix,jy,kz,ks,nrelpointer)+xmass1(i,ks)/rhoi*w
          end do
        endif

        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=wx*(1.-wy)
          do ks=1,nspec
            init_cond(ix,jyp,kz,ks,nrelpointer)= &
                 init_cond(ix,jyp,kz,ks,nrelpointer)+xmass1(i,ks)/rhoi*w
          end do
        endif
      endif


      if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=(1.-wx)*(1.-wy)
          do ks=1,nspec
            init_cond(ixp,jyp,kz,ks,nrelpointer)= &
                 init_cond(ixp,jyp,kz,ks,nrelpointer)+xmass1(i,ks)/rhoi*w
          end do
        endif

        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          w=(1.-wx)*wy
          do ks=1,nspec
            init_cond(ixp,jy,kz,ks,nrelpointer)= &
                 init_cond(ixp,jy,kz,ks,nrelpointer)+xmass1(i,ks)/rhoi*w
          end do
        endif
      endif
    endif

  endif

end subroutine initial_cond_calc
