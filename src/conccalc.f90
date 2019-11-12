subroutine conccalc(itime,weight)
  !                      i     i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the concentrations on a regular grid using volume       *
  !     sampling                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1996                                                            *
  !                                                                            *
  !     April 2000: Update to calculate age spectra                            *
  !                 Bug fix to avoid negative conc. at the domain boundaries,  *
  !                 as suggested by Petra Seibert                              *
  !                                                                            *
  !     2 July 2002: re-order if-statements in order to optimize CPU time      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nspeciesdim     = nspec for forward runs, 1 for backward runs              *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: itime,itage,i,ix,jy,ixp,jyp,kz,ks,n,nage
  integer :: il,ind,indz,indzp,nrelpointer
  real :: rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: weight,hx,hy,hz,h,xd,yd,zd,xkern,r2,c(maxspec),ddx,ddy
  real :: rhoprof(2),rhoi
  real :: xl,yl,wx,wy,w
  real,parameter :: factor=.596831, hxmax=6.0, hymax=4.0, hzmax=150.
!  integer xscav_count

  ! For forward simulations, make a loop over the number of species;
  ! for backward simulations, make an additional loop over the
  ! releasepoints
  !***************************************************************************
!  xscav_count=0
  do i=1,numpart
    if (itra1(i).ne.itime) goto 20

  ! Determine age class of the particle
    itage=abs(itra1(i)-itramem(i))
    do nage=1,nageclass
      if (itage.lt.lage(nage)) goto 33
    end do
33   continue

!  if (xscav_frac1(i,1).lt.0) xscav_count=xscav_count+1
           
  ! For special runs, interpolate the air density to the particle position
  !************************************************************************
  !***********************************************************************
  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af    NOTE that in backward simulations the release of particles takes place
  !Af    at the receptor and the sampling at the source.
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1="mass"
  !Af          2="mass mixing ratio"

  !Af switches for the conccalcfile:
  !AF IND_SAMP =  0 : xmass * 1
  !Af IND_SAMP = -1 : xmass / rho

  !Af ind_samp is defined in readcommand.f

    if ( ind_samp .eq. -1 ) then

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

! eso: Temporary fix for particle exactly at north pole
      if (jyp >= nymax) then
      !  write(*,*) 'WARNING: conccalc.f90 jyp >= nymax'
        jyp=jyp-1
      end if

      do il=2,nz
        if (height(il).gt.ztra1(i)) then
          indz=il-1
          indzp=il
          goto 6
        endif
      end do
6     continue

      dz1=ztra1(i)-height(indz)
      dz2=height(indzp)-ztra1(i)
      dz=1./(dz1+dz2)

  ! Take density from 2nd wind field in memory (accurate enough, no time interpolation needed)
  !*****************************************************************************
      do ind=indz,indzp
        rhoprof(ind-indz+1)=p1*rho(ix ,jy ,ind,memind(2)) &
             +p2*rho(ixp,jy ,ind,2) &
             +p3*rho(ix ,jyp,ind,2) &
             +p4*rho(ixp,jyp,ind,2)
      end do
      rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz
   elseif (ind_samp.eq.0) then
      rhoi = 1.
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


  !********************************
  ! Do everything for mother domain
  !********************************

      xl=(xtra1(i)*dx+xoutshift)/dxout
      yl=(ytra1(i)*dy+youtshift)/dyout
      ix=int(xl)
      if (xl.lt.0.) ix=ix-1
      jy=int(yl)
      if (yl.lt.0.) jy=jy-1



  ! For particles aged less than 3 hours, attribute particle mass to grid cell
  ! it resides in rather than use the kernel, in order to avoid its smoothing effect.
  ! For older particles, use the uniform kernel.
  ! If a particle is close to the domain boundary, do not use the kernel either.
  !*****************************************************************************

      if ((.not.lusekerneloutput).or.(itage.lt.10800).or. &
           (xl.lt.0.5).or.(yl.lt.0.5).or. &
           (xl.gt.real(numxgrid-1)-0.5).or. &
           (yl.gt.real(numygrid-1)-0.5)) then             ! no kernel, direct attribution to grid cell
        if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
             (jy.le.numygrid-1)) then
          if (DRYBKDEP.or.WETBKDEP) then
            do ks=1,nspec
              gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
            end do
          else
            if (lparticlecountoutput) then
              do ks=1,nspec
                gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+1
              end do
            else
              do ks=1,nspec
                gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight
              end do
            end if
          endif
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
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w
               end do
            endif
          endif

          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=wx*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
              do ks=1,nspec
                 gridunc(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
               end do
             else
              do ks=1,nspec
                 gridunc(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w
               end do
             endif
          endif
        endif !ix ge 0


        if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=(1.-wx)*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w
               end do
            endif
          endif

          if ((jy.ge.0).and.(jy.le.numygrid-1)) then
            w=(1.-wx)*wy
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*w
               end do
            endif
          endif
        endif !ixp ge 0
     endif

  !************************************
  ! Do everything for the nested domain
  !************************************

      if (nested_output.eq.1) then
        xl=(xtra1(i)*dx+xoutshiftn)/dxoutn
        yl=(ytra1(i)*dy+youtshiftn)/dyoutn
        ix=int(xl)
        if (xl.lt.0.) ix=ix-1
        jy=int(yl)
        if (yl.lt.0.) jy=jy-1


  ! For particles aged less than 3 hours, attribute particle mass to grid cell
  ! it resides in rather than use the kernel, in order to avoid its smoothing effect.
  ! For older particles, use the uniform kernel.
  ! If a particle is close to the domain boundary, do not use the kernel either.
  !*****************************************************************************

        if ((itage.lt.10800).or.(xl.lt.0.5).or.(yl.lt.0.5).or. &
             (xl.gt.real(numxgridn-1)-0.5).or. &
             (yl.gt.real(numygridn-1)-0.5).or.((.not.lusekerneloutput))) then
! no kernel, direct attribution to grid cell
          if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
               (jy.le.numygridn-1)) then
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                   griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                   xmass1(i,ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
              if (lparticlecountoutput) then
                do ks=1,nspec
                  griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+1
                end do
              else
                do ks=1,nspec
                  griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                       xmass1(i,ks)/rhoi*weight
                end do
              endif
            endif
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

          if ((ix.ge.0).and.(ix.le.numxgridn-1)) then
            if ((jy.ge.0).and.(jy.le.numygridn-1)) then
              w=wx*wy
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                do ks=1,nspec
                   griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w
                 end do
              endif
            endif

            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=wx*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                   griduncn(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w
                 end do
              endif
            endif
          endif


          if ((ixp.ge.0).and.(ixp.le.numxgridn-1)) then
            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=(1.-wx)*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                   griduncn(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w
                 end do
              endif
            endif

            if ((jy.ge.0).and.(jy.le.numygridn-1)) then
              w=(1.-wx)*wy
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                    griduncn(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,nclass(i),nage)+ &
                     xmass1(i,ks)/rhoi*weight*w
                 end do
              endif
            endif
          endif
        endif
      endif
    endif
20  continue
  end do
!  write(*,*) 'xscav count:',xscav_count

  !***********************************************************************
  ! 2. Evaluate concentrations at receptor points, using the kernel method
  !***********************************************************************

  do n=1,numreceptor


  ! Reset concentrations
  !*********************

    do ks=1,nspec
      c(ks)=0.
    end do


  ! Estimate concentration at receptor
  !***********************************

    do i=1,numpart

      if (itra1(i).ne.itime) goto 40
      itage=abs(itra1(i)-itramem(i))

      hz=min(50.+0.3*sqrt(real(itage)),hzmax)
      zd=ztra1(i)/hz
      if (zd.gt.1.) goto 40          ! save computing time, leave loop

      hx=min((0.29+2.222e-3*sqrt(real(itage)))*dx+ &
           real(itage)*1.2e-5,hxmax)                     ! 80 km/day
      xd=(xtra1(i)-xreceptor(n))/hx
      if (xd*xd.gt.1.) goto 40       ! save computing time, leave loop

      hy=min((0.18+1.389e-3*sqrt(real(itage)))*dy+ &
           real(itage)*7.5e-6,hymax)                     ! 80 km/day
      yd=(ytra1(i)-yreceptor(n))/hy
      if (yd*yd.gt.1.) goto 40       ! save computing time, leave loop
      h=hx*hy*hz

      r2=xd*xd+yd*yd+zd*zd
      if (r2.lt.1.) then
        xkern=factor*(1.-r2)
        do ks=1,nspec
          c(ks)=c(ks)+xmass1(i,ks)*xkern/h
        end do
      endif
40    continue
    end do

    do ks=1,nspec
      creceptor(n,ks)=creceptor(n,ks)+2.*weight*c(ks)/receptorarea(n)
    end do
  end do

end subroutine conccalc
