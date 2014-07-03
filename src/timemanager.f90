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

subroutine timemanager

  !*****************************************************************************
  !                                                                            *
  ! Handles the computation of trajectories, i.e. determines which             *
  ! trajectories have to be computed at what time.                             *
  ! Manages dry+wet deposition routines, radioactive decay and the computation *
  ! of concentrations.                                                         *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 May 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !        Call of convmix when new windfield is read                          *
  !------------------------------------                                        *
  !  Changes Petra Seibert, Sept 2002                                          *
  !     fix wet scavenging problem                                             *
  !     Code may not be correct for decay of deposition!                       *
  !  Changes Petra Seibert, Nov 2002                                           *
  !     call convection BEFORE new fields are read in BWD mode                 *
  !  Changes Caroline Forster, Feb 2005                                        *
  !new interface between flexpart and convection scheme                        *
  !Emanuel's latest subroutine convect43c.f is used                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! DEP                .true. if either wet or dry deposition is switched on   *
  ! decay(maxspec) [1/s] decay constant for radioactive decay                  *
  ! DRYDEP             .true. if dry deposition is switched on                 *
  ! ideltas [s]        modelling period                                        *
  ! itime [s]          actual temporal position of calculation                 *
  ! ldeltat [s]        time since computation of radioact. decay of depositions*
  ! loutaver [s]       averaging period for concentration calculations         *
  ! loutend [s]        end of averaging for concentration calculations         *
  ! loutnext [s]       next time at which output fields shall be centered      *
  ! loutsample [s]     sampling interval for averaging of concentrations       *
  ! loutstart [s]      start of averaging for concentration calculations       *
  ! loutstep [s]       time interval for which concentrations shall be         *
  !                    calculated                                              *
  ! npoint(maxpart)    index, which starting point the trajectory has          *
  !                    starting positions of trajectories                      *
  ! nstop              serves as indicator for fate of particles               *
  !                    in the particle loop                                    *
  ! nstop1             serves as indicator for wind fields (see getfields)     *
  ! outnum             number of samples for each concentration calculation    *
  ! outnum             number of samples for each concentration calculation    *
  ! prob               probability of absorption at ground due to dry          *
  !                    deposition                                              *
  ! WETDEP             .true. if wet deposition is switched on                 *
  ! weight             weight for each concentration sample (1/2 or 1)         *
  ! uap(maxpart),ucp(maxpart),uzp(maxpart) = random velocities due to          *
  !                    turbulence                                              *
  ! us(maxpart),vs(maxpart),ws(maxpart) = random velocities due to inter-      *
  !                    polation                                                *
  ! xtra1(maxpart), ytra1(maxpart), ztra1(maxpart) =                           *
  !                    spatial positions of trajectories                       *
  !                                                                            *
  ! Constants:                                                                 *
  ! maxpart            maximum number of trajectories                          *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use point_mod
  use xmass_mod
  use flux_mod
  use outg_mod
  use oh_mod
  use par_mod
  use com_mod

  implicit none

  integer :: j,ks,kp,l,n,itime,nstop,nstop1
! integer :: ksp
  integer :: loutnext,loutstart,loutend
  integer :: ix,jy,ldeltat,itage,nage
  real :: outnum,weight,prob(maxspec)
  real :: uap(maxpart),ucp(maxpart),uzp(maxpart),decfact
  real :: us(maxpart),vs(maxpart),ws(maxpart)
  integer(kind=2) :: cbt(maxpart)
  real :: drydeposit(maxspec),gridtotalunc,wetgridtotalunc
  real :: drygridtotalunc,xold,yold,zold,xmassfract
  !double precision xm(maxspec,maxpointspec_act),
  !    +                 xm_depw(maxspec,maxpointspec_act),
  !    +                 xm_depd(maxspec,maxpointspec_act)


  !open(88,file='TEST.dat')

  ! First output for time 0
  !************************

  loutnext=loutstep/2
  outnum=0.
  loutstart=loutnext-loutaver/2
  loutend=loutnext+loutaver/2

  !  open(127,file=path(2)(1:length(2))//'depostat.dat'
  !    +  ,form='unformatted')
  !write (*,*) 'writing deposition statistics depostat.dat!'

  !**********************************************************************
  ! Loop over the whole modelling period in time steps of mintime seconds
  !**********************************************************************


  !write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,drygridtotalunc
  write(*,46) float(itime)/3600,itime,numpart
  if (verbosity.gt.0) then
    write (*,*) 'timemanager> starting simulation'
    if (verbosity.gt.1) then
      CALL SYSTEM_CLOCK(count_clock)
      WRITE(*,*) 'timemanager> SYSTEM CLOCK',(count_clock - count_clock0)/real(count_rate)
    endif     
  endif

  do itime=0,ideltas,lsynctime


  ! Computation of wet deposition, OH reaction and mass transfer
  ! between two species every lsynctime seconds
  ! maybe wet depo frequency can be relaxed later but better be on safe side
  ! wetdepo must be called BEFORE new fields are read in but should not
  ! be called in the very beginning before any fields are loaded, or
  ! before particles are in the system
  ! Code may not be correct for decay of deposition
  ! changed by Petra Seibert 9/02
  !********************************************************************

    if (WETDEP .and. itime .ne. 0 .and. numpart .gt. 0) then
        if (verbosity.gt.0) then
           write (*,*) 'timemanager> call wetdepo'
        endif     
         call wetdepo(itime,lsynctime,loutnext)
    endif

    if (OHREA .and. itime .ne. 0 .and. numpart .gt. 0) &
         call ohreaction(itime,lsynctime,loutnext)

    if (ASSSPEC .and. itime .ne. 0 .and. numpart .gt. 0) then
       stop 'associated species not yet implemented!'
  !     call transferspec(itime,lsynctime,loutnext)
    endif

  ! compute convection for backward runs
  !*************************************

   if ((ldirect.eq.-1).and.(lconvection.eq.1).and.(itime.lt.0)) then
        if (verbosity.gt.0) then
           write (*,*) 'timemanager> call convmix -- backward'
        endif         
      call convmix(itime)
        if (verbosity.gt.1) then
          !CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
          CALL SYSTEM_CLOCK(count_clock)
          WRITE(*,*) 'timemanager> SYSTEM CLOCK',(count_clock - count_clock0)/real(count_rate)
        endif 
   endif

  ! Get necessary wind fields if not available
  !*******************************************
    if (verbosity.gt.0) then
           write (*,*) 'timemanager> call getfields'
    endif 
    call getfields(itime,nstop1)
        if (verbosity.gt.1) then
          CALL SYSTEM_CLOCK(count_clock)
          WRITE(*,*) 'timemanager> SYSTEM CLOCK',(count_clock - count_clock0)/real(count_rate)
        endif 
    if (nstop1.gt.1) stop 'NO METEO FIELDS AVAILABLE'

  ! Release particles
  !******************

    if (verbosity.gt.0) then
           write (*,*) 'timemanager>  Release particles'
    endif 

    if (mdomainfill.ge.1) then
      if (itime.eq.0) then
        if (verbosity.gt.0) then
          write (*,*) 'timemanager>  call init_domainfill'
        endif       
        call init_domainfill
      else
        if (verbosity.gt.0) then
          write (*,*) 'timemanager>  call boundcond_domainfill'
        endif   
        call boundcond_domainfill(itime,loutend)
      endif
    else
      if (verbosity.gt.0) then
        print*,'call releaseparticles'  
      endif
      call releaseparticles(itime)
      if (verbosity.gt.1) then
        CALL SYSTEM_CLOCK(count_clock)
        WRITE(*,*) 'timemanager> SYSTEM CLOCK',(count_clock - count_clock0)/real(count_rate)
      endif 
    endif


  ! Compute convective mixing for forward runs
  ! for backward runs it is done before next windfield is read in
  !**************************************************************

   if ((ldirect.eq.1).and.(lconvection.eq.1)) then
     if (verbosity.gt.0) then
       write (*,*) 'timemanager> call convmix -- forward'
     endif    
     call convmix(itime)
   endif

  ! If middle of averaging period of output fields is reached, accumulated
  ! deposited mass radioactively decays
  !***********************************************************************

    if (DEP.and.(itime.eq.loutnext).and.(ldirect.gt.0)) then
      do ks=1,nspec
      do kp=1,maxpointspec_act
        if (decay(ks).gt.0.) then
          do nage=1,nageclass
            do l=1,nclassunc
  ! Mother output grid
              do jy=0,numygrid-1
                do ix=0,numxgrid-1
                  wetgridunc(ix,jy,ks,kp,l,nage)= &
                       wetgridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                  drygridunc(ix,jy,ks,kp,l,nage)= &
                       drygridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                end do
              end do
  ! Nested output grid
              if (nested_output.eq.1) then
                do jy=0,numygridn-1
                  do ix=0,numxgridn-1
                    wetgriduncn(ix,jy,ks,kp,l,nage)= &
                         wetgriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                    drygriduncn(ix,jy,ks,kp,l,nage)= &
                         drygriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                  end do
                end do
              endif
            end do
          end do
        endif
      end do
      end do
    endif

  !!! CHANGE: These lines may be switched on to check the conservation
  !!! of mass within FLEXPART
  !   if (itime.eq.loutnext) then
  !   do 247 ksp=1, nspec
  !   do 247 kp=1, maxpointspec_act
  !47         xm(ksp,kp)=0.

  !   do 249 ksp=1, nspec
  !     do 249 j=1,numpart
  !          if (ioutputforeachrelease.eq.1) then
  !            kp=npoint(j)
  !          else
  !            kp=1
  !          endif
  !       if (itra1(j).eq.itime) then
  !          xm(ksp,kp)=xm(ksp,kp)+xmass1(j,ksp)
  !         write(*,*) 'xmass: ',xmass1(j,ksp),j,ksp,nspec
  !       endif
  !49     continue
  !  do 248 ksp=1,nspec
  !  do 248 kp=1,maxpointspec_act
  !  xm_depw(ksp,kp)=0.
  !  xm_depd(ksp,kp)=0.
  !     do 248 nage=1,nageclass
  !       do 248 ix=0,numxgrid-1
  !         do 248 jy=0,numygrid-1
  !           do 248 l=1,nclassunc
  !              xm_depw(ksp,kp)=xm_depw(ksp,kp)
  !    +                  +wetgridunc(ix,jy,ksp,kp,l,nage)
  !48                 xm_depd(ksp,kp)=xm_depd(ksp,kp)
  !    +                  +drygridunc(ix,jy,ksp,kp,l,nage)
  !             do 246 ksp=1,nspec
  !46                    write(88,'(2i10,3e12.3)')
  !    +              itime,ksp,(xm(ksp,kp),kp=1,maxpointspec_act),
  !    +                (xm_depw(ksp,kp),kp=1,maxpointspec_act),
  !    +                (xm_depd(ksp,kp),kp=1,maxpointspec_act)
  !  endif
  !!! CHANGE



  ! Check whether concentrations are to be calculated
  !**************************************************

    if ((ldirect*itime.ge.ldirect*loutstart).and. &
         (ldirect*itime.le.ldirect*loutend)) then ! add to grid
      if (mod(itime-loutstart,loutsample).eq.0) then

  ! If we are exactly at the start or end of the concentration averaging interval,
  ! give only half the weight to this sample
  !*****************************************************************************

        if ((itime.eq.loutstart).or.(itime.eq.loutend)) then
          weight=0.5
        else
          weight=1.0
        endif
        outnum=outnum+weight
        call conccalc(itime,weight)
      endif


      if ((mquasilag.eq.1).and.(itime.eq.(loutstart+loutend)/2)) &
           call partoutput_short(itime)    ! dump particle positions in extremely compressed format


  ! Output and reinitialization of grid
  ! If necessary, first sample of new grid is also taken
  !*****************************************************

      if ((itime.eq.loutend).and.(outnum.gt.0.)) then
        if ((iout.le.3.).or.(iout.eq.5)) then
          if (surf_only.ne.1) then 
          call concoutput(itime,outnum,gridtotalunc, &
               wetgridtotalunc,drygridtotalunc)
          else  
  if (verbosity.eq.1) then
     print*,'call concoutput_surf '
     CALL SYSTEM_CLOCK(count_clock)
     WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif
          call concoutput_surf(itime,outnum,gridtotalunc, &
               wetgridtotalunc,drygridtotalunc)
  if (verbosity.eq.1) then
     print*,'called concoutput_surf '
     CALL SYSTEM_CLOCK(count_clock)
     WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif
          endif

          if ((nested_output.eq.1).and.(surf_only.ne.1)) call concoutput_nest(itime,outnum)
          if ((nested_output.eq.1).and.(surf_only.eq.1)) call concoutput_surf_nest(itime,outnum)
          outnum=0.
        endif
        if ((iout.eq.4).or.(iout.eq.5)) call plumetraj(itime)
        if (iflux.eq.1) call fluxoutput(itime)
        !write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,drygridtotalunc
        write(*,46) float(itime)/3600,itime,numpart
45      format(i9,' SECONDS SIMULATED: ',i8, ' PARTICLES:    Uncertainty: ',3f7.3)
46      format(' Simulated ',f7.1,' hours (',i9,' s), ',i8, ' particles')
        if (ipout.ge.1) call partoutput(itime)    ! dump particle positions
        loutnext=loutnext+loutstep
        loutstart=loutnext-loutaver/2
        loutend=loutnext+loutaver/2
        if (itime.eq.loutstart) then
          weight=0.5
          outnum=outnum+weight
          call conccalc(itime,weight)
        endif


  ! Check, whether particles are to be split:
  ! If so, create new particles and attribute all information from the old
  ! particles also to the new ones; old and new particles both get half the
  ! mass of the old ones
  !************************************************************************

        if (ldirect*itime.ge.ldirect*itsplit) then
          n=numpart
          do j=1,numpart
            if (ldirect*itime.ge.ldirect*itrasplit(j)) then
              if (n.lt.maxpart) then
                n=n+1
                itrasplit(j)=2*(itrasplit(j)-itramem(j))+itramem(j)
                itrasplit(n)=itrasplit(j)
                itramem(n)=itramem(j)
                itra1(n)=itra1(j)
                idt(n)=idt(j)
                npoint(n)=npoint(j)
                nclass(n)=nclass(j)
                xtra1(n)=xtra1(j)
                ytra1(n)=ytra1(j)
                ztra1(n)=ztra1(j)
                uap(n)=uap(j)
                ucp(n)=ucp(j)
                uzp(n)=uzp(j)
                us(n)=us(j)
                vs(n)=vs(j)
                ws(n)=ws(j)
                cbt(n)=cbt(j)
                do ks=1,nspec
                  xmass1(j,ks)=xmass1(j,ks)/2.
                  xmass1(n,ks)=xmass1(j,ks)
                end do
              endif
            endif
          end do
          numpart=n
        endif
      endif
    endif


    if (itime.eq.ideltas) exit         ! almost finished

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

    if (itime.lt.loutnext) then
      ldeltat=itime-(loutnext-loutstep)
    else                                  ! first half of next interval
      ldeltat=itime-loutnext
    endif


  ! Loop over all particles
  !************************

    do j=1,numpart


  ! If integration step is due, do it
  !**********************************

      if (itra1(j).eq.itime) then

        if (ioutputforeachrelease.eq.1) then
            kp=npoint(j)
        else
            kp=1
        endif
  ! Determine age class of the particle
        itage=abs(itra1(j)-itramem(j))
        do nage=1,nageclass
          if (itage.lt.lage(nage)) exit
        end do

  ! Initialize newly released particle
  !***********************************

        if ((itramem(j).eq.itime).or.(itime.eq.0)) &
             call initialize(itime,idt(j),uap(j),ucp(j),uzp(j), &
             us(j),vs(j),ws(j),xtra1(j),ytra1(j),ztra1(j),cbt(j))

  ! Memorize particle positions
  !****************************

        xold=xtra1(j)
        yold=ytra1(j)
        zold=ztra1(j)

  ! Integrate Lagevin equation for lsynctime seconds
  !*************************************************

        call advance(itime,npoint(j),idt(j),uap(j),ucp(j),uzp(j), &
             us(j),vs(j),ws(j),nstop,xtra1(j),ytra1(j),ztra1(j),prob, &
             cbt(j))

  ! Calculate the gross fluxes across layer interfaces
  !***************************************************

        if (iflux.eq.1) call calcfluxes(nage,j,xold,yold,zold)


  ! Determine, when next time step is due
  ! If trajectory is terminated, mark it
  !**************************************

        if (nstop.gt.1) then
          if (linit_cond.ge.1) call initial_cond_calc(itime,j)
          itra1(j)=-999999999
        else
          itra1(j)=itime+lsynctime


  ! Dry deposition and radioactive decay for each species
  ! Also check maximum (of all species) of initial mass remaining on the particle;
  ! if it is below a threshold value, terminate particle
  !*****************************************************************************

          xmassfract=0.
          do ks=1,nspec
            if (decay(ks).gt.0.) then             ! radioactive decay
              decfact=exp(-real(abs(lsynctime))*decay(ks))
            else
              decfact=1.
            endif

            if (DRYDEPSPEC(ks)) then        ! dry deposition
              drydeposit(ks)=xmass1(j,ks)*prob(ks)*decfact
              xmass1(j,ks)=xmass1(j,ks)*(1.-prob(ks))*decfact
              if (decay(ks).gt.0.) then   ! correct for decay (see wetdepo)
                drydeposit(ks)=drydeposit(ks)* &
                     exp(real(abs(ldeltat))*decay(ks))
              endif
            else                           ! no dry deposition
              xmass1(j,ks)=xmass1(j,ks)*decfact
            endif


            if (mdomainfill.eq.0) then
              if (xmass(npoint(j),ks).gt.0.) &
                   xmassfract=max(xmassfract,real(npart(npoint(j)))* &
                   xmass1(j,ks)/xmass(npoint(j),ks))
            else
              xmassfract=1.
            endif
          end do

          if (xmassfract.lt.0.0001) then   ! terminate all particles carrying less mass
            itra1(j)=-999999999
          endif

  !        Sabine Eckhardt, June 2008
  !        don't create depofield for backward runs
          if (DRYDEP.AND.(ldirect.eq.1)) then
            call drydepokernel(nclass(j),drydeposit,real(xtra1(j)), &
                 real(ytra1(j)),nage,kp)
            if (nested_output.eq.1) call drydepokernel_nest( &
                 nclass(j),drydeposit,real(xtra1(j)),real(ytra1(j)), &
                 nage,kp)
          endif

  ! Terminate trajectories that are older than maximum allowed age
  !***************************************************************

          if (abs(itra1(j)-itramem(j)).ge.lage(nageclass)) then
            if (linit_cond.ge.1) &
                 call initial_cond_calc(itime+lsynctime,j)
            itra1(j)=-999999999
          endif
        endif

      endif

    end do

  end do


  ! Complete the calculation of initial conditions for particles not yet terminated
  !*****************************************************************************

  do j=1,numpart
    if (linit_cond.ge.1) call initial_cond_calc(itime,j)
  end do

  if (ipout.eq.2) call partoutput(itime)     ! dump particle positions

  if (linit_cond.ge.1) call initial_cond_output(itime)   ! dump initial cond. field

  close(104)

  ! De-allocate memory and end
  !***************************

  if (iflux.eq.1) then
      deallocate(flux)
  endif
  if (OHREA.eqv..TRUE.) then
      deallocate(OH_field,OH_field_height)
  endif
  if (ldirect.gt.0) then
  deallocate(drygridunc,wetgridunc)
  endif
  deallocate(gridunc)
  deallocate(xpoint1,xpoint2,ypoint1,ypoint2,zpoint1,zpoint2,xmass)
  deallocate(ireleasestart,ireleaseend,npart,kindz)
  deallocate(xmasssave)
  if (nested_output.eq.1) then
     deallocate(orooutn, arean, volumen)
     if (ldirect.gt.0) then
     deallocate(griduncn,drygriduncn,wetgriduncn)
     endif
  endif
  deallocate(outheight,outheighthalf)
  deallocate(oroout, area, volume)

end subroutine timemanager

