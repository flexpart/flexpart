subroutine concoutput_inversion_nest(itime,outnum)
  !                        i     i
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !
  !     January 2017,  Separate files by release but include all timesteps     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use point_mod
  use outg_mod
  use par_mod
  use com_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  !real densityoutgrid(0:numxgrid-1,0:numygrid-1,numzgrid),
  !    +grid(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec,maxpointspec_act,
  !    +    maxageclass)
  !real wetgrid(0:numxgrid-1,0:numygrid-1,maxspec,maxpointspec_act,
  !    +       maxageclass)
  !real drygrid(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +       maxpointspec_act,maxageclass)
  !real gridsigma(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec,
  !    +       maxpointspec_act,maxageclass),
  !    +     drygridsigma(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +     maxpointspec_act,maxageclass),
  !    +     wetgridsigma(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +     maxpointspec_act,maxageclass)
  !real factor(0:numxgrid-1,0:numygrid-1,numzgrid)
  !real sparse_dump_r(numxgrid*numygrid*numzgrid)
  !integer sparse_dump_i(numxgrid*numygrid*numzgrid)

  !real sparse_dump_u(numxgrid*numygrid*numzgrid)
  real(dep_prec) :: auxgrid(nclassunc)
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  logical,save :: lnstart=.true.
  logical,save,allocatable,dimension(:) :: lnstartrel
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist
  character :: areldate*8,areltime*6

  if(lnstart) then
    allocate(lnstartrel(maxpointspec_act))
    lnstartrel(:)=.true.
  endif
  print*, 'lnstartrel = ',lnstartrel

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  print*, 'outnum:',outnum
  print*, 'datetime:',adate//atime

  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


    if (ldirect.eq.1) then
       do ks=1,nspec
         do kp=1,maxpointspec_act
           tot_mu(ks,kp)=1
         end do
       end do
   else
      do ks=1,nspec
             do kp=1,maxpointspec_act
               tot_mu(ks,kp)=xmass(kp,ks)
             end do
      end do
    endif


  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46   kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,2)*dz1+ &
             rho_dry(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,2)
! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,2)
  end do

! RLT
! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
    do kz=1,numzgrid
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          if (ldirect.eq.1) then
            factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
          else
            factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
          endif
        end do
      end do
    end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

  write(anspec,'(i3.3)') ks

    do kp=1,maxpointspec_act

      print*, 'itime = ',itime
      print*, 'lage(1) = ',lage(1)
      print*, 'ireleasestart(kp) = ',ireleasestart(kp)
      print*, 'ireleaseend(kp) = ',ireleaseend(kp)

      ! check itime is within release and backward trajectory length
      if (nageclass.eq.1) then
        if ((itime.gt.ireleaseend(kp)).or.(itime.lt.(ireleasestart(kp)-lage(1)))) then
          go to 10
        endif
      endif

      ! calculate date of release
      jul=bdate+real(ireleasestart(kp),kind=dp)/86400._dp    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)
      write(areldate,'(i8.8)') jjjjmmdd
      write(areltime,'(i6.6)') ihmmss
      print*, areldate//areltime

      ! calculate date of field
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
      print*, adate//atime

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
          ! concentrations
          inquire(file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lnstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
          endif
        else
          ! residence times
          inquire(file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lnstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
          endif
        endif
        write(unitoutgrid) jjjjmmdd
        write(unitoutgrid) ihmmss
      endif

      if ((iout.eq.2).or.(iout.eq.3)) then
        ! mixing ratio
        inquire(file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
                areltime//'_'//anspec,exist=lexist)
        if(lexist.and..not.lnstartrel(kp)) then
          ! open and append to existing file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
        else
          ! open new file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='replace',action='write')
        endif
        write(unitoutgridppt) jjjjmmdd
        write(unitoutgridppt) ihmmss
      endif

      lnstartrel(kp)=.false.

      do nage=1,nageclass

        do jy=0,numygridn-1
          do ix=0,numxgridn-1

!  ! WET DEPOSITION
!            if ((WETDEP).and.(ldirect.gt.0)) then
!              do l=1,nclassunc
!                auxgrid(l)=wetgriduncn(ix,jy,ks,kp,l,nage)
!              end do
!              call mean(auxgrid,wetgrid(ix,jy), &
!                   wetgridsigma(ix,jy),nclassunc)
!  ! Multiply by number of classes to get total concentration
!              wetgrid(ix,jy)=wetgrid(ix,jy) &
!                   *nclassunc
!  ! Calculate standard deviation of the mean
!              wetgridsigma(ix,jy)= &
!                   wetgridsigma(ix,jy)* &
!                   sqrt(real(nclassunc))
!            endif

!  ! DRY DEPOSITION
!            if ((DRYDEP).and.(ldirect.gt.0)) then
!              do l=1,nclassunc
!                auxgrid(l)=drygriduncn(ix,jy,ks,kp,l,nage)
!              end do
!              call mean(auxgrid,drygrid(ix,jy), &
!                   drygridsigma(ix,jy),nclassunc)
!  ! Multiply by number of classes to get total concentration
!              drygrid(ix,jy)=drygrid(ix,jy)* &
!                   nclassunc
!  ! Calculate standard deviation of the mean
!              drygridsigma(ix,jy)= &
!                   drygridsigma(ix,jy)* &
!                   sqrt(real(nclassunc))
!            endif

  ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                     gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*nclassunc
  ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
            end do
          end do
        end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Concentration output
  !*********************

        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

!  ! Wet deposition
!          sp_count_i=0
!          sp_count_r=0
!          sp_fact=-1.
!          sp_zer=.true.
!          if ((ldirect.eq.1).and.(WETDEP)) then
!          do jy=0,numygridn-1
!            do ix=0,numxgridn-1
!  ! concentration greater zero
!              if (wetgrid(ix,jy).gt.smallnum) then
!                 if (sp_zer.eqv..true.) then ! first non zero value
!                    sp_count_i=sp_count_i+1
!                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
!                    sp_zer=.false.
!                    sp_fact=sp_fact*(-1.)
!                 endif
!                 sp_count_r=sp_count_r+1
!                 sparse_dump_r(sp_count_r)= &
!                      sp_fact*1.e12*wetgrid(ix,jy)/arean(ix,jy)
!                 sparse_dump_u(sp_count_r)= &
!                      1.e12*wetgridsigma(ix,jy)/area(ix,jy)
!              else ! concentration is zero
!                  sp_zer=.true.
!              endif
!            end do
!         end do
!         else
!            sp_count_i=0
!            sp_count_r=0
!         endif
!         write(unitoutgrid) sp_count_i
!         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
!         write(unitoutgrid) sp_count_r
!         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
!         write(unitoutgrid) sp_count_r
!         write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)

!  ! Dry deposition
!         sp_count_i=0
!         sp_count_r=0
!         sp_fact=-1.
!         sp_zer=.true.
!         if ((ldirect.eq.1).and.(DRYDEP)) then
!          do jy=0,numygridn-1
!            do ix=0,numxgridn-1
!              if (drygrid(ix,jy).gt.smallnum) then
!                 if (sp_zer.eqv..true.) then ! first non zero value
!                    sp_count_i=sp_count_i+1
!                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
!                    sp_zer=.false.
!                    sp_fact=sp_fact*(-1.)
!                 endif
!                 sp_count_r=sp_count_r+1
!                 sparse_dump_r(sp_count_r)= &
!                      sp_fact* &
!                      1.e12*drygrid(ix,jy)/arean(ix,jy)
!                 sparse_dump_u(sp_count_r)= &
!                      1.e12*drygridsigma(ix,jy)/area(ix,jy)
!              else ! concentration is zero
!                  sp_zer=.true.
!              endif
!            end do
!          end do
!         else
!            sp_count_i=0
!            sp_count_r=0
!         endif
!         write(unitoutgrid) sp_count_i
!         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
!         write(unitoutgrid) sp_count_r
!         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
!         write(unitoutgrid) sp_count_r
!         write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)
!

  ! Concentrations

  ! surf_only write only 1st layer 

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
  !                 if ((factor(ix,jy,kz)/tot_mu(ks,kp)).eq.0)
  !    +              write (*,*) factor(ix,jy,kz),tot_mu(ks,kp),ks,kp
                   sparse_dump_u(sp_count_r)= &
                        gridsigma(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
!         write(unitoutgrid) sp_count_r
!         write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)

      endif !  concentration output

  ! Mixing ratio output
  !********************

      if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

!  ! Wet deposition
!         sp_count_i=0
!         sp_count_r=0
!         sp_fact=-1.
!         sp_zer=.true.
!         if ((ldirect.eq.1).and.(WETDEP)) then
!          do jy=0,numygridn-1
!            do ix=0,numxgridn-1
!                if (wetgrid(ix,jy).gt.smallnum) then
!                  if (sp_zer.eqv..true.) then ! first non zero value
!                    sp_count_i=sp_count_i+1
!                    sparse_dump_i(sp_count_i)= &
!                         ix+jy*numxgridn
!                    sp_zer=.false.
!                    sp_fact=sp_fact*(-1.)
!                 endif
!                 sp_count_r=sp_count_r+1
!                 sparse_dump_r(sp_count_r)= &
!                      sp_fact* &
!                      1.e12*wetgrid(ix,jy)/arean(ix,jy)
!                 sparse_dump_u(sp_count_r)= &
!                      1.e12*wetgridsigma(ix,jy)/area(ix,jy)
!              else ! concentration is zero
!                  sp_zer=.true.
!              endif
!            end do
!          end do
!         else
!           sp_count_i=0
!           sp_count_r=0
!         endif
!         write(unitoutgridppt) sp_count_i
!         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
!         write(unitoutgridppt) sp_count_r
!         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
!         write(unitoutgridppt) sp_count_r
!         write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)
!

!  ! Dry deposition
!         sp_count_i=0
!         sp_count_r=0
!         sp_fact=-1.
!         sp_zer=.true.
!         if ((ldirect.eq.1).and.(DRYDEP)) then
!          do jy=0,numygridn-1
!            do ix=0,numxgridn-1
!                if (drygrid(ix,jy).gt.smallnum) then
!                  if (sp_zer.eqv..true.) then ! first non zero value
!                    sp_count_i=sp_count_i+1
!                    sparse_dump_i(sp_count_i)= &
!                         ix+jy*numxgridn
!                    sp_zer=.false.
!                    sp_fact=sp_fact*(-1)
!                 endif
!                 sp_count_r=sp_count_r+1
!                 sparse_dump_r(sp_count_r)= &
!                      sp_fact* &
!                      1.e12*drygrid(ix,jy)/arean(ix,jy)
!                 sparse_dump_u(sp_count_r)= &
!                      1.e12*drygridsigma(ix,jy)/area(ix,jy)
!              else ! concentration is zero
!                  sp_zer=.true.
!              endif
!            end do
!          end do
!         else
!           sp_count_i=0
!           sp_count_r=0
!         endif
!         write(unitoutgridppt) sp_count_i
!         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
!         write(unitoutgridppt) sp_count_r
!         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
!         write(unitoutgridppt) sp_count_r
!         write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)
!

  ! Mixing ratios

    ! surf_only write only 1st layer 

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volumen(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*gridsigma(ix,jy,kz)/volumen(ix,jy,kz)/ &
                      outnum*weightair/weightmolar(ks)/ &
                      densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
!          write(unitoutgridppt) sp_count_r
!          write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)

        endif ! output for ppt
 
      end do ! nageclass

      close(unitoutgridppt)
      close(unitoutgrid)

      ! itime is outside range
10    continue

    end do ! maxpointspec_act

  end do ! nspec


! RLT Aug 2017
! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid_nest',exist=lexist)
  if (lexist.and..not.lnstart) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='replace',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,1
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgridn+kz*numxgridn*numygridn
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)

  ! reset lnstart
  if (lnstart) then
    lnstart=.false.
  endif

  ! Reinitialization of grid
  !*************************

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
      creceptor(i,ks)=0.
    end do
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        do l=1,nclassunc
          do nage=1,nageclass
            do kz=1,numzgrid
              griduncn(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do


end subroutine concoutput_inversion_nest

