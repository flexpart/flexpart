subroutine fluxoutput(itime)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Output of the gridded fluxes.                                          *
  !     Eastward, westward, northward, southward, upward and downward gross    *
  !     fluxes are written to output file in either sparse matrix or grid dump *
  !     format, whichever is more efficient.                                   *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     04 April 2000                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncellse         number of cells with non-zero values for eastward fluxes   *
  ! sparsee         .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use flux_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,ix,jy,kz,k,nage,jjjjmmdd,ihmmss,kp,i
  integer :: ncellse(maxspec,maxageclass),ncellsw(maxspec,maxageclass)
  integer :: ncellss(maxspec,maxageclass),ncellsn(maxspec,maxageclass)
  integer :: ncellsu(maxspec,maxageclass),ncellsd(maxspec,maxageclass)
  logical :: sparsee(maxspec,maxageclass),sparsew(maxspec,maxageclass)
  logical :: sparses(maxspec,maxageclass),sparsen(maxspec,maxageclass)
  logical :: sparseu(maxspec,maxageclass),sparsed(maxspec,maxageclass)
  character :: adate*8,atime*6


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  open(unitflux,file=path(2)(1:length(2))//'grid_flux_'//adate// &
       atime,form='unformatted')

  !**************************************************************
  ! Check, whether output of full grid or sparse matrix format is
  ! more efficient in terms of storage space. This is checked for
  ! every species and for every age class
  !**************************************************************

  do k=1,nspec
    do nage=1,nageclass
      ncellse(k,nage)=0
      ncellsw(k,nage)=0
      ncellsn(k,nage)=0
      ncellss(k,nage)=0
      ncellsu(k,nage)=0
      ncellsd(k,nage)=0
    end do
  end do

  do k=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do kz=1,numzgrid
            if (flux(2,ix,jy,kz,k,kp,nage).gt.0) ncellse(k,nage)= &
                 ncellse(k,nage)+1
            if (flux(1,ix,jy,kz,k,kp,nage).gt.0) ncellsw(k,nage)= &
                 ncellsw(k,nage)+1
            if (flux(4,ix,jy,kz,k,kp,nage).gt.0) ncellsn(k,nage)= &
                 ncellsn(k,nage)+1
            if (flux(3,ix,jy,kz,k,kp,nage).gt.0) ncellss(k,nage)= &
                 ncellss(k,nage)+1
            if (flux(5,ix,jy,kz,k,kp,nage).gt.0) ncellsu(k,nage)= &
                 ncellsu(k,nage)+1
            if (flux(6,ix,jy,kz,k,kp,nage).gt.0) ncellsd(k,nage)= &
                 ncellsd(k,nage)+1
          end do
        end do
      end do
    end do
  end do
  end do

  ! Output in sparse matrix format more efficient, if less than
  ! 2/5 of all cells contains concentrations>0
  !************************************************************

  do k=1,nspec
    do nage=1,nageclass
      if (4*ncellse(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsee(k,nage)=.true.
      else
        sparsee(k,nage)=.false.
      endif
      if (4*ncellsw(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsew(k,nage)=.true.
      else
        sparsew(k,nage)=.false.
      endif
      if (4*ncellsn(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsen(k,nage)=.true.
      else
        sparsen(k,nage)=.false.
      endif
      if (4*ncellss(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparses(k,nage)=.true.
      else
        sparses(k,nage)=.false.
      endif
      if (4*ncellsu(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparseu(k,nage)=.true.
      else
        sparseu(k,nage)=.false.
      endif
      if (4*ncellsd(k,nage).lt.numxgrid*numygrid*numzgrid) then
        sparsed(k,nage)=.true.
      else
        sparsed(k,nage)=.false.
      endif
    end do
  end do



  ! Flux output: divide by area and time to get flux in ng/m2/s
  !************************************************************

  write(unitflux) itime
  do k=1,nspec
  do kp=1,maxpointspec_act
    do nage=1,nageclass

      if (sparsee(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(2,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(2,ix,jy,kz,k,kp,nage)/areaeast(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(2,ix,jy,kz,k,kp,nage)/ &
                 areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsew(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(1,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(1,ix,jy,kz,k,kp,nage)/areaeast(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(1,ix,jy,kz,k,kp,nage)/ &
                 areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparses(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(3,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(3,ix,jy,kz,k,kp,nage)/areanorth(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(3,ix,jy,kz,k,kp,nage)/ &
                 areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsen(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1 ! north
              if (flux(4,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(4,ix,jy,kz,k,kp,nage)/areanorth(ix,jy,kz)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(4,ix,jy,kz,k,kp,nage)/ &
                 areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparseu(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(5,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(5,ix,jy,kz,k,kp,nage)/area(ix,jy)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(5,ix,jy,kz,k,kp,nage)/ &
                 area(ix,jy)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

      if (sparsed(k,nage)) then
        write(unitflux) 1
        do kz=1,numzgrid
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              if (flux(6,ix,jy,kz,k,kp,nage).gt.0.) write(unitflux) &
                   ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12* &
                   flux(6,ix,jy,kz,k,kp,nage)/area(ix,jy)/outstep
            end do
          end do
        end do
        write(unitflux) -999,999.
      else
        write(unitflux) 2
        do kz=1,numzgrid
          do ix=0,numxgrid-1
            write(unitflux) (1.e12*flux(6,ix,jy,kz,k,kp,nage)/ &
                 area(ix,jy)/outstep,jy=0,numygrid-1)
          end do
        end do
      endif

    end do
  end do
  end do


  close(unitflux)


  ! Reinitialization of grid
  !*************************

  do k=1,nspec
  do kp=1,maxpointspec_act
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
          do kz=1,numzgrid
            do nage=1,nageclass
              do i=1,6
                flux(i,ix,jy,kz,k,kp,nage)=0.
              end do
            end do
          end do
      end do
    end do
  end do
  end do


end subroutine fluxoutput
