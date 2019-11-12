subroutine readoutgrid

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the output grid.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     4 June 1996                                                            *
  !     HSO, 1 July 2014
  !     Added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dxout,dyout          grid distance                                         *
  ! numxgrid,numygrid,numzgrid    grid dimensions                              *
  ! outlon0,outlat0      lower left corner of grid                             *
  ! outheight(maxzgrid)  height levels of output grid [m]                      *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitoutgrid          unit connected to file OUTGRID                        *
  !                                                                            *
  !*****************************************************************************

  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i,j,stat
  real :: outhelp,xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4

  ! namelist variables
  integer, parameter :: maxoutlev=500
  integer :: readerror
  real,allocatable, dimension (:) :: outheights

  ! declare namelist
  namelist /outgrid/ &
    outlon0,outlat0, &
    numxgrid,numygrid, &
    dxout,dyout, &
    outheights

  ! allocate large array for reading input
  allocate(outheights(maxoutlev),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheights'

  ! helps identifying failed namelist input
  dxout=-1.0
  outheights=-1.0

  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID',status='old',form='formatted',err=999)

  ! try namelist input
  read(unitoutgrid,outgrid,iostat=readerror)
  close(unitoutgrid)

  if ((dxout.le.0).or.(readerror.ne.0)) then

    readerror=1

    open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID',status='old',err=999)

    call skplin(5,unitoutgrid)

    ! 1.  Read horizontal grid specifications
    !****************************************

    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlon0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlat0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numxgrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numygrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dxout
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dyout

  endif

  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0+real(numxgrid)*dxout
  yr=outlat0+real(numygrid)*dyout
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0+eps.lt.xlon0).or.(outlat0+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) outlon0,outlat0
    write(*,*) xr1,yr1,xlon0,ylat0,xr,yr,dxout,dyout
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') path(1)(1:length(1))
    stop
  endif

  ! 2. Count Vertical levels of output grid
  !****************************************

  if (readerror.ne.0) then
    j=0
100 j=j+1
    do i=1,3
      read(unitoutgrid,*,end=99)
    end do
    read(unitoutgrid,'(4x,f7.1)',end=99) outhelp
    if (outhelp.eq.0.) goto 99
    goto 100
99  numzgrid=j-1
  else
    do i=1,maxoutlev
      if (outheights(i).lt.0) exit
    end do
    numzgrid=i-1
  end if

  allocate(outheight(numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheight'
  allocate(outheighthalf(numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheighthalf'

  ! 2. Vertical levels of output grid
  !**********************************

  if (readerror.ne.0) then

    rewind(unitoutgrid)
    call skplin(29,unitoutgrid)

    do j=1,numzgrid
      do i=1,3
        read(unitoutgrid,*)
      end do
      read(unitoutgrid,'(4x,f7.1)') outhelp
      outheight(j)=outhelp
      outheights(j)=outhelp
    end do
    close(unitoutgrid)

  else

    do j=1,numzgrid
      outheight(j)=outheights(j)
    end do

  endif

  ! write outgrid file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    ! reallocate outheights with actually required dimension for namelist writing
    deallocate(outheights)
    allocate(outheights(numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outheights'

    do j=1,numzgrid
      outheights(j)=outheight(j)
    end do

    open(unitoutgrid,file=path(2)(1:length(2))//'OUTGRID.namelist',err=1000)
    write(unitoutgrid,nml=outgrid)
    close(unitoutgrid)
  endif

  ! Check whether vertical levels are specified in ascending order
  !***************************************************************

  do j=2,numzgrid
    if (outheight(j).le.outheight(j-1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### '
    write(*,*) ' #### OF OUTPUT LEVELS IS CORRUPT AT LEVEL    #### '
    write(*,*) ' #### ',j,'                              #### '
    write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### '
    endif
  end do

  ! Determine the half levels, i.e. middle levels of the output grid
  !*****************************************************************

  outheighthalf(1)=outheight(1)/2.
  do j=2,numzgrid
    outheighthalf(j)=(outheight(j-1)+outheight(j))/2.
  end do

  xoutshift=xlon0-outlon0
  youtshift=ylat0-outlat0

  allocate(oroout(0:numxgrid-1,0:numygrid-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate oroout'
  allocate(area(0:numxgrid-1,0:numygrid-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate area'
  allocate(volume(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate volume'
  allocate(areaeast(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areaeast'
  allocate(areanorth(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areanorth'
  return


999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop

end subroutine readoutgrid
