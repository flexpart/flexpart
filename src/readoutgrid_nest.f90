subroutine readoutgrid_nest

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the output nest.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     4 June 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dxoutn,dyoutn        grid distances of output nest                         *
  ! numxgridn,numygridn,numzgrid    nest dimensions                            *
  ! outlon0n,outlat0n    lower left corner of nest                             *
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

  integer :: stat
  real :: xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4

  integer :: readerror

  ! declare namelist
  namelist /outgridn/ &
    outlon0n,outlat0n, &
    numxgridn,numygridn, &
    dxoutn,dyoutn

  ! helps identifying failed namelist input
  dxoutn=-1.0

  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID_NEST',form='formatted',status='old',err=999)

  ! try namelist input
  read(unitoutgrid,outgridn,iostat=readerror)
  close(unitoutgrid)

  if ((dxoutn.le.0).or.(readerror.ne.0)) then

    open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID_NEST',status='old',err=999)
    call skplin(5,unitoutgrid)

    ! 1.  Read horizontal grid specifications
    !****************************************

    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlon0n
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlat0n
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numxgridn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numygridn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dxoutn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dyoutn

    close(unitoutgrid)
  endif

  ! write outgrid_nest file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitoutgrid,file=path(2)(1:length(2))//'OUTGRID_NEST.namelist',err=1000)
    write(unitoutgrid,nml=outgridn)
    close(unitoutgrid)
  endif

  allocate(orooutn(0:numxgridn-1,0:numygridn-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate orooutn'
  allocate(arean(0:numxgridn-1,0:numygridn-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate arean'
  allocate(volumen(0:numxgridn-1,0:numygridn-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate volumen'

  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0n+real(numxgridn)*dxoutn
  yr=outlat0n+real(numygridn)*dyoutn
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0n+eps.lt.xlon0).or.(outlat0n+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### NEST IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') path(1)(1:length(1))
    stop
  endif

  xoutshiftn=xlon0-outlon0n
  youtshiftn=ylat0-outlat0n
  return

999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop

end subroutine readoutgrid_nest
