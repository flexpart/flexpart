! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine readreceptors

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the receptor points.    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !     1 August 1996                                                          *
  !     HSO, 14 August 2013
  !     Added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! receptorarea(maxreceptor)  area of dx*dy at location of receptor           *
  ! receptorname(maxreceptor)  names of receptors                              *
  ! xreceptor,yreceptor  coordinates of receptor points                        *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitreceptor         unit connected to file RECEPTORS                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: j
  real :: x,y,xm,ym
  character(len=16) :: receptor

  integer :: readerror
  real :: lon,lat   ! for namelist input, lon/lat are used instead of x,y
  integer,parameter :: unitreceptorout=2

  ! declare namelist
  namelist /receptors/ &
    receptor, lon, lat

  lon=-999.9

  ! For backward runs, do not allow receptor output. Thus, set number of receptors to zero
  !*****************************************************************************

  if (ldirect.lt.0) then
    numreceptor=0
    return
  endif


  ! Open the RECEPTORS file and read output grid specifications
  !************************************************************

  open(unitreceptor,file=path(1)(1:length(1))//'RECEPTORS',form='formatted',status='old',err=999)

  ! try namelist input
  read(unitreceptor,receptors,iostat=readerror)

  ! prepare namelist output if requested
  if (nmlout.and.lroot) then
    open(unitreceptorout,file=path(2)(1:length(2))//'RECEPTORS.namelist',&
         &access='append',status='replace',err=1000)
  endif

  if ((lon.lt.-900).or.(readerror.ne.0)) then

    close(unitreceptor)
    open(unitreceptor,file=path(1)(1:length(1))//'RECEPTORS',status='old',err=999)
    call skplin(5,unitreceptor)

    ! Read the names and coordinates of the receptors
    !************************************************

    j=0
100 j=j+1
    read(unitreceptor,*,end=99)
    read(unitreceptor,*,end=99)
    read(unitreceptor,*,end=99)
    read(unitreceptor,'(4x,a16)',end=99) receptor
    call skplin(3,unitreceptor)
    read(unitreceptor,'(4x,f11.4)',end=99) x
    call skplin(3,unitreceptor)
    read(unitreceptor,'(4x,f11.4)',end=99) y
    if ((x.eq.0.).and.(y.eq.0.).and. &
         (receptor.eq.'                ')) then
      j=j-1
      goto 100
    endif
    if (j.gt.maxreceptor) then
      write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY RECEPTOR #### '
      write(*,*) ' #### POINTS ARE GIVEN.                       #### '
      write(*,*) ' #### MAXIMUM NUMBER IS ',maxreceptor,'       #### '
      write(*,*) ' #### PLEASE MAKE CHANGES IN FILE RECEPTORS   #### '
    endif
    receptorname(j)=receptor
    xreceptor(j)=(x-xlon0)/dx       ! transform to grid coordinates
    yreceptor(j)=(y-ylat0)/dy
    xm=r_earth*cos(y*pi/180.)*dx/180.*pi
    ym=r_earth*dy/180.*pi
    receptorarea(j)=xm*ym

    ! write receptors file in namelist format to output directory if requested
    if (nmlout.and.lroot) then
      lon=x
      lat=y
      write(unitreceptorout,nml=receptors)
    endif

    goto 100

99  numreceptor=j-1

  else ! continue with namelist input

    j=0
    do while (readerror.eq.0)
      j=j+1
      lon=-999.9
      read(unitreceptor,receptors,iostat=readerror)
      if ((lon.lt.-900).or.(readerror.ne.0)) then
        readerror=1
      else
        if (j.gt.maxreceptor) then
          write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY RECEPTOR #### '
          write(*,*) ' #### POINTS ARE GIVEN.                       #### '
          write(*,*) ' #### MAXIMUM NUMBER IS ',maxreceptor,'       #### '
          write(*,*) ' #### PLEASE MAKE CHANGES IN FILE RECEPTORS   #### '
        endif
        receptorname(j)=receptor
        xreceptor(j)=(lon-xlon0)/dx       ! transform to grid coordinates
        yreceptor(j)=(lat-ylat0)/dy
        xm=r_earth*cos(lat*pi/180.)*dx/180.*pi
        ym=r_earth*dy/180.*pi
        receptorarea(j)=xm*ym
      endif

      ! write receptors file in namelist format to output directory if requested
      if (nmlout.and.lroot) then
        write(unitreceptorout,nml=receptors)
      endif

    end do
    numreceptor=j-1
    close(unitreceptor)

  endif

  if (nmlout.and.lroot) then
    close(unitreceptorout)
  endif

  return


999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "RECEPTORS"  #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "RECEPTORS"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  stop

end subroutine readreceptors
