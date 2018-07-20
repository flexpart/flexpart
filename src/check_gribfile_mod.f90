module check_gribfile_mod

!*
! Valid-License-Identifier:    	GPL-3.0-or-later
! Copyright (c) 2018 Petra Seibert (petra seibert at univie ac at)
!
! Prepared for use in FLEXPART (see flexpart.eu) version 10+
! 1. provide centre ids for ECMWF and NCEP
! 2. obtain centre id from a given grib file
! intended to replace class_gribfile from ctbto project
! requires grib_api_fortran90 and grib_api
! either from eccodes or grib_api

! I am using the convention to put an abbreviated module name (here: cg)
! in front of public entities. If they are integer, then icg.
! subroutines / functions from external libs with upper first letter

  implicit none

  integer, parameter :: icg_id_ncep = 7, icg_id_ecmwf = 98
!!     official centre codes 
  integer :: icentre ! centre ID found

contains

  subroutine cg_get_centre(pfname, icentre)
  
! this subroutine returns a code for the centre which produced the gribfile  

  use grib_api

  integer, intent(out) :: icentre
  integer              :: ifile ! grib file handle
  integer              :: iret  ! status code
  integer              :: igrib ! message handle
  character (len=*), intent(in)  :: pfname ! path+filenmae

  call Grib_open_file(ifile,pfname,'r',iret)
  if (iret == 0) then
  
    call Grib_new_from_file(ifile,igrib)    ! load first message
    call Grib_get(igrib,'centre',icentre)  ! read centre ID
    call Grib_close_file(ifile)              ! done

    if (icentre .ne. icg_id_ecmwf .and. &
        icentre .ne. icg_id_ncep) then
!! centre not foreseen in Fp
      write (*,*) ' #### FLEXPART MODEL ERROR! Met input file'
      write (*,*) ' #### '//trim(pfname)
      write (*,*) ' #### is neither ECMWF nor NCEP grib file'
      stop 'error in check_gribfile'
    endif

  else

!! reading gribfile failed
    write (*,*) ' #### FLEXPART MODEL ERROR! Met input file'
    write (*,*) ' #### '//trim(pfname)
    write (*,*) ' #### cannot be opened with grib_api'
    stop 'error in check_gribfile'

  endif

  end subroutine cg_get_centre
  
end module check_gribfile_mod

! what needs to be done to replace the current code with this one:
! see email
