module file_netcdf

! ==== Procedures for reading and writing netcdf 3/4 files, written by Sourish Basu (summer 2010) ====
!
! 1. Reading a netcdf file:
!     (a) Open the file with nc_id = nc_open(filename, mode)
!     (b) Get the group (if relevant, for netcdf-4) with grp_id = nc_get_group(nc_id, group_name)
!     (c) Get the variable with var = nc_read_var(nc_id, var_name) or var = nc_read_var(grp_id, var_name)
!     (d) Close the file with call nc_close(nc_id)
! 2. Writing to a netcdf file:
!     (a) Open the file with nc_id = nc_open(filename, mode)
!     (b) Create a group (if relevant, for netcdf-4) with grp_id = nc_create_group(nc_id, var_name)
!     (c) Set global attributes with call nc_set_attrs(nc_id, attr_name, attr_value) or call nc_set_attrs(grp_id, attr_name, attr_value)
!     (d) Create dimensions with call nc_create_dim(nc_id, dim_name, dim_length) or call nc_create_dim(nc_id, dim_name, dim_length)
!     (e) Create and fill a variable call nc_dump_var(nc_id, var_name, dim_list, var_value, attr_names, attr_values) or
!         nc_dump_var(grp_id, var_name, dim_list, var_value, attr_names, attr_values)
!     (f) Close the file with call nc_close(nc_id)
!
! All the nc_* procedures take an optional intent(out) integer argument io_status. If there is no error in the procedure,
! then io_status = 0 and the control returns to the calling procedure. If there is an error in the procedure and io_status
! is not present, nc_error() is called and the program stops. If there is an error and io_status is present, control
! is tranferred back to the calling procedure with io_status >= 1.

use netcdf

implicit none

private

public :: nc_open, nc_create_var, nc_create_dim, nc_set_var, nc_dump_var, nc_close, nc_set_attrs
public :: t_ncVar, assignment(=), nc_uni_var, nc_uni_att
public :: nc_read_var, nc_get_dim, nc_get_group, nc_create_group, nc_get_attr, nc_grp_exists
public :: nc_variables_deflate, nc_deflate_level

integer :: nc_status
! for some reason, on alarik (lunarc.lu.se), calling this variable "NF90_ENOGRP" creates a conflict (it is already defined in the netcdf module), but using it directly from the netcdf module makes the software crach. Renaming it works though ...
integer, parameter :: NF90_ENOGRP_LOCAL = -125
logical :: nc_variables_deflate = .false.
integer :: nc_shuffle_filter = 1
integer :: nc_deflate_level = 5

type t_ncVar
    integer                         :: nc_id
    integer                         :: var_id
    character(len=80)               :: var_name, var_type
    character(len=80), allocatable  :: attr_names(:)
    character(len=256), allocatable :: attr_values(:)
end type t_ncVar

type nc_uni_var
    real(8), allocatable    :: r8_1d(:), r8_2d(:,:), r8_3d(:,:,:), r8_4d(:,:,:,:), r8_5d(:,:,:,:,:), r8_6d(:,:,:,:,:,:), r8_7d(:,:,:,:,:,:,:)
    real(4), allocatable    :: r4_1d(:), r4_2d(:,:), r4_3d(:,:,:), r4_4d(:,:,:,:), r4_5d(:,:,:,:,:), r4_6d(:,:,:,:,:,:), r4_7d(:,:,:,:,:,:,:)
    integer(1), allocatable :: i1_1d(:), i1_2d(:,:), i1_3d(:,:,:), i1_4d(:,:,:,:), i1_5d(:,:,:,:,:), i1_6d(:,:,:,:,:,:), i1_7d(:,:,:,:,:,:,:)
    integer(2), allocatable :: i2_1d(:), i2_2d(:,:), i2_3d(:,:,:), i2_4d(:,:,:,:), i2_5d(:,:,:,:,:), i2_6d(:,:,:,:,:,:), i2_7d(:,:,:,:,:,:,:)
    integer(4), allocatable :: i4_1d(:), i4_2d(:,:), i4_3d(:,:,:), i4_4d(:,:,:,:), i4_5d(:,:,:,:,:), i4_6d(:,:,:,:,:,:), i4_7d(:,:,:,:,:,:,:)
end type nc_uni_var

type nc_uni_att
    character(len=256)      :: c
    integer(1), allocatable :: i1_1d(:)
    integer(2), allocatable :: i2_1d(:)
    integer(4), allocatable :: i4_1d(:)
    real(4), allocatable    :: r4_1d(:)
    real(8), allocatable    :: r8_1d(:)
end type nc_uni_att

interface nc_set_var
    module procedure nc_set_var_i1_1d
    module procedure nc_set_var_i1_2d
    module procedure nc_set_var_i1_3d
    module procedure nc_set_var_i1_4d
    module procedure nc_set_var_i1_5d
    module procedure nc_set_var_i1_6d
    module procedure nc_set_var_i1_7d
    module procedure nc_set_var_i2_1d
    module procedure nc_set_var_i2_2d
    module procedure nc_set_var_i2_3d
    module procedure nc_set_var_i2_4d
    module procedure nc_set_var_i2_5d
    module procedure nc_set_var_i2_6d
    module procedure nc_set_var_i2_7d
    module procedure nc_set_var_i4_1d
    module procedure nc_set_var_i4_2d
    module procedure nc_set_var_i4_3d
    module procedure nc_set_var_i4_4d
    module procedure nc_set_var_i4_5d
    module procedure nc_set_var_i4_6d
    module procedure nc_set_var_i4_7d
    module procedure nc_set_var_r4_1d
    module procedure nc_set_var_r4_2d
    module procedure nc_set_var_r4_3d
    module procedure nc_set_var_r4_4d
    module procedure nc_set_var_r4_5d
    module procedure nc_set_var_r4_6d
    module procedure nc_set_var_r4_7d
    module procedure nc_set_var_r8_1d
    module procedure nc_set_var_r8_2d
    module procedure nc_set_var_r8_3d
    module procedure nc_set_var_r8_4d
    module procedure nc_set_var_r8_5d
    module procedure nc_set_var_r8_6d
    module procedure nc_set_var_r8_7d
end interface nc_set_var

interface nc_dump_var
    module procedure nc_dump_var_i1_1d
    module procedure nc_dump_var_i1_2d
    module procedure nc_dump_var_i1_3d
    module procedure nc_dump_var_i1_4d
    module procedure nc_dump_var_i1_5d
    module procedure nc_dump_var_i1_6d
    module procedure nc_dump_var_i1_7d
    module procedure nc_dump_var_i2_1d
    module procedure nc_dump_var_i2_2d
    module procedure nc_dump_var_i2_3d
    module procedure nc_dump_var_i2_4d
    module procedure nc_dump_var_i2_5d
    module procedure nc_dump_var_i2_6d
    module procedure nc_dump_var_i2_7d
    module procedure nc_dump_var_i4_1d
    module procedure nc_dump_var_i4_2d
    module procedure nc_dump_var_i4_3d
    module procedure nc_dump_var_i4_4d
    module procedure nc_dump_var_i4_5d
    module procedure nc_dump_var_i4_6d
    module procedure nc_dump_var_i4_7d
    module procedure nc_dump_var_r4_1d
    module procedure nc_dump_var_r4_2d
    module procedure nc_dump_var_r4_3d
    module procedure nc_dump_var_r4_4d
    module procedure nc_dump_var_r4_5d
    module procedure nc_dump_var_r4_6d
    module procedure nc_dump_var_r4_7d
    module procedure nc_dump_var_r8_1d
    module procedure nc_dump_var_r8_2d
    module procedure nc_dump_var_r8_3d
    module procedure nc_dump_var_r8_4d
    module procedure nc_dump_var_r8_5d
    module procedure nc_dump_var_r8_6d
    module procedure nc_dump_var_r8_7d
end interface nc_dump_var

interface nc_set_attrs
    module procedure nc_set_attrs_var
    module procedure nc_set_attrs_glb
end interface nc_set_attrs

interface nc_get_attr
    module procedure nc_get_attr_glb
    module procedure nc_get_attr_var
    module procedure nc_get_attr_varid
end interface nc_get_attr

! nc_set_attrs() is called to set file, group or variable attributes. To set file or group
! attributes, call it as:
!   call nc_set_attrs(nc_id, attr_name, attr_value)
! where nc_id can be a group or file ID, and attr_value must be a string. To set variable
! attributes, call:
!   call nc_set_attrs(var_id, attr_names, attr_values)
! where attr_names and attr_values are arrays of strings. There is no good reason for these
! to be arrays, and perhaps I will change it later to be consistent with the call to set
! a global attribute.

interface assignment (=)
    module procedure nc_var_to_array_i1_1d
    module procedure nc_var_to_array_i1_2d
    module procedure nc_var_to_array_i1_3d
    module procedure nc_var_to_array_i1_4d
    module procedure nc_var_to_array_i1_5d
    module procedure nc_var_to_array_i1_6d
    module procedure nc_var_to_array_i1_7d
    module procedure nc_var_to_array_i2_1d
    module procedure nc_var_to_array_i2_2d
    module procedure nc_var_to_array_i2_3d
    module procedure nc_var_to_array_i2_4d
    module procedure nc_var_to_array_i2_5d
    module procedure nc_var_to_array_i2_6d
    module procedure nc_var_to_array_i2_7d
    module procedure nc_var_to_array_i4_1d
    module procedure nc_var_to_array_i4_2d
    module procedure nc_var_to_array_i4_3d
    module procedure nc_var_to_array_i4_4d
    module procedure nc_var_to_array_i4_5d
    module procedure nc_var_to_array_i4_6d
    module procedure nc_var_to_array_i4_7d
    module procedure nc_var_to_array_r4_1d
    module procedure nc_var_to_array_r4_2d
    module procedure nc_var_to_array_r4_3d
    module procedure nc_var_to_array_r4_4d
    module procedure nc_var_to_array_r4_5d
    module procedure nc_var_to_array_r4_6d
    module procedure nc_var_to_array_r4_7d
    module procedure nc_var_to_array_r8_1d
    module procedure nc_var_to_array_r8_2d
    module procedure nc_var_to_array_r8_3d
    module procedure nc_var_to_array_r8_4d
    module procedure nc_var_to_array_r8_5d
    module procedure nc_var_to_array_r8_6d
    module procedure nc_var_to_array_r8_7d
    module procedure nc_attr_to_scalar_c
    module procedure nc_attr_to_scalar_i1
    module procedure nc_attr_to_scalar_i2
    module procedure nc_attr_to_scalar_i4
    module procedure nc_attr_to_scalar_r4
    module procedure nc_attr_to_scalar_r8
    module procedure nc_attr_to_array_i1
    module procedure nc_attr_to_array_i2
    module procedure nc_attr_to_array_i4
    module procedure nc_attr_to_array_r4
    module procedure nc_attr_to_array_r8
end interface

contains

subroutine nc_error(msg)
    ! Subroutine for printing an error message and stopping the program.
    ! The most common use is:
    ! if (nc_status /= NF90_NOERR) call nc_error(nf90_strerror(nc_status))
    character(len=*), intent(in) :: msg

    write(0,*) trim(adjustl(msg))
    stop

end subroutine nc_error

function nc_get_group(nc_id, group_name, io_status) result (group_id)
    ! Gets a group_id from an open nc_id and a group name.
    character(len=*), intent(in)   :: group_name
    integer, intent(in)            :: nc_id
    integer, optional, intent(out) :: io_status
    integer                        :: group_id, dummy_len, dummy_int
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_inq_ncid(nc_id, group_name, group_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Failure reading group ", a, " from parent ", a)') trim(group_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end function nc_get_group

function nc_grp_exists(nc_id, group_name, io_status) result (grp_exists)
    ! returns a boolean depending on whether a group exists
    ! or not in the ncid
    character(len=*), intent(in)   :: group_name
    integer, intent(in)            :: nc_id
    integer, optional, intent(out) :: io_status
    logical                        :: grp_exists
    integer                        :: dummy_id, dummy_len, dummy_int
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_inq_ncid(nc_id, group_name, dummy_id)
    if (nc_status == NF90_NOERR) then
        grp_exists = .true.
    else if (nc_status == NF90_ENOGRP_LOCAL) then
        grp_exists = .false.
    else
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not check if group ", a, " exists in parent ", a)') trim(group_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
    end if

end function nc_grp_exists

function nc_create_group(nc_id, group_name, io_status) result (group_id)
    ! Creates a group with name group_name in the file/group nc_id, and
    ! returns an integer group_id
    character(len=*), intent(in)   :: group_name
    integer, intent(in)            :: nc_id
    integer, optional, intent(out) :: io_status
    integer                        :: group_id, dummy_len, dummy_int
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_redef(nc_id)
    nc_status = nf90_def_grp(nc_id, trim(adjustl(group_name)), group_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not create group ", a, " in parent ", a)') trim(group_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end function nc_create_group

function nc_open(filename, mode, io_status) result (output_ncid)
    ! Opens a netcdf-3/4 file. Valid modes are:
    ! 'r' : read-only
    ! 'w' : create a netcdf-3 file for writing, overwrite existing file
    ! 'c' : create a netcdf-4 file for writing, overwrite existing file
    ! 'a' : read-write (append), netcdf-3 only (I think!)
    character(len=*), intent(in)   :: filename
    character(len=1), intent(in)   :: mode
    integer, optional, intent(out) :: io_status
    integer                        :: output_ncid

    if (present(io_status)) io_status = 0
    select case (mode)
        case ('r')
            nc_status = nf90_open(filename, NF90_NOWRITE, output_ncid)
        case ('w')
            nc_status = nf90_create(filename, NF90_CLOBBER, output_ncid)
        case ('c')
            nc_status = nf90_create(filename, or(NF90_CLOBBER,NF90_NETCDF4), output_ncid)
        case ('a')
            nc_status = nf90_open(filename, NF90_WRITE, output_ncid)
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            write(0,'("Could not open file ", a, " in mode ", a)') trim(filename), mode
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end function nc_open

subroutine nc_close(nc_id, io_status)
    ! Close a netcdf file given its file_ID
    integer, intent(in)            :: nc_id
    integer, optional, intent(out) :: io_status

    if (present(io_status)) io_status = 0
    nc_status = nf90_close(nc_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            write(0,*) 'Error in nc_close, almost unthinkable!'
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end subroutine nc_close

function nc_get_dim(nc_id, dim_name, io_status) result (dim_len)
    ! Given a file/group ID and a dimension name, get the integer
    ! dimension ID
    character(len=*), intent(in)   :: dim_name
    integer, intent(in)            :: nc_id
    integer, optional, intent(out) :: io_status
    integer                        :: dim_len, dim_id, dummy_len, dummy_int
    character(len=NF90_MAX_NAME)   :: useless_name
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_inq_dimid(nc_id, dim_name, dim_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not find dimension ", a, " in group ", a)') dim_name, dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
        return
    end if
    nc_status = nf90_inquire_dimension(nc_id, dim_id, useless_name, dim_len)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not get length of dimension ", a, " from group ", a)') dim_name, dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end function nc_get_dim

subroutine nc_create_dim(nc_id, dim_name, dim_length, io_status)
    ! Create a dimension with given name and length inside a
    ! file/group specified by nc_id
    character(len=*), intent(in)   :: dim_name
    integer, intent(in)            :: nc_id, dim_length
    integer, optional, intent(out) :: io_status
    integer                        :: output_dimid, dummy_len, dummy_int
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_redef(nc_id)
    if (dim_length /= 0) then
        nc_status = nf90_def_dim(nc_id, trim(adjustl(dim_name)), dim_length, output_dimid)
    else
        nc_status = nf90_def_dim(nc_id, trim(adjustl(dim_name)), NF90_UNLIMITED, output_dimid)
    end if
    nc_status = nf90_enddef(nc_id)

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not create dimension ", a, " in group ", a)') trim(adjustl(dim_name)), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end subroutine nc_create_dim

function nc_create_var(nc_id, var_name, dim_list, var_type, io_status) result (output_varid)
    ! Create a variable of type var_type with name var_name and dimensions dim_list
    ! inside a file/group given by nc_id. This routine is called internally by
    ! nc_dump_var, and is made public for people who want a more flexible option (than
    ! nc_dump_var) to write a variable. Type specification is one of 'byte', 'short',
    ! 'int', 'float' or 'double'. Upon successful creation, an integer variable ID is
    ! returned.
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    integer, allocatable           :: dim_id_list(:)
    character(len=*), intent(in)   :: dim_list(:)
    integer                        :: i, j, parent_id, group_id, dummy_len, dummy_int
    character(len=*), intent(in)   :: var_type
    type(t_ncVar)                  :: output_varid
    integer, optional, intent(out) :: io_status
    character(len=512)             :: dummy_name

    allocate(dim_id_list(size(dim_list)))
    nc_status = nf90_redef(nc_id)

    do i = 1, size(dim_list)
        group_id = nc_id
        do while (nf90_inq_dimid(group_id, trim(adjustl(dim_list(i))), j) == NF90_EBADDIM)
            nc_status = nf90_inq_grp_parent(group_id, parent_id)
            if (nc_status == NF90_ENOGRP_LOCAL) call nc_error('Root group reached, dimension '//trim(adjustl(dim_list(i)))//' still not found')
            group_id = parent_id
        end do
        dim_id_list(i) = j
    end do

    if (present(io_status)) io_status = 0
    select case (var_type)
        case ('byte')
            nc_status = nf90_def_var(nc_id, trim(adjustl(var_name)), NF90_BYTE, dim_id_list, j)
        case ('short')
            nc_status = nf90_def_var(nc_id, trim(adjustl(var_name)), NF90_SHORT, dim_id_list, j)
        case ('int')
            nc_status = nf90_def_var(nc_id, trim(adjustl(var_name)), NF90_INT, dim_id_list, j)
        case ('float')
            nc_status = nf90_def_var(nc_id, trim(adjustl(var_name)), NF90_FLOAT, dim_id_list, j)
        case ('double')
            nc_status = nf90_def_var(nc_id, trim(adjustl(var_name)), NF90_DOUBLE, dim_id_list, j)
        case default
            write(*,*) 'Error from nc_create_var : ', trim(nf90_strerror(nc_status))
            stop
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error creating variable ", a, " in group ", a)') trim(adjustl(var_name)), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

    if (nc_variables_deflate) then
        nc_status = nf90_def_var_deflate(nc_id, j, nc_shuffle_filter, 1, nc_deflate_level)
        if (nc_status /= NF90_NOERR) then
            if (.not. present(io_status)) then
                dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
                write(0,'("Error setting variable compression for ", a, " in group ", a)') trim(adjustl(var_name)), dummy_name(1:dummy_len)
                call nc_error(nf90_strerror(nc_status))
            end if
            io_status = 1
        end if
    end if

    output_varid%nc_id = nc_id
    output_varid%var_id = j
    output_varid%var_name = var_name
    output_varid%var_type = var_type

end function nc_create_var

subroutine nc_set_attrs_var(nc_var, attr_names, attr_values, io_status)

    type(t_ncVar), intent(inout)   :: nc_var
    character(len=*), intent(in)   :: attr_names(:), attr_values(:)
    integer                        :: i, dummy_int, dummy_len
    integer, optional, intent(out) :: io_status
    character(len=512)             :: fq_varname, fq_grpname

    allocate(nc_var%attr_names(size(attr_names)))
    allocate(nc_var%attr_values(size(attr_values)))

    nc_var%attr_names = attr_names(:)
    nc_var%attr_values = attr_values(:)

    nc_status = nf90_redef(nc_var%nc_id)
    if (present(io_status)) io_status = 0
    do i = 1, size(nc_var%attr_names)
        nc_status = nf90_put_att(nc_var%nc_id, nc_var%var_id, trim(adjustl(nc_var%attr_names(i))), trim(adjustl(nc_var%attr_values(i))))
        if (nc_status /= NF90_NOERR) then
            if (.not. present(io_status)) then
                dummy_int = nf90_inq_grpname_full(nc_var%nc_id, dummy_len, fq_grpname)
                dummy_int = nf90_inquire_variable(nc_var%nc_id, nc_var%var_id, fq_varname)
                write(0,'("Error setting attribute ", a, " for variable ", a, " in group ", a)') trim(adjustl(nc_var%attr_names(i))), trim(fq_varname), fq_grpname(1:dummy_len)
                call nc_error(nf90_strerror(nc_status))
            end if
            io_status = io_status + 1
            return
        end if
    end do

end subroutine nc_set_attrs_var

subroutine nc_set_attrs_glb(nc_id, attr_name, attr_value, io_status)

    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: attr_name, attr_value
    integer, optional, intent(out) :: io_status
    character(len=512)             :: dummy_name
    integer                        :: dummy_len, dummy_int

    if (present(io_status)) io_status = 0
    nc_status = nf90_redef(nc_id)
    nc_status = nf90_put_att(nc_id, NF90_GLOBAL, trim(adjustl(attr_name)), trim(adjustl(attr_value)))
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error setting attribute ", a, " for group ", a)') trim(adjustl(attr_name)), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status=1
    end if

end subroutine nc_set_attrs_glb

function nc_get_attr_glb(nc_id, attr_name, io_status) result (out_attr)
    ! Read a global attribute into the type nc_uni_att. The assignment
    ! operator is overloaded, so that
    ! attrib = nc_get_attr(nc_id, 'attribute name')
    ! automatically assigns attrib the correct value(s). Of course,
    ! nc_get_attr is also interfaced to provide reading global and variable
    ! attributes. The difference between the two is the extra argument
    ! 'var_name' found in the routine nc_get_attr_var.
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: attr_name
    integer, optional, intent(out) :: io_status
    type(nc_uni_att)               :: out_attr

    integer                        :: attr_type, attr_len, dummy_len, dummy_int
    character(len=512)             :: dummy_name

    if (present(io_status)) io_status = 0
    nc_status = nf90_inquire_attribute(nc_id, NF90_GLOBAL, attr_name, attr_type, attr_len)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error inquiring information about attribute ", a, " from group ", a)') trim(attr_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if

    select case (attr_type)
        case (NF90_BYTE)
            allocate(out_attr%i1_1d(attr_len))
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%i1_1d)
        case (NF90_CHAR)
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%c)
        case (NF90_SHORT)
            allocate(out_attr%i2_1d(attr_len))
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%i2_1d)
        case (NF90_INT)
            allocate(out_attr%i4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%i4_1d)
        case (NF90_FLOAT)
            allocate(out_attr%r4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%r4_1d)
        case (NF90_DOUBLE)
            allocate(out_attr%r8_1d(attr_len))
            nc_status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, out_attr%r8_1d)
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error reading attribute ", a, " from group ", a)') trim(attr_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
    end if

end function nc_get_attr_glb

function nc_get_attr_var(nc_id, var_name, attr_name, io_status) result (out_attr)
    ! Read a global attribute into the type nc_uni_att. The assignment
    ! operator is overloaded, so that
    ! attrib = nc_get_attr(nc_id, 'attribute name')
    ! automatically assigns attrib the correct value(s). Of course,
    ! nc_get_attr is also interfaced to provide reading global and variable
    ! attributes. The difference between the two is the extra argument
    ! 'var_name' found in the routine nc_get_attr_var.
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: attr_name, var_name
    integer, optional, intent(out) :: io_status
    type(nc_uni_att)               :: out_attr
    character(len=512)             :: dummy_name
    integer                        :: attr_type, attr_len, var_id, dummy_len, dummy_int

    if (present(io_status)) io_status = 0
    nc_status = nf90_inq_varid(nc_id, var_name, var_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not find variable ", a, " in group ", a)') var_name, dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if
    nc_status = nf90_inquire_attribute(nc_id, var_id, attr_name, attr_type, attr_len)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not find attribute ", a, " of variable ", a, " in group ", a)') attr_name, var_name, dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if

    select case (attr_type)
        case (NF90_BYTE)
            allocate(out_attr%i1_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i1_1d)
        case (NF90_CHAR)
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%c)
        case (NF90_SHORT)
            allocate(out_attr%i2_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i2_1d)
        case (NF90_INT)
            allocate(out_attr%i4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i4_1d)
        case (NF90_FLOAT)
            allocate(out_attr%r4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%r4_1d)
        case (NF90_DOUBLE)
            allocate(out_attr%r8_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%r8_1d)
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Could not read value of attribute ", a, " of variable ", a, " in group ", a)') attr_name, var_name, dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
    end if

end function nc_get_attr_var

function nc_get_attr_varid(nc_id, var_id, attr_name, io_status) result (out_attr)
    ! Read a global attribute into the type nc_uni_att. The assignment
    ! operator is overloaded, so that
    ! attrib = nc_get_attr(nc_id, 'attribute name')
    ! automatically assigns attrib the correct value(s). Of course,
    ! nc_get_attr is also interfaced to provide reading global and variable
    ! attributes. The difference between the two is the extra argument
    ! 'var_name' found in the routine nc_get_attr_var, or the argument
    ! 'var_id' found in nc_get_attr_varid.
    integer, intent(in)            :: nc_id, var_id
    character(len=*), intent(in)   :: attr_name
    integer, optional, intent(out) :: io_status
    type(nc_uni_att)               :: out_attr
    character(len=512)             :: fq_grpname, fq_varname
    integer                        :: attr_type, attr_len, dummy_int, dummy_len

    if (present(io_status)) io_status = 0
    nc_status = nf90_inquire_attribute(nc_id, var_id, attr_name, attr_type, attr_len)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, fq_grpname)
            dummy_int = nf90_inquire_variable(nc_id, var_id, fq_varname)
            write(0,'("Attribute ", a, " not found for variable ", a, " in group ", a)') trim(attr_name), trim(fq_varname), fq_grpname(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if

    select case (attr_type)
        case (NF90_BYTE)
            allocate(out_attr%i1_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i1_1d)
        case (NF90_CHAR)
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%c)
        case (NF90_SHORT)
            allocate(out_attr%i2_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i2_1d)
        case (NF90_INT)
            allocate(out_attr%i4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%i4_1d)
        case (NF90_FLOAT)
            allocate(out_attr%r4_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%r4_1d)
        case (NF90_DOUBLE)
            allocate(out_attr%r8_1d(attr_len))
            nc_status = nf90_get_att(nc_id, var_id, attr_name, out_attr%r8_1d)
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, fq_grpname)
            dummy_int = nf90_inquire_variable(nc_id, var_id, fq_varname)
            write(0,'("Error reading attribute ", a, " of variable ", a, " from group ", a)') trim(attr_name), trim(fq_varname), fq_grpname(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
    end if

end function nc_get_attr_varid

subroutine nc_attr_to_scalar_i1(out_var, in_nc_att)

    type(nc_uni_att), intent(in) :: in_nc_att
    integer(1), intent(out)      :: out_var

    if (allocated(in_nc_att%i1_1d) .and. size(in_nc_att%i1_1d) == 1) out_var = in_nc_att%i1_1d(1)

end subroutine nc_attr_to_scalar_i1

subroutine nc_attr_to_array_i1(out_var, in_nc_att)

    type(nc_uni_att), intent(in)         :: in_nc_att
    integer(1), allocatable, intent(out) :: out_var(:)
    integer                              :: arr_len

    if (allocated(in_nc_att%i1_1d)) then
        arr_len = size(in_nc_att%i1_1d)
        if (allocated(out_var)) deallocate(out_var)
        allocate(out_var(arr_len))
        out_var(:) = in_nc_att%i1_1d(:)
    end if

end subroutine nc_attr_to_array_i1

subroutine nc_attr_to_scalar_c(out_var, in_nc_att)

    type(nc_uni_att), intent(in)  :: in_nc_att
    character(len=*), intent(out) :: out_var

    out_var = trim(adjustl(in_nc_att%c))

end subroutine nc_attr_to_scalar_c

subroutine nc_attr_to_scalar_i2(out_var, in_nc_att)

    type(nc_uni_att), intent(in) :: in_nc_att
    integer(2), intent(out)      :: out_var

    if (allocated(in_nc_att%i2_1d) .and. size(in_nc_att%i2_1d) == 1) out_var = in_nc_att%i2_1d(1)

end subroutine nc_attr_to_scalar_i2

subroutine nc_attr_to_array_i2(out_var, in_nc_att)

    type(nc_uni_att), intent(in)         :: in_nc_att
    integer(2), allocatable, intent(out) :: out_var(:)
    integer                              :: arr_len

    if (allocated(in_nc_att%i2_1d)) then
        arr_len = size(in_nc_att%i2_1d)
        if (allocated(out_var)) deallocate(out_var)
        allocate(out_var(arr_len))
        out_var(:) = in_nc_att%i2_1d(:)
    end if

end subroutine nc_attr_to_array_i2

subroutine nc_attr_to_scalar_i4(out_var, in_nc_att)

    type(nc_uni_att), intent(in) :: in_nc_att
    integer(4), intent(out)      :: out_var

    if (allocated(in_nc_att%i4_1d) .and. size(in_nc_att%i4_1d) == 1) out_var = in_nc_att%i4_1d(1)

end subroutine nc_attr_to_scalar_i4

subroutine nc_attr_to_array_i4(out_var, in_nc_att)

    type(nc_uni_att), intent(in)         :: in_nc_att
    integer(4), allocatable, intent(out) :: out_var(:)
    integer                              :: arr_len

    if (allocated(in_nc_att%i4_1d)) then
        arr_len = size(in_nc_att%i4_1d)
        if (allocated(out_var)) deallocate(out_var)
        allocate(out_var(arr_len))
        out_var(:) = in_nc_att%i4_1d(:)
    end if

end subroutine nc_attr_to_array_i4

subroutine nc_attr_to_scalar_r4(out_var, in_nc_att)

    type(nc_uni_att), intent(in) :: in_nc_att
    real(4), intent(out)         :: out_var

    if (allocated(in_nc_att%r4_1d) .and. size(in_nc_att%r4_1d) == 1) out_var = in_nc_att%r4_1d(1)

end subroutine nc_attr_to_scalar_r4

subroutine nc_attr_to_array_r4(out_var, in_nc_att)

    type(nc_uni_att), intent(in)      :: in_nc_att
    real(4), allocatable, intent(out) :: out_var(:)
    integer                           :: arr_len

    if (allocated(in_nc_att%r4_1d)) then
        arr_len = size(in_nc_att%r4_1d)
        if (allocated(out_var)) deallocate(out_var)
        allocate(out_var(arr_len))
        out_var(:) = in_nc_att%r4_1d(:)
    end if

end subroutine nc_attr_to_array_r4

subroutine nc_attr_to_scalar_r8(out_var, in_nc_att)

    type(nc_uni_att), intent(in) :: in_nc_att
    real(8), intent(out)         :: out_var

    if (allocated(in_nc_att%r8_1d) .and. size(in_nc_att%r8_1d) == 1) out_var = in_nc_att%r8_1d(1)

end subroutine nc_attr_to_scalar_r8

subroutine nc_attr_to_array_r8(out_var, in_nc_att)

    type(nc_uni_att), intent(in)      :: in_nc_att
    real(8), allocatable, intent(out) :: out_var(:)
    integer                           :: arr_len

    if (allocated(in_nc_att%r8_1d)) then
        arr_len = size(in_nc_att%r8_1d)
        if (allocated(out_var)) deallocate(out_var)
        allocate(out_var(arr_len))
        out_var(:) = in_nc_att%r8_1d(:)
    end if

end subroutine nc_attr_to_array_r8

function nc_read_var(nc_id, var_name, io_status) result(out_var)
    ! Read a variable directly into an array. The call nc_read_var() returns a
    ! variable of type nc_uni_var. However, the assignment (=) is overloaded so
    ! that x = y, where x is an allocatable array and y is a nc_uni_var,
    ! automatically copies the relevant element of nc_uni_var into x, depending
    ! on the shape and type of x. In summary, to read a variable 'var_name' into
    ! an array var of the correct shape and type, execute:
    !   var = nc_read_var(nc_id, var_name)
    ! where nc_id can also be a group ID. Note that nc_uni_var, being an
    ! intent(in) variable for the assignment(=), is not deallocated. I don't
    ! know if this can result in memory leaks. Needs testing.
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: var_name
    integer, optional, intent(out) :: io_status
    type(nc_uni_var)               :: out_var
    integer                        :: var_id, n_atts, var_type, n_dims, dimids(NF90_MAX_VAR_DIMS), var_dims(NF90_MAX_VAR_DIMS), i
    character(len=NF90_MAX_NAME)   :: useless_name
    character(len=512)             :: dummy_name
    integer                        :: dummy_int, dummy_len

    if (present(io_status)) io_status = 0
    nc_status = nf90_inq_varid(nc_id, var_name, var_id)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Variable ", a, " not found in group ", a)') trim(var_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if
    nc_status = nf90_inquire_variable(nc_id, var_id, useless_name, var_type, n_dims, dimids, n_atts)
    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error reading variable information for variable ", a, " from group ", a)') trim(var_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
        return
    end if

    do i=1,n_dims
        nc_status = nf90_inquire_dimension(nc_id, dimids(i), useless_name, var_dims(i))
    end do

    select case (var_type)
        case (NF90_DOUBLE)
            select case (n_dims)
                case (1)
                    allocate(out_var%r8_1d(var_dims(1)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_1d)
                case (2)
                    allocate(out_var%r8_2d(var_dims(1),var_dims(2)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_2d)
                case (3)
                    allocate(out_var%r8_3d(var_dims(1),var_dims(2),var_dims(3)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_3d)
                case (4)
                    allocate(out_var%r8_4d(var_dims(1),var_dims(2),var_dims(3),var_dims(4)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_4d)
                case (5)
                    allocate(out_var%r8_5d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_5d)
                case (6)
                    allocate(out_var%r8_6d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_6d)
                case (7)
                    allocate(out_var%r8_7d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6),var_dims(7)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r8_7d)
            end select
        case (NF90_FLOAT)
            select case (n_dims)
                case (1)
                    allocate(out_var%r4_1d(var_dims(1)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_1d)
                case (2)
                    allocate(out_var%r4_2d(var_dims(1),var_dims(2)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_2d)
                case (3)
                    allocate(out_var%r4_3d(var_dims(1),var_dims(2),var_dims(3)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_3d)
                case (4)
                    allocate(out_var%r4_4d(var_dims(1),var_dims(2),var_dims(3),var_dims(4)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_4d)
                case (5)
                    allocate(out_var%r4_5d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_5d)
                case (6)
                    allocate(out_var%r4_6d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_6d)
                case (7)
                    allocate(out_var%r4_7d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6),var_dims(7)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%r4_7d)
            end select
        case (NF90_INT)
            select case (n_dims)
                case (1)
                    allocate(out_var%i4_1d(var_dims(1)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_1d)
                case (2)
                    allocate(out_var%i4_2d(var_dims(1),var_dims(2)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_2d)
                case (3)
                    allocate(out_var%i4_3d(var_dims(1),var_dims(2),var_dims(3)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_3d)
                case (4)
                    allocate(out_var%i4_4d(var_dims(1),var_dims(2),var_dims(3),var_dims(4)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_4d)
                case (5)
                    allocate(out_var%i4_5d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_5d)
                case (6)
                    allocate(out_var%i4_6d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_6d)
                case (7)
                    allocate(out_var%i4_7d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6),var_dims(7)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i4_7d)
            end select
        case (NF90_SHORT)
            select case (n_dims)
                case (1)
                    allocate(out_var%i2_1d(var_dims(1)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_1d)
                case (2)
                    allocate(out_var%i2_2d(var_dims(1),var_dims(2)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_2d)
                case (3)
                    allocate(out_var%i2_3d(var_dims(1),var_dims(2),var_dims(3)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_3d)
                case (4)
                    allocate(out_var%i2_4d(var_dims(1),var_dims(2),var_dims(3),var_dims(4)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_4d)
                case (5)
                    allocate(out_var%i2_5d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_5d)
                case (6)
                    allocate(out_var%i2_6d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_6d)
                case (7)
                    allocate(out_var%i2_7d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6),var_dims(7)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i2_7d)
            end select
        case (NF90_BYTE)
            select case (n_dims)
                case (1)
                    allocate(out_var%i1_1d(var_dims(1)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_1d)
                case (2)
                    allocate(out_var%i1_2d(var_dims(1),var_dims(2)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_2d)
                case (3)
                    allocate(out_var%i1_3d(var_dims(1),var_dims(2),var_dims(3)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_3d)
                case (4)
                    allocate(out_var%i1_4d(var_dims(1),var_dims(2),var_dims(3),var_dims(4)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_4d)
                case (5)
                    allocate(out_var%i1_5d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_5d)
                case (6)
                    allocate(out_var%i1_6d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_6d)
                case (7)
                    allocate(out_var%i1_7d(var_dims(1),var_dims(2),var_dims(3),var_dims(4),var_dims(5),var_dims(6),var_dims(7)))
                    nc_status = nf90_get_var(nc_id, var_id, out_var%i1_7d)
            end select
    end select

    if (nc_status /= NF90_NOERR) then
        if (.not. present(io_status)) then
            dummy_int = nf90_inq_grpname_full(nc_id, dummy_len, dummy_name)
            write(0,'("Error reading value of variable ", a, " from group ", a)') trim(var_name), dummy_name(1:dummy_len)
            call nc_error(nf90_strerror(nc_status))
        end if
        io_status = 1
    end if

end function nc_read_var

subroutine nc_var_to_array_i1_1d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:)
    integer                              :: arr_shape(1)

    if (allocated(in_nc_var%i1_1d)) then
        arr_shape = shape(in_nc_var%i1_1d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1)))
        out_array(:) = in_nc_var%i1_1d(:)
    end if

end subroutine nc_var_to_array_i1_1d

subroutine nc_var_to_array_i1_2d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:)
    integer                              :: arr_shape(2)

    if (allocated(in_nc_var%i1_2d)) then
        arr_shape = shape(in_nc_var%i1_2d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2)))
        out_array(:,:) = in_nc_var%i1_2d(:,:)
    end if

end subroutine nc_var_to_array_i1_2d

subroutine nc_var_to_array_i1_3d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:,:)
    integer                              :: arr_shape(3)

    if (allocated(in_nc_var%i1_3d)) then
        arr_shape = shape(in_nc_var%i1_3d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3)))
        out_array(:,:,:) = in_nc_var%i1_3d(:,:,:)
    end if

end subroutine nc_var_to_array_i1_3d

subroutine nc_var_to_array_i1_4d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:,:,:)
    integer                              :: arr_shape(4)

    if (allocated(in_nc_var%i1_4d)) then
        arr_shape = shape(in_nc_var%i1_4d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4)))
        out_array(:,:,:,:) = in_nc_var%i1_4d(:,:,:,:)
    end if

end subroutine nc_var_to_array_i1_4d

subroutine nc_var_to_array_i1_5d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:,:,:,:)
    integer                              :: arr_shape(5)

    if (allocated(in_nc_var%i1_5d)) then
        arr_shape = shape(in_nc_var%i1_5d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5)))
        out_array(:,:,:,:,:) = in_nc_var%i1_5d(:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i1_5d

subroutine nc_var_to_array_i1_6d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:,:,:,:,:)
    integer                              :: arr_shape(6)

    if (allocated(in_nc_var%i1_6d)) then
        arr_shape = shape(in_nc_var%i1_6d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6)))
        out_array(:,:,:,:,:,:) = in_nc_var%i1_6d(:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i1_6d

subroutine nc_var_to_array_i1_7d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(1), allocatable, intent(out) :: out_array(:,:,:,:,:,:,:)
    integer                              :: arr_shape(7)

    if (allocated(in_nc_var%i1_7d)) then
        arr_shape = shape(in_nc_var%i1_7d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6),arr_shape(7)))
        out_array(:,:,:,:,:,:,:) = in_nc_var%i1_7d(:,:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i1_7d

subroutine nc_var_to_array_i2_1d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:)
    integer                              :: arr_shape(1)

    if (allocated(in_nc_var%i2_1d)) then
        arr_shape = shape(in_nc_var%i2_1d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1)))
        out_array(:) = in_nc_var%i2_1d(:)
    end if

end subroutine nc_var_to_array_i2_1d

subroutine nc_var_to_array_i2_2d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:)
    integer                              :: arr_shape(2)

    if (allocated(in_nc_var%i2_2d)) then
        arr_shape = shape(in_nc_var%i2_2d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2)))
        out_array(:,:) = in_nc_var%i2_2d(:,:)
    end if

end subroutine nc_var_to_array_i2_2d

subroutine nc_var_to_array_i2_3d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:,:)
    integer                              :: arr_shape(3)

    if (allocated(in_nc_var%i2_3d)) then
        arr_shape = shape(in_nc_var%i2_3d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3)))
        out_array(:,:,:) = in_nc_var%i2_3d(:,:,:)
    end if

end subroutine nc_var_to_array_i2_3d

subroutine nc_var_to_array_i2_4d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:,:,:)
    integer                              :: arr_shape(4)

    if (allocated(in_nc_var%i2_4d)) then
        arr_shape = shape(in_nc_var%i2_4d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4)))
        out_array(:,:,:,:) = in_nc_var%i2_4d(:,:,:,:)
    end if

end subroutine nc_var_to_array_i2_4d

subroutine nc_var_to_array_i2_5d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:,:,:,:)
    integer                              :: arr_shape(5)

    if (allocated(in_nc_var%i2_5d)) then
        arr_shape = shape(in_nc_var%i2_5d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5)))
        out_array(:,:,:,:,:) = in_nc_var%i2_5d(:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i2_5d

subroutine nc_var_to_array_i2_6d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:,:,:,:,:)
    integer                              :: arr_shape(6)

    if (allocated(in_nc_var%i2_6d)) then
        arr_shape = shape(in_nc_var%i2_6d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6)))
        out_array(:,:,:,:,:,:) = in_nc_var%i2_6d(:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i2_6d

subroutine nc_var_to_array_i2_7d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(2), allocatable, intent(out) :: out_array(:,:,:,:,:,:,:)
    integer                              :: arr_shape(7)

    if (allocated(in_nc_var%i2_7d)) then
        arr_shape = shape(in_nc_var%i2_7d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6),arr_shape(7)))
        out_array(:,:,:,:,:,:,:) = in_nc_var%i2_7d(:,:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i2_7d

subroutine nc_var_to_array_i4_1d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:)
    integer                              :: arr_shape(1)

    if (allocated(in_nc_var%i4_1d)) then
        arr_shape = shape(in_nc_var%i4_1d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1)))
        out_array(:) = in_nc_var%i4_1d(:)
    end if

end subroutine nc_var_to_array_i4_1d

subroutine nc_var_to_array_i4_2d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:)
    integer                              :: arr_shape(2)

    if (allocated(in_nc_var%i4_2d)) then
        arr_shape = shape(in_nc_var%i4_2d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2)))
        out_array(:,:) = in_nc_var%i4_2d(:,:)
    end if

end subroutine nc_var_to_array_i4_2d

subroutine nc_var_to_array_i4_3d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:,:)
    integer                              :: arr_shape(3)

    if (allocated(in_nc_var%i4_3d)) then
        arr_shape = shape(in_nc_var%i4_3d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3)))
        out_array(:,:,:) = in_nc_var%i4_3d(:,:,:)
    end if

end subroutine nc_var_to_array_i4_3d

subroutine nc_var_to_array_i4_4d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:,:,:)
    integer                              :: arr_shape(4)

    if (allocated(in_nc_var%i4_4d)) then
        arr_shape = shape(in_nc_var%i4_4d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4)))
        out_array(:,:,:,:) = in_nc_var%i4_4d(:,:,:,:)
    end if

end subroutine nc_var_to_array_i4_4d

subroutine nc_var_to_array_i4_5d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:,:,:,:)
    integer                              :: arr_shape(5)

    if (allocated(in_nc_var%i4_5d)) then
        arr_shape = shape(in_nc_var%i4_5d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5)))
        out_array(:,:,:,:,:) = in_nc_var%i4_5d(:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i4_5d

subroutine nc_var_to_array_i4_6d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:,:,:,:,:)
    integer                              :: arr_shape(6)

    if (allocated(in_nc_var%i4_6d)) then
        arr_shape = shape(in_nc_var%i4_6d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6)))
        out_array(:,:,:,:,:,:) = in_nc_var%i4_6d(:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i4_6d

subroutine nc_var_to_array_i4_7d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    integer(4), allocatable, intent(out) :: out_array(:,:,:,:,:,:,:)
    integer                              :: arr_shape(7)

    if (allocated(in_nc_var%i4_7d)) then
        arr_shape = shape(in_nc_var%i4_7d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6),arr_shape(7)))
        out_array(:,:,:,:,:,:,:) = in_nc_var%i4_7d(:,:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_i4_7d

subroutine nc_var_to_array_r4_1d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:)
    integer                              :: arr_shape(1)

    if (allocated(in_nc_var%r4_1d)) then
        arr_shape = shape(in_nc_var%r4_1d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1)))
        out_array(:) = in_nc_var%r4_1d(:)
    end if

end subroutine nc_var_to_array_r4_1d

subroutine nc_var_to_array_r4_2d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:)
    integer                              :: arr_shape(2)

    if (allocated(in_nc_var%r4_2d)) then
        arr_shape = shape(in_nc_var%r4_2d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2)))
        out_array(:,:) = in_nc_var%r4_2d(:,:)
    end if

end subroutine nc_var_to_array_r4_2d

subroutine nc_var_to_array_r4_3d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:,:)
    integer                              :: arr_shape(3)

    if (allocated(in_nc_var%r4_3d)) then
        arr_shape = shape(in_nc_var%r4_3d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3)))
        out_array(:,:,:) = in_nc_var%r4_3d(:,:,:)
    end if

end subroutine nc_var_to_array_r4_3d

subroutine nc_var_to_array_r4_4d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:,:,:)
    integer                              :: arr_shape(4)

    if (allocated(in_nc_var%r4_4d)) then
        arr_shape = shape(in_nc_var%r4_4d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4)))
        out_array(:,:,:,:) = in_nc_var%r4_4d(:,:,:,:)
    end if

end subroutine nc_var_to_array_r4_4d

subroutine nc_var_to_array_r4_5d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:,:,:,:)
    integer                              :: arr_shape(5)

    if (allocated(in_nc_var%r4_5d)) then
        arr_shape = shape(in_nc_var%r4_5d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5)))
        out_array(:,:,:,:,:) = in_nc_var%r4_5d(:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r4_5d

subroutine nc_var_to_array_r4_6d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:,:,:,:,:)
    integer                              :: arr_shape(6)

    if (allocated(in_nc_var%r4_6d)) then
        arr_shape = shape(in_nc_var%r4_6d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6)))
        out_array(:,:,:,:,:,:) = in_nc_var%r4_6d(:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r4_6d

subroutine nc_var_to_array_r4_7d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(4), allocatable, intent(out) :: out_array(:,:,:,:,:,:,:)
    integer                              :: arr_shape(7)

    if (allocated(in_nc_var%r4_7d)) then
        arr_shape = shape(in_nc_var%r4_7d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6),arr_shape(7)))
        out_array(:,:,:,:,:,:,:) = in_nc_var%r4_7d(:,:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r4_7d

subroutine nc_var_to_array_r8_1d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:)
    integer                              :: arr_shape(1)

    if (allocated(in_nc_var%r8_1d)) then
        arr_shape = shape(in_nc_var%r8_1d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1)))
        out_array(:) = in_nc_var%r8_1d(:)
    end if

end subroutine nc_var_to_array_r8_1d

subroutine nc_var_to_array_r8_2d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:)
    integer                              :: arr_shape(2)

    if (allocated(in_nc_var%r8_2d)) then
        arr_shape = shape(in_nc_var%r8_2d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2)))
        out_array(:,:) = in_nc_var%r8_2d(:,:)
    end if

end subroutine nc_var_to_array_r8_2d

subroutine nc_var_to_array_r8_3d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:,:)
    integer                              :: arr_shape(3)

    if (allocated(in_nc_var%r8_3d)) then
        arr_shape = shape(in_nc_var%r8_3d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3)))
        out_array(:,:,:) = in_nc_var%r8_3d(:,:,:)
    end if

end subroutine nc_var_to_array_r8_3d

subroutine nc_var_to_array_r8_4d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:,:,:)
    integer                              :: arr_shape(4)

    if (allocated(in_nc_var%r8_4d)) then
        arr_shape = shape(in_nc_var%r8_4d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4)))
        out_array(:,:,:,:) = in_nc_var%r8_4d(:,:,:,:)
    end if

end subroutine nc_var_to_array_r8_4d

subroutine nc_var_to_array_r8_5d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:,:,:,:)
    integer                              :: arr_shape(5)

    if (allocated(in_nc_var%r8_5d)) then
        arr_shape = shape(in_nc_var%r8_5d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5)))
        out_array(:,:,:,:,:) = in_nc_var%r8_5d(:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r8_5d

subroutine nc_var_to_array_r8_6d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:,:,:,:,:)
    integer                              :: arr_shape(6)

    if (allocated(in_nc_var%r8_6d)) then
        arr_shape = shape(in_nc_var%r8_6d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6)))
        out_array(:,:,:,:,:,:) = in_nc_var%r8_6d(:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r8_6d

subroutine nc_var_to_array_r8_7d(out_array, in_nc_var)

    type(nc_uni_var), intent(in)         :: in_nc_var
    real(8), allocatable, intent(out) :: out_array(:,:,:,:,:,:,:)
    integer                              :: arr_shape(7)

    if (allocated(in_nc_var%r8_7d)) then
        arr_shape = shape(in_nc_var%r8_7d)
        if (allocated(out_array)) deallocate(out_array)
        allocate(out_array(arr_shape(1),arr_shape(2),arr_shape(3),arr_shape(4),arr_shape(5),arr_shape(6),arr_shape(7)))
        out_array(:,:,:,:,:,:,:) = in_nc_var%r8_7d(:,:,:,:,:,:,:)
    end if

end subroutine nc_var_to_array_r8_7d

subroutine nc_dump_var_i1_1d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:))

end subroutine nc_dump_var_i1_1d

subroutine nc_dump_var_i1_2d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:))

end subroutine nc_dump_var_i1_2d

subroutine nc_dump_var_i1_3d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:))

end subroutine nc_dump_var_i1_3d

subroutine nc_dump_var_i1_4d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:))

end subroutine nc_dump_var_i1_4d

subroutine nc_dump_var_i1_5d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:))

end subroutine nc_dump_var_i1_5d

subroutine nc_dump_var_i1_6d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:))

end subroutine nc_dump_var_i1_6d

subroutine nc_dump_var_i1_7d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(1), intent(in)         :: var_value(:,:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'byte')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:,:))

end subroutine nc_dump_var_i1_7d

subroutine nc_dump_var_i2_1d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:))

end subroutine nc_dump_var_i2_1d

subroutine nc_dump_var_i2_2d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:))

end subroutine nc_dump_var_i2_2d

subroutine nc_dump_var_i2_3d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:))

end subroutine nc_dump_var_i2_3d

subroutine nc_dump_var_i2_4d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:))

end subroutine nc_dump_var_i2_4d

subroutine nc_dump_var_i2_5d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:))

end subroutine nc_dump_var_i2_5d

subroutine nc_dump_var_i2_6d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:))

end subroutine nc_dump_var_i2_6d

subroutine nc_dump_var_i2_7d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(2), intent(in)         :: var_value(:,:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'short')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:,:))

end subroutine nc_dump_var_i2_7d

subroutine nc_dump_var_i4_1d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:))

end subroutine nc_dump_var_i4_1d

subroutine nc_dump_var_i4_2d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:))

end subroutine nc_dump_var_i4_2d

subroutine nc_dump_var_i4_3d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:))

end subroutine nc_dump_var_i4_3d

subroutine nc_dump_var_i4_4d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:))

end subroutine nc_dump_var_i4_4d

subroutine nc_dump_var_i4_5d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:))

end subroutine nc_dump_var_i4_5d

subroutine nc_dump_var_i4_6d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:))

end subroutine nc_dump_var_i4_6d

subroutine nc_dump_var_i4_7d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    integer(4), intent(in)         :: var_value(:,:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'int')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:,:))

end subroutine nc_dump_var_i4_7d

subroutine nc_dump_var_r4_1d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:))

end subroutine nc_dump_var_r4_1d

subroutine nc_dump_var_r4_2d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:))

end subroutine nc_dump_var_r4_2d

subroutine nc_dump_var_r4_3d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:))

end subroutine nc_dump_var_r4_3d

subroutine nc_dump_var_r4_4d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:))

end subroutine nc_dump_var_r4_4d

subroutine nc_dump_var_r4_5d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:))

end subroutine nc_dump_var_r4_5d

subroutine nc_dump_var_r4_6d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:))

end subroutine nc_dump_var_r4_6d

subroutine nc_dump_var_r4_7d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(4), intent(in)            :: var_value(:,:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'float')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:,:))

end subroutine nc_dump_var_r4_7d

subroutine nc_dump_var_r8_1d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:))

end subroutine nc_dump_var_r8_1d

subroutine nc_dump_var_r8_2d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:))

end subroutine nc_dump_var_r8_2d

subroutine nc_dump_var_r8_3d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:))

end subroutine nc_dump_var_r8_3d

subroutine nc_dump_var_r8_4d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:))

end subroutine nc_dump_var_r8_4d

subroutine nc_dump_var_r8_5d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:))

end subroutine nc_dump_var_r8_5d

subroutine nc_dump_var_r8_6d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:))

end subroutine nc_dump_var_r8_6d

subroutine nc_dump_var_r8_7d(nc_id, var_name, dim_list, var_value, attr_names, attr_values)

    real(8), intent(in)            :: var_value(:,:,:,:,:,:,:)
    character(len=*), intent(in)   :: var_name
    integer, intent(in)            :: nc_id
    character(len=*), intent(in)   :: dim_list(:)
    type(t_ncVar)                  :: nc_var
    character(len=*), optional, intent(in)   :: attr_names(:), attr_values(:)

    nc_status = nf90_redef(nc_id)
    nc_var = nc_create_var(nc_id, var_name, dim_list, 'double')
    if (present(attr_names) .and. present(attr_values)) call nc_set_attrs_var(nc_var, attr_names, attr_values)
    call nc_set_var(nc_var, var_value(:,:,:,:,:,:,:))

end subroutine nc_dump_var_r8_7d

subroutine nc_set_var_i1_1d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_1d

subroutine nc_set_var_i1_2d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_2d

subroutine nc_set_var_i1_3d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_3d

subroutine nc_set_var_i1_4d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_4d

subroutine nc_set_var_i1_5d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_5d

subroutine nc_set_var_i1_6d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_6d

subroutine nc_set_var_i1_7d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(1), intent(in)     :: var_value(:,:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i1_7d

subroutine nc_set_var_i2_1d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_1d

subroutine nc_set_var_i2_2d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_2d

subroutine nc_set_var_i2_3d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_3d

subroutine nc_set_var_i2_4d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_4d

subroutine nc_set_var_i2_5d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_5d

subroutine nc_set_var_i2_6d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_6d

subroutine nc_set_var_i2_7d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(2), intent(in)     :: var_value(:,:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i2_7d

subroutine nc_set_var_i4_1d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_1d

subroutine nc_set_var_i4_2d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_2d

subroutine nc_set_var_i4_3d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_3d

subroutine nc_set_var_i4_4d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_4d

subroutine nc_set_var_i4_5d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_5d

subroutine nc_set_var_i4_6d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_6d

subroutine nc_set_var_i4_7d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    integer(4), intent(in)     :: var_value(:,:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_i4_7d

subroutine nc_set_var_r4_1d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_1d

subroutine nc_set_var_r4_2d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_2d

subroutine nc_set_var_r4_3d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_3d

subroutine nc_set_var_r4_4d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_4d

subroutine nc_set_var_r4_5d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_5d

subroutine nc_set_var_r4_6d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_6d

subroutine nc_set_var_r4_7d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(4), intent(in)     :: var_value(:,:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r4_7d

subroutine nc_set_var_r8_1d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_1d

subroutine nc_set_var_r8_2d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_2d

subroutine nc_set_var_r8_3d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_3d

subroutine nc_set_var_r8_4d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_4d

subroutine nc_set_var_r8_5d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_5d

subroutine nc_set_var_r8_6d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_6d

subroutine nc_set_var_r8_7d(nc_var, var_value)
    type(t_ncVar), intent(in)  :: nc_var
    real(8), intent(in)     :: var_value(:,:,:,:,:,:,:)
    nc_status = nf90_enddef(nc_var%nc_id)
    nc_status = nf90_put_var(nc_var%nc_id, nc_var%var_id, var_value(:,:,:,:,:,:,:))
    nc_status = nf90_sync(nc_var%nc_id)
end subroutine nc_set_var_r8_7d

!function pad(input_str, length) result (output_str)
!
!    character(len=*), intent(in)  :: input_str
!    character(len=length)         :: output_str
!    integer, intent(in)           :: length
!
!    output_str = input_str
!
!end function pad

end module file_netcdf
