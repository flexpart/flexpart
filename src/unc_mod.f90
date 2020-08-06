! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

! DJM - 2017-05-09 - added #ifdef USE_MPIINPLACE cpp directive to     *
! enable declaration of a gridunc0 array if required by MPI code in   *
! mpi_mod.f90                                                         *
!                                                                     *
!**********************************************************************

module unc_mod

  use par_mod, only:dep_prec

  implicit none

  real,allocatable, dimension (:,:,:,:,:,:,:) :: gridunc
#ifdef USE_MPIINPLACE
#else
  ! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
  ! then an aux array is needed for parallel grid reduction
  real,allocatable, dimension (:,:,:,:,:,:,:) :: gridunc0
  real,allocatable, dimension (:,:,:,:,:,:,:) :: griduncn0
#endif
  real,allocatable, dimension (:,:,:,:,:,:,:) :: griduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn

! For sum of individual contributions, used for the MPI version
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn0

  real,allocatable, dimension (:,:,:,:,:) :: init_cond

end module unc_mod
