! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module point_mod

  implicit none

  integer, allocatable, dimension (:) :: ireleasestart
  integer, allocatable, dimension (:) :: ireleaseend
  integer, allocatable, dimension (:) :: npart
  integer*2, allocatable, dimension (:) :: kindz

  real,allocatable, dimension (:) :: xpoint1
  real,allocatable, dimension (:) :: xpoint2
  real,allocatable, dimension (:) :: ypoint1
  real,allocatable, dimension (:) :: ypoint2
  real,allocatable, dimension (:) :: zpoint1
  real,allocatable, dimension (:) :: zpoint2

  real,allocatable, dimension (:,:) :: xmass
  real,allocatable, dimension (:) :: rho_rel

end module point_mod
