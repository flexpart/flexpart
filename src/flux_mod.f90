! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module flux_mod

  ! flux eastward, westward, northward, southward, upward and downward
  ! fluxes of all species and all ageclasses
  ! areaeast,areanorth [m2] side areas of each grid cell

  implicit none

  real,allocatable, dimension (:,:,:,:,:,:,:) :: flux

  !1 fluxw west - east
  !2 fluxe east - west
  !3 fluxs south - north
  !4 fluxn north - south
  !5 fluxu upward
  !6 fluxd downward
  !real,allocatable, dimension (:,:,:) :: areanorth
  !real,allocatable, dimension (:,:,:) :: areaeast

end module flux_mod
