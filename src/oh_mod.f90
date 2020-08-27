! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module oh_mod

  !includes OH concentration field as well as the height information
  !for this field

  implicit none

  integer :: nxOH,nyOH,nzOH
  real, allocatable, dimension(:) :: lonOH,latOH,altOH
  real, allocatable, dimension(:,:,:,:) :: OH_hourly
  real, allocatable, dimension (:,:,:,:) :: OH_field
  real, dimension(2) :: memOHtime
  real, dimension(360,180,12) :: jrate_average
  real, dimension(360) :: lonjr
  real, dimension(180) :: latjr

end module oh_mod
