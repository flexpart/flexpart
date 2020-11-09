! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine readlanduse

  !*****************************************************************************
  !                                                                            *
  !      Reads the landuse inventory into memory and relates it to Leaf Area   *
  !      Index and roughness length.                                           *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 10 January 1994                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! i                       loop indices                                       *
  ! landinvent(1200,600,13) area fractions of 13 landuse categories            *
  ! LENGTH(numpath)         length of the path names                           *
  ! PATH(numpath)           contains the path names                            *
  ! unitland                unit connected with landuse inventory              *
  !                                                                            *
  ! -----                                                                      *
  ! Sabine Eckhardt, Dec 06 - new landuse inventary                            *
  ! after                                                                      *
  ! Belward, A.S., Estes, J.E., and Kline, K.D., 1999,                         *
  ! The IGBP-DIS 1-Km Land-Cover Data Set DISCover:                            *
  ! A Project Overview: Photogrammetric Engineering and Remote Sensing,        *
  ! v. 65, no. 9, p. 1013-1020                                                 *
  !                                                                            *
  ! LANDUSE CATEGORIES:                                                        *
  !                                                                            *
  ! 1   Urban land                                                             *
  ! 2   Agricultural land                                  *
  ! 3   Range land                                         *
  ! 4   Deciduous forest                                           *
  ! 5   Coniferous forest                                              *
  ! 6   Mixed forest including wetland                                         *
  ! 7   water, both salt and fresh                                             *
  ! 8   barren land mostly desert                                              *
  ! 9   nonforested wetland                                                *
  ! 10  mixed agricultural and range land                              *
  ! 11  rocky open areas with low growing shrubs                               *
  ! 12  ice                                                                    *
  ! 13  rainforest                                                             *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: ix,jy,i,k,lu_cat,lu_perc
  integer(kind=1) :: ilr
  integer(kind=1) :: ilr_buffer(2160000)
  integer :: il,irecread
  real :: rlr, r2lr


  ! Read landuse inventory
  !***********************
  ! The landuse information is saved in a compressed format and written
  ! out by records of the length of 1 BYTE. Each grid cell consists of 3
  ! Bytes, which include 3 landuse categories (val 1-13 and 16 percentage
  ! categories) So one half byte is used to store the Landusecat the other
  ! for the percentageclass in 6.25 steps (100/6.25=16)
  ! e.g.
  ! 4 3  percentage 4 = 4*6.25 => 25% landuse class 3
  ! 2 1  percentage 2 = 2*6.25 => 13% landuse class 1
  ! 1 12 percentage 1 = 1*6.26 => 6.25% landuse class 12

  open(unitland,file=path(1)(1:length(1)) &
       //'IGBP_int1.dat',status='old', &
  !    +form='UNFORMATTED', err=998)
       form='UNFORMATTED', err=998, convert='little_endian')
!  print*,unitland
  read (unitland) (ilr_buffer(i),i=1,2160000)
  close(unitland)

  irecread=1
  do ix=1,1200
    do jy=1,600
  ! the 3 most abundant landuse categories in the inventory
  ! first half byte contains the landuse class
  ! second half byte contains the respective percentage
      do k=1,3
  ! 1 byte is read
        ilr=ilr_buffer(irecread)
  !      ilr=0
        irecread=irecread+1
  ! as only signed integer values exist an unsigned value is constructed
        if (ilr.lt.0) then
           il=ilr+256
        else
           il=ilr
        endif
  ! dividing by 16 has the effect to get rid of the right half of the byte
  ! so just the left half remains, this corresponds to a shift right of 4
  ! bits
        rlr=real(il)/16.
        lu_cat=int(rlr)
  ! the left half of the byte is substracted from the whole in order to
  ! get only the right half of the byte
        r2lr=rlr-int(rlr)
  ! shift left by 4
        lu_perc=r2lr*16.
        landinvent(ix,jy,k)=lu_cat
        landinvent(ix,jy,k+3)=lu_perc
  !       if ((jy.lt.10).and.(ix.lt.10)) write(*,*) 'reading: ' , ix, jy, lu_cat, lu_perc
      end do
    end do
  end do

  ! Read relation landuse,z0
  !*****************************

  open(unitsurfdata,file=path(1)(1:length(1))//'surfdata.t', &
       status='old',err=999)

  do i=1,4
    read(unitsurfdata,*)
  end do
  do i=1,numclass
    read(unitsurfdata,'(45x,f15.3)') z0(i)
  end do
  close(unitsurfdata)

  return

  ! Issue error messages
  !*********************

998   write(*,*) ' #### FLEXPART ERROR! FILE CONTAINING          ####'
  write(*,*) ' #### LANDUSE INVENTORY DOES NOT EXIST         ####'
  stop

999   write(*,*) ' #### FLEXPART ERROR! FILE CONTAINING          ####'
  write(*,*) ' #### RELATION LANDUSE,z0 DOES NOT EXIST       ####'
  stop

end subroutine readlanduse
