! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine readdepo

  !*****************************************************************************
  !                                                                            *
  !  Reads dry deposition parameters needed by the procedure of Wesely (1989). *
  !  Wesely (1989): Parameterization of surface resistances to gaseous         *
  !  dry deposition in regional-scale numerical models.                        *
  !  Atmos. Environ. 23, 1293-1304.                                            *
  !                                                                            *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 19 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! rcl(maxspec,5,9) [s/m] Lower canopy resistance                             *
  ! rgs(maxspec,5,9) [s/m] Ground resistance                                   *
  ! rlu(maxspec,5,9) [s/m] Leaf cuticular resistance                           *
  ! rm(maxspec) [s/m]      Mesophyll resistance, set in readreleases           *
  ! ri(maxspec) [s/m]      Stomatal resistance                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  ! FOR THIS SUBROUTINE, numclass=9 IS ASSUMED
  !*******************************************

  real :: rluh(5,numclass),rgssh(5,numclass),rgsoh(5,numclass)
  real :: rclsh(5,numclass),rcloh(5,numclass)
  integer :: i,j,ic


  ! Read deposition constants related with landuse and seasonal category
  !*********************************************************************
  open(unitwesely,file=path(1)(1:length(1))//'surfdepo.t', &
       status='old',err=999)

  do i=1,16
    read(unitwesely,*)
  end do
  do i=1,5
    read(unitwesely,*)
    read(unitwesely,'(8x,13f8.0)') (ri(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rluh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rac(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rgssh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rgsoh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rclsh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rcloh(i,j),j=1,numclass)
  end do

  ! TEST
  ! do 31 i=1,5
  !   ri(i,13)=ri(i,5)
  !   rluh(i,13)=rluh(i,5)
  !   rac(i,13)=rac(i,5)
  !   rgssh(i,13)=rgssh(i,5)
  !   rgsoh(i,13)=rgsoh(i,5)
  !   rclsh(i,13)=rclsh(i,5)
  !   rcloh(i,13)=rcloh(i,5)
  !31             continue
  ! TEST
  ! Sabine Eckhardt, Dec 06, set resistances of 9999 to 'infinite' (1E25)
  do i=1,5
    do j=1,numclass
      if    (ri(i,j).eq.9999.)    ri(i,j)=1.E25
      if  (rluh(i,j).eq.9999.)  rluh(i,j)=1.E25
      if   (rac(i,j).eq.9999.)   rac(i,j)=1.E25
      if (rgssh(i,j).eq.9999.) rgssh(i,j)=1.E25
      if (rgsoh(i,j).eq.9999.) rgsoh(i,j)=1.E25
      if (rclsh(i,j).eq.9999.) rclsh(i,j)=1.E25
      if (rcloh(i,j).eq.9999.) rcloh(i,j)=1.E25
    end do
  end do



  do i=1,5
    do j=1,numclass
      ri(i,j)=max(ri(i,j),0.001)
      rluh(i,j)=max(rluh(i,j),0.001)
      rac(i,j)=max(rac(i,j),0.001)
      rgssh(i,j)=max(rgssh(i,j),0.001)
      rgsoh(i,j)=max(rgsoh(i,j),0.001)
      rclsh(i,j)=max(rclsh(i,j),0.001)
      rcloh(i,j)=max(rcloh(i,j),0.001)
    end do
  end do
  close(unitwesely)


  ! Compute additional parameters
  !******************************

  do ic=1,nspec
    if (reldiff(ic).gt.0.) then     ! gas is dry deposited
      do i=1,5
        do j=1,numclass
          rlu(ic,i,j)=rluh(i,j)/(1.e-5*henry(ic)+f0(ic))
          rgs(ic,i,j)=1./(henry(ic)/(10.e5*rgssh(i,j))+f0(ic)/ &
               rgsoh(i,j))
          rcl(ic,i,j)=1./(henry(ic)/(10.e5*rclsh(i,j))+f0(ic)/ &
               rcloh(i,j))
        end do
      end do
    endif
  end do


  return


999   write(*,*) '### FLEXPART ERROR! FILE              ###'
  write(*,*) '### surfdepo.t DOES NOT EXIST.        ###'
  stop

end subroutine readdepo
