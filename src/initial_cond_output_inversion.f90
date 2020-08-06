! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine initial_cond_output_inversion(itime)
  !                                 i
  !*****************************************************************************
  !                                                                            *
  !     Output of the initial condition sensitivity field.                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: itime,i,ix,jy,kz,ks,kp,sp_count_i,sp_count_r
  integer :: jjjjmmdd, ihmmss
  real(kind=dp) :: jul
  real :: sp_fact,fact_recept
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  logical :: sp_zer,lexist
  logical,save :: listart=.true.
  logical,save,allocatable,dimension(:) :: listartrel
  character :: adate*8,atime*6
  character :: areldate*8,areltime*6
  character(len=3) :: anspec

  if(listart) then
    allocate(listartrel(maxpointspec_act))
    listartrel(:)=.true.
  endif
  print*, 'listartrel = ',listartrel

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

    write(anspec,'(i3.3)') ks

    do kp=1,maxpointspec_act

      ! calculate date of release
      jul=bdate+real(ireleasestart(kp),kind=dp)/86400._dp    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)
      write(areldate,'(i8.8)') jjjjmmdd
      write(areltime,'(i6.6)') ihmmss
      print*, areldate//areltime

      ! calculate date of field
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
      print*, adate//atime

      inquire(file=path(2)(1:length(2))//'grid_initial_'//areldate// &
         areltime//'_'//anspec,exist=lexist)
      if(lexist.and..not.listartrel(kp)) then
        ! open and append to existing file
        open(97,file=path(2)(1:length(2))//'grid_initial_'//areldate// &
             areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append') 
      else
        ! open new file
        open(97,file=path(2)(1:length(2))//'grid_initial_'//areldate// &
             areltime//'_'//anspec,form='unformatted',status='replace',action='write')
      endif
      write(97) jjjjmmdd
      write(97) ihmmss

      listartrel(kp)=.false.

      if (ind_rel.eq.1) then
        fact_recept=rho_rel(kp)
      else
        fact_recept=1.
      endif

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Write out dummy "wet and dry deposition" fields, to keep same format
  ! as for concentration output
!      sp_count_i=0
!      sp_count_r=0
!      write(97) sp_count_i
!      write(97) (sparse_dump_i(i),i=1,sp_count_i)
!      write(97) sp_count_r
!      write(97) (sparse_dump_r(i),i=1,sp_count_r)
!      write(97) sp_count_i
!      write(97) (sparse_dump_i(i),i=1,sp_count_i)
!      write(97) sp_count_r
!      write(97) (sparse_dump_r(i),i=1,sp_count_r)


  ! Write out sensitivity to initial conditions
      sp_count_i=0
      sp_count_r=0
      sp_fact=-1.
      sp_zer=.true.
      do kz=1,numzgrid
        do jy=0,numygrid-1
          do ix=0,numxgrid-1
            if (init_cond(ix,jy,kz,ks,kp).gt.smallnum) then
              if (sp_zer.eqv..true.) then ! first non zero value
                sp_count_i=sp_count_i+1
                sparse_dump_i(sp_count_i)= &
                     ix+jy*numxgrid+kz*numxgrid*numygrid
                sp_zer=.false.
                sp_fact=sp_fact*(-1.)
              endif
              sp_count_r=sp_count_r+1
              sparse_dump_r(sp_count_r)=sp_fact* &
                   init_cond(ix,jy,kz,ks,kp)/xmass(kp,ks)*fact_recept
            else ! concentration is zero
              sp_zer=.true.
            endif
          end do
        end do
      end do
      write(97) sp_count_i
      write(97) (sparse_dump_i(i),i=1,sp_count_i)
      write(97) sp_count_r
      write(97) (sparse_dump_r(i),i=1,sp_count_r)

      close(97)

    end do

  end do

  ! reset listart
  if (listart) then
    listart=.false.
  endif

end subroutine initial_cond_output_inversion
