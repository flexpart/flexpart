subroutine shift_field(field,nxf,nyf,nzfmax,nzf,nmax,n)
  !                        i/o   i   i    i     i   i   i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine shifts global fields by nxshift grid cells, in order to   *
  !  facilitate all sorts of nested wind fields, or output grids, which,       *
  !  without shifting, would overlap with the domain "boundary".               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    3 July 2002                                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  integer :: nxf,nyf,nzf,n,ix,jy,kz,ixs,nzfmax,nmax
  real :: field(0:nxmax-1,0:nymax-1,nzfmax,nmax),xshiftaux(0:nxmax-1)

  ! Loop over y and z
  !******************

  do kz=1,nzf
    do jy=0,nyf-1

  ! Shift the data
  !***************

      if (nxshift.ne.0) then
        do ix=0,nxf-1
          if (ix.ge.nxshift) then
            ixs=ix-nxshift
          else
            ixs=nxf-nxshift+ix
          endif
          xshiftaux(ixs)=field(ix,jy,kz,n)
        end do
        do ix=0,nxf-1
          field(ix,jy,kz,n)=xshiftaux(ix)
        end do
      endif

  ! Repeat the westernmost grid cells at the easternmost domain "boundary"
  !***********************************************************************

      field(nxf,jy,kz,n)=field(0,jy,kz,n)
    end do
  end do

end subroutine shift_field
