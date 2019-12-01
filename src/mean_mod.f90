module mean_mod
  public

! Interface to select single or double precision version of the 'mean'
! function from type of input arguments ("function overloading")
!********************************************************************
  interface mean
    module procedure mean_sp
    module procedure mean_dp
    module procedure mean_mixed_dss
    module procedure mean_mixed_dsd

  end interface mean

contains

  subroutine mean_sp(x_sp,xm,xs,number)

!*****************************************************************************
!                                                                            *
!  This subroutine calculates mean and standard deviation of a given element.*
!                                                                            *
!      AUTHOR: Andreas Stohl, 25 January 1994                                *
!                                                                            *
!      Single precision version ESO 2016                                     *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! x_sp(number)        field of input data                                    *
! xm                  mean                                                   *
! xs                  standard deviation                                     *
! number              number of elements of field x_sp                       *
!                                                                            *
! Constants:                                                                 *
! eps                 tiny number                                            *
!                                                                            *
!*****************************************************************************

    use par_mod, only: sp

    implicit none

    ! integer :: number,i
    ! real(sp) :: x_sp(number),xm,xs,xl,xq,xaux
    ! real(sp),parameter :: eps=1.0e-30

    integer,intent(in) :: number
    real(sp), intent(in) :: x_sp(number)
    real(sp), intent(out) ::xm,xs
    real(sp) :: xl,xq,xaux
    real(sp),parameter :: eps=1.0e-30
    integer :: i

    xl=0.
    xq=0.
    do i=1,number
      xl=xl+x_sp(i)
      xq=xq+x_sp(i)*x_sp(i)
    end do

    xm=xl/real(number,kind=sp)

    xaux=xq-xl*xl/real(number,kind=sp)

    if (xaux.lt.eps) then
      xs=0.
    else
      xs=sqrt(xaux/real(number-1,kind=sp))
    endif

  end subroutine mean_sp

  subroutine mean_dp(x_dp,xm,xs,number)

!*****************************************************************************
!                                                                            *
!  This subroutine calculates mean and standard deviation of a given element.*
!                                                                            *
!      AUTHOR: Andreas Stohl, 25 January 1994                                *
!                                                                            *
!      Double precision version ESO 2016                                     *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! x_dp(number)        field of input data                                    *
! xm                  mean                                                   *
! xs                  standard deviation                                     *
! number              number of elements of field x_dp                       *
!                                                                            *
! Constants:                                                                 *
! eps                 tiny number                                            *
!                                                                            *
!*****************************************************************************

    use par_mod, only: dp

    implicit none

    integer,intent(in) :: number
    real(dp), intent(in) :: x_dp(number)
    real(dp), intent(out) ::xm,xs
    real(dp) :: xl,xq,xaux
    real(dp),parameter :: eps=1.0e-30
    integer :: i

    xl=0._dp
    xq=0._dp
    do i=1,number
      xl=xl+x_dp(i)
      xq=xq+x_dp(i)*x_dp(i)
    end do

    xm=xl/real(number,kind=dp)

    xaux=xq-xl*xl/real(number,kind=dp)

    if (xaux.lt.eps) then
      xs=0._dp
    else
      xs=sqrt(xaux/real(number-1,kind=dp))
    endif

  end subroutine mean_dp

  subroutine mean_mixed_dss(x_dp,xm,xs,number)

!*****************************************************************************
!                                                                            *
!  This subroutine calculates mean and standard deviation of a given element.*
!                                                                            *
!      AUTHOR: Andreas Stohl, 25 January 1994                                *
!                                                                            *
!      Mixed precision version ESO 2016 (dp in, sp out, sp out)              *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! x_dp(number)        field of input data                                    *
! xm                  mean                                                   *
! xs                  standard deviation                                     *
! number              number of elements of field x_dp                       *
!                                                                            *
! Constants:                                                                 *
! eps                 tiny number                                            *
!                                                                            *
!*****************************************************************************

    use par_mod, only: sp,dp

    implicit none

    integer,intent(in) :: number
    real(dp), intent(in) :: x_dp(number)
    real(sp), intent(out) ::xm,xs
    real(sp) :: xl,xq,xaux
    real(sp),parameter :: eps=1.0e-30
    integer :: i

    xl=0._sp
    xq=0._sp
    do i=1,number
      xl=xl+real(x_dp(i),kind=sp)
      xq=xq+real(x_dp(i),kind=sp)*real(x_dp(i),kind=sp)
    end do

    xm=xl/real(number,kind=sp)

    xaux=xq-xl*xl/real(number,kind=sp)

    if (xaux.lt.eps) then
      xs=0._sp
    else
      xs=sqrt(xaux/real(number-1,kind=sp))
    endif

  end subroutine mean_mixed_dss

  subroutine mean_mixed_dsd(x_dp,xm,xs_dp,number)

!*****************************************************************************
!                                                                            *
!  This subroutine calculates mean and standard deviation of a given element.*
!                                                                            *
!      AUTHOR: Andreas Stohl, 25 January 1994                                *
!                                                                            *
!      Mixed precision version ESO 2016 (dp in, sp out, dp out)              *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! x_dp(number)        field of input data                                    *
! xm                  mean                                                   *
! xs_dp               standard deviation                                     *
! number              number of elements of field x_dp                       *
!                                                                            *
! Constants:                                                                 *
! eps                 tiny number                                            *
!                                                                            *
!*****************************************************************************

    use par_mod, only: sp,dp

    implicit none

    integer,intent(in) :: number
    real(dp), intent(in) :: x_dp(number)
    real(sp), intent(out) ::xm
    real(dp), intent(out) ::xs_dp
    real(dp) :: xl,xq,xaux
    real(dp),parameter :: eps=1.0e-30_dp
    integer :: i

    xl=0._dp
    xq=0._dp
    do i=1,number
      xl=xl+x_dp(i)
      xq=xq+x_dp(i)*x_dp(i)
    end do

    xm=real(xl,kind=sp)/real(number,kind=sp)

    xaux=xq-xl*xl/real(number,kind=dp)

    if (xaux.lt.eps) then
      xs_dp=0._dp
    else
      xs_dp=sqrt(xaux/real(number-1,kind=dp))
    endif

  end subroutine mean_mixed_dsd

end module mean_mod
