! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine calcmatrix(lconv,delt,cbmf,metdata_format)
  !                        o    i    o
  !*****************************************************************************
  !                                                                            *
  !  This subroutine calculates the matrix describing convective               *
  !  redistribution of mass in a grid column, using the subroutine             *
  !  convect43c.f provided by Kerry Emanuel.                                   *
  !                                                                            *
  !  Petra Seibert, Bernd C. Krueger, 2000-2001                                *
  !                                                                            *
  !*****************************************************************************
  ! Changes:                                                                   *
  !  changed by C. Forster, November 2003 - February 2004                      *
  !  array fmassfrac(nconvlevmax,nconvlevmax) represents                       *
  !  the convective redistribution matrix for the particles                    *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Merged calcmatrix and calcmatrix_gfs into one routine using if-then  *
  !       for meteo-type dependent code                                        *
  !*****************************************************************************
  !                                                                            *
  !  lconv        indicates whether there is convection in this cell, or not   *
  !  delt         time step for convection [s]                                 *
  !  cbmf         cloud base mass flux                                         *
  !  metdata_format format of metdata (ecmwf/gfs)                              *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use conv_mod
  use class_gribfile

  implicit none

  real :: rlevmass,summe
  integer :: metdata_format

  integer :: iflag, k, kk, kuvz

  !1-d variables for convection
  !variables for redistribution matrix
  real :: cbmfold, precip, qprime
  real :: tprime, wd, f_qvsat
  real :: delt,cbmf
  logical :: lconv

  lconv = .false.


  ! calculate pressure at eta levels for use in convect
  ! and assign temp & spec. hum. to 1D workspace
  ! -------------------------------------------------------

  ! pconv(1) is the pressure at the first level above ground
  ! phconv(k) is the pressure between levels k-1 and k
  ! dpr(k) is the pressure difference "around" tconv(k)
  ! phconv(kmax) must also be defined 1/2 level above pconv(kmax)
  ! Therefore, we define k = kuvz-1 and let kuvz start from 2
  ! top layer cannot be used for convection because p at top of this layer is
  ! not given


  phconv(1) = psconv
  ! Emanuel subroutine needs pressure in hPa, therefore convert all pressures
  do kuvz = 2,nuvz
    k = kuvz-1
    if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    pconv(k) = (akz(kuvz) + bkz(kuvz)*psconv)
    phconv(kuvz) = (akm(kuvz) + bkm(kuvz)*psconv)
    else
      phconv(kuvz) =  0.5*(pconv(kuvz)+pconv(k))
    endif
    dpr(k) = phconv(k) - phconv(kuvz)
    qsconv(k) = f_qvsat( pconv(k), tconv(k) )

  ! initialize mass fractions
    do kk=1,nconvlev
      fmassfrac(k,kk)=0.
    end do
  end do


  !note that Emanuel says it is important
  !a. to set this =0. every grid point
  !b. to keep this value in the calling programme in the iteration

  ! CALL CONVECTION
  !******************

    cbmfold = cbmf
  ! Convert pressures to hPa, as required by Emanuel scheme
  !********************************************************
!!$    do k=1,nconvlev     !old
    do k=1,nconvlev+1      !bugfix
      pconv_hpa(k)=pconv(k)/100.
      phconv_hpa(k)=phconv(k)/100.
    end do
    phconv_hpa(nconvlev+1)=phconv(nconvlev+1)/100.
    call convect(nconvlevmax, nconvlev, delt, iflag, &
         precip, wd, tprime, qprime, cbmf)

  ! do not update fmassfrac and cloudbase massflux
  ! if no convection takes place or
  ! if a CFL criterion is violated in convect43c.f
   if (iflag .ne. 1 .and. iflag .ne. 4) then
     cbmf=cbmfold
     goto 200
   endif

  ! do not update fmassfrac and cloudbase massflux
  ! if the old and the new cloud base mass
  ! fluxes are zero
   if (cbmf.le.0..and.cbmfold.le.0.) then
     cbmf=cbmfold
     goto 200
   endif

  ! Update fmassfrac
  ! account for mass displaced from level k to level k

   lconv = .true.
    do k=1,nconvtop
      rlevmass = dpr(k)/ga
      summe = 0.
      do kk=1,nconvtop
        fmassfrac(k,kk) = delt*fmass(k,kk)
        summe = summe + fmassfrac(k,kk)
      end do
      fmassfrac(k,k)=fmassfrac(k,k) + rlevmass - summe
    end do

200   continue

end subroutine calcmatrix
