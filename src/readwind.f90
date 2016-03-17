!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann, Ignacio Pisso                                      *
!                                                                     *
! This file is part of FLEXPART-NorESM                                *
!                                                                     *
! FLEXPART-NorESM is free software: you can redistribute it           *
! and/or modify                                                       *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART-NorESM is distributed in the hope that it will be useful,  *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART-NorESM.                                         *
!  If not, see <http://www.gnu.org/licenses/>.                        * 
!**********************************************************************
      


      Subroutine readwind(indj,n,uuh,vvh,wwh)


!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!*                                                                             *
!*     Read all fields needed to run FLEXPART                                  *
!*     note: part of netcdf reading strucutre adapted from the routines by     *
!*     Fast and Easter in FLEXPART-WRF                                         *
!*     if methodw is 1 it call a ruotine to obtain etadotdpdeta                * 
!*                                                                             *
!*                                                                             *
!*     Author:                                                                 *
!*     M. Cassiani  2016                                                       *
!*                                                                             *
!*                                                                             *
!*                                                                             *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


      
      use noresm_variables  
      use par_mod
      use com_mod   
      use conv_mod
      use cmapf_mod, only: stlmbr,stcm2p
      
      implicit none
      
      include 'netcdf.inc'
      
      !integer maxdim,maxvar,maxtime !moved to par_mod.f90
      real(kind=4) :: eps
      
      !parameter(maxdim=9,maxvar=60,maxtime=12, !moved in in common block par_mod
      parameter(eps=0.0001)      
      
      integer :: gotGrid     
      real(kind=4) :: sizesouth,sizenorth,xauxa,pint
    
!c---------------------------------------------
      real(kind=4) :: xaux1,xaux2,yaux1,yaux2
      real(kind=dp) :: xaux1in,xaux2in,yaux1in,yaux2in

      integer :: nyfieldin,nxfieldin
    
!C-------------- dichiarazioni per netcdf
      real(kind=4) :: dewpoint ! this is a function
      real(kind=4) :: duma,twomdewpoint
	  real(kind=4), allocatable, dimension(:) :: duma_alloc
      integer :: date_aid(maxtime),idumb,itime ,dimtimenum,indextime, &
      time_interval
      real(kind=dp) :: jul,juldate,dumb  !variable jul used for date in days, juldate is a function
      
      integer :: ierr    !error code message
      integer :: idiagaa !flag
      integer :: id_var,id_dim(maxdim)
      integer :: nvar_exp_in_grid_atm_nc, &
      nvar_exp_in_meteo_field
      
   
      integer :: i,ii,j,k, iatt, idimid_unlim, idum, iret, ivtype
      integer :: ix,jy,kz
      integer :: lenatt, lendim(maxdim)
      integer :: natts_tot, ncid, ndims_tot, nvars_tot
      integer :: sizetype
      character*110 :: fnamenc
      character*80 :: dimname(maxdim)
      character*80 :: attname
      character*160 :: varname,vartype
      character*160 :: varnamev(maxvar)
      character*160 :: units(maxvar)
      character*160 :: vartypev(maxvar)
      integer :: ndimsv(maxvar)
      integer :: dimidsv(maxvar,maxdim)
      character*1000 :: dumch1000
 
      !integer istart(maxdim),icount(maxdim)
      integer :: xtype,xtypev(maxvar)
      integer :: ndims
      integer :: dimids(maxdim),dimidsm(maxvar,maxdim)
      integer :: LENDIM_EXP(maxdim),LENDIM_MAX(maxdim)
      integer :: natts
      integer :: varid
      
      real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
      real(kind=4) :: nsss(0:nxmax-1,0:nymax-1),&
      ewss(0:nxmax-1,0:nymax-1)
       real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      logical :: hflswitch,strswitch
      integer :: indj,n !i,j,k,n,levdiff2,ifield,iumax,iwmax
      character*7 :: stringwftime
      

!***************************************************************************************************
9100  format( / '** read_noresmout_gridinfo -- ', a /&
      'file = ', a )
9110  format( / '** read_noresmout_gridinfo -- ', a, 1x, i8 /&
      'file = ', a )
9120  format( / '** read_noresmout_gridinfo -- ', a, 2(1x,i8) /&
      'file = ', a )

9030  format( a, 2i6, 2(2x,a) )
9031  format( 1i4,1a20,2i4,1a30)
   
     
      
      write(stringwftime,'(I7.7)')abs(wftime(indj))
      
      idiagaa=0 !diagnostic on reading and opening files if 1 write more info      
      if (idiagaa.eq.1) then 
      !open(71,file='..\windtrans\data_type.txt') !  for testing: mc
      !open(72,file='..\windtrans\list_windfield.txt') !  for testing: mc
      open(unitdiagnostic1,file='list_windfield.txt') !  for testing: mc
      !open(73,file='..\windtrans\list_global_att.txt') !  for testing: mc
      open(unitdiagnostic2,file='list_global_att.txt') !  for testing: mc
      !open(74,file='..\windtrans\list_variable_att.txt') !  for testing: mc
      open(unitdiagnostic3,file='list_variable_att.txt') !  for testing: mc
      !open(75,file='..\windtrans\seq_diagnostict.txt') !  for testing: mc
      !open(75,file='..\windtrans\seq_diagnostict.txt') !  for testing: mc
      !open(76,file='..\windtrans\vertical_wind_nov2_prima'//stringwftime//'.dat') !  for testing: mc
      !open(77,file='..\windtrans\vertical_wind_nov2_dopo'//stringwftime//'.dat') !  for testing: mc
      end if
       
      gotgrid=0
      ierr=0

!c-------------- open 4D (3D plus dtime dimension) wind meteo file 
      fnamenc=path(3)(1:length(3))&
      //trim(wfname(indj))
       
      iret = nf_open( fnamenc, NF_NOWRITE, ncid )
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error doing open', fnamenc
        ierr = -1
        stop
      end if
!c-------------- get information on dimension
      iret = nf_inq( ncid, &
      ndims_tot, nvars_tot, natts_tot, idimid_unlim )
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error inquiring dimensions', fnamenc
        ierr = -2
        stop
      end if
       
!c-------------inquiring abouot dimensions name -------
      dimtimenum=0
      do i = 1, min(ndims_tot,maxdim)
        iret = nf_inq_dim( ncid, i, dimname(i), lendim(i) )
        if (iret .ne. nf_noerr) then
          write(*,9110) 'error inquiring dimensions for dim#', i, fnamenc
          ierr = -2
          stop
        end if
        if (dimname(i).eq.'time') dimtimenum=i 
      end do
!c-----------------------------------------------------------      
!C------------inquiring about global attributes ---------

      if (idiagaa .gt. 0) then
        write(unitdiagnostic2,*)
        write(unitdiagnostic2,*) 'attribute #, name, type, value'
      end if
      do 3401 iatt = 1, natts_tot
        iret = nf_inq_attname( ncid, nf_global, iatt, attname)
        if (iret .ne. nf_noerr) goto 3601
        iret = nf_inq_att( ncid, nf_global, attname, ivtype, lenatt )
        if (iret .ne. nf_noerr) goto 3601
        if (ivtype .eq. 2) then
          iret = nf_get_att_text( ncid, nf_global, attname, dumch1000 )
          if (iret .ne. nf_noerr) goto 3601
          i = max(1,min(1000,lenatt))
                   
          if (idiagaa .gt. 0) write(unitdiagnostic2,91010) &
            iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
        else if (ivtype .eq. 4) then
          iret = nf_get_att_int( ncid, nf_global, attname, idum )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) write(unitdiagnostic2,91020) &
          iatt, attname(1:40), ivtype, lenatt, idum
        else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
          iret = nf_get_att_real( ncid, nf_global, attname, duma )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) write(unitdiagnostic2,91030) &
            iatt, attname(1:40), ivtype, lenatt, duma
        else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
          allocate( duma_alloc(lenatt) )
          iret = nf_get_att_real( ncid, nf_global, attname, duma_alloc )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) then
            write(unitdiagnostic2,91010) iatt, attname(1:40), ivtype, lenatt
            write(unitdiagnostic2,91040) (duma_alloc(i), i=1,lenatt)
           end if
          
          deallocate( duma_alloc )
        else
          if (idiagaa .gt. 0) write(unitdiagnostic2,'(i4,1x,a,2(1x,i6))') &
          iatt, attname(1:40), ivtype, lenatt
          goto 3401
        endif
    
3401  continue
91010 format( i4, 1x, a, 2(1x,i6), 1x, a )
91020 format( i4, 1x, a, 2(1x,i6), 1x, i10 )
91030 format( i4, 1x, a, 2(1x,i6), 1x, 1pe12.4 )
91040 format(( 12x, 5(1pe12.4) ))

91050 format(a,1x,2(1x,i5), 12(1x,e15.8))
      goto 3901

3601  write(*,9110) 'error inquiring global attribute', iatt, fnamenc
      !write(73,9110) 'error inquiring global attribute', iatt, fnamenc
      stop

3901  continue

    
  !c----------------some checks on the varibles--------------------    
      do ii=1, nvars_tot
        varid=ii 
        iret = nf_inq_var ( ncid, varid, varname, xtype, &
        ndims, dimids, natts)
        if (iret .ne. nf_noerr) then
          write (*,*) 'error in nf_inq_var'
          stop
        end if
  
        do iatt = 1, natts
          iret = nf_inq_attname( ncid, varid, iatt, attname)
          if (iret .ne. nf_noerr) goto 3602
          iret = nf_inq_att( ncid, varid, attname, ivtype, lenatt )
          if (iret .ne. nf_noerr) goto 3602
          if (ivtype .eq. 2) then
            dumch1000=''
            iret = nf_get_att_text( ncid, varid, attname, dumch1000 )
            if (iret .ne. nf_noerr) goto 3602
            i = max(1,min(1000,lenatt))
            if (attname.eq.'units') then 
              units(ii)=dumch1000
            end if
          
            if (idiagaa .gt. 0) write(unitdiagnostic3,91010) &
            iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
          else if (ivtype .eq. 4) then
            iret = nf_get_att_int( ncid, varid, attname, idum )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) write(unitdiagnostic3,91020) &
              iatt, attname(1:40), ivtype, lenatt, idum
          else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
            iret = nf_get_att_real( ncid, varid, attname, duma )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) write(unitdiagnostic3,91030) &
            iatt, attname(1:40), ivtype, lenatt, duma
          else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
            allocate( duma_alloc(lenatt) )
            iret = nf_get_att_real( ncid, varid, attname, duma_alloc )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) then
              write(unitdiagnostic3,91010) iatt, attname(1:40), ivtype, lenatt
               write(unitdiagnostic3,91040) (duma_alloc(i), i=1,lenatt)
            end if
          
            deallocate( duma_alloc )
          else
            if (idiagaa .gt. 0) write(unitdiagnostic3,'(i4,1x,a,2(1x,i6))') &
            iatt, attname(1:40), ivtype, lenatt
            exit
          endif
     
        end do
        goto 3902
3602    write(*,9110) 'error inquiring variable attribute', varname, &
        iatt, attname, fnamenc
        !write(74,9110) 'error inquiring variable attribute', varname, &
        !iatt, attname, fnamenc
        stop   
3902    continue
      
        !write(75,*)'start',ncid, varid, varname, xtype, & !more diagnostic for testing mc
        !ndims, dimids, natts,'fine'
  
        xtypev(ii)=xtype
        varnamev(ii)=varname
        ndimsv(ii)=ndims
        do j=1,maxdim
          dimidsm(ii,j)=dimids(j)
        end do
        !write(71,9031)i,varnamev(ii),xtypev(ii),ndims,units(ii) !WRITE diagnostic to file BY M C
      
      
      
      end do

      !C-------------- find the time --------------
          
      varname='date' !read date
      call check_variable(varname,fnamenc,maxdim,nf_int, &
      id_var,ndims,id_dim,ierr,ncid)     
      if (ierr.lt.0) goto 100
      lendim_exp=0
      if (ndims.eq.0) then  !this shoudl happen if there is only one time per file
        ndims=1
        lendim_exp(1)=1
      else
        do j = 1, ndims
          lendim_exp(j) = lendim(id_dim(j))            
        end do
      end if
      call allocatedumarray(ndims,lendim_exp,maxdim,nf_int) 
      do j=1, ndims
        istart(j) = 1 !ndims
        icount(j) = lendim_exp(j)
      end do 
      !here we assume taht ndims=1 for date variable
      iret = &
      nf_get_vara_int( ncid, id_var, istart, icount, dumarray1D_int)
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error inquiring var value ' // varname, fnamenc
        ierr = -5 
        goto 100 
      end if  
      do i=1,lendim_exp(1) !one dimension only expected!
          date_aid(i)=dumarray1D_int(i)
      end do
      
      
      varname='datesec'  !read date  in seconds
      call check_variable(varname,fnamenc,maxdim,nf_int, &
      id_var,ndims,id_dim,ierr,ncid)     
      if (ierr.lt.0) goto 100
        lendim_exp=0
      if (ndims.eq.0) then 
        ndims=1
        lendim_exp(1)=1
      else
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))            
        end do
      end if
      call allocatedumarray(ndims,lendim_exp,maxdim,nf_int) 
      do j=1, ndims
        istart(j) = 1 !ndims
        icount(j) = lendim_exp(j)
      end do 
      iret = &
      nf_get_vara_int( ncid, id_var, istart, icount, dumarray1D_int)
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error inquiring var value ' // varname, fnamenc
        ierr = -5 
        goto 100 
      end if  
      
      do i=1,lendim_exp(1) !one dimension only expected
          jul=juldate(date_aid(i),0)
          
          jul=jul+real(dumarray1D_int(i),kind=dp)/86400._dp
          dumb= (jul-bdate)*86400._dp
          idumb=nint(dumb)
          if (idumb.eq.wftime(indj)) exit
      end do     
      itime=i   !until here works.... it select right file and record....
      !time_interval=wftime(indj)-wftime(indj-1)
      !c--------  finish find time --------------
      if(idiagaa.eq.1) then
        write (unitdiagnostic1,*)wfname(indj),wftime(indj),itime,n  !write diagnostic to file for test reason by, by mc
      end if
      
      !c---------------- Inquiring about varnames and real fields ---------------------
      strswitch=.false.  ! shoud become .true. below if wind file has the right data
      hflswitch=.false.  ! shoud become .true. below if wind file has the right data
      nvar_exp_in_meteo_field=17
      varnamev(1)='U'  !horizontal wind    m/s
      varnamev(2)='V'   !horizontal wind   m/s 
      varnamev(3)='OMEGA'  ! vertical wind   Pa/s
      varnamev(4)='T'   !temperature   k
      varnamev(5)='Q'   !     Specific humidity  kg/kg
      varnamev(6)='PS'  !surface pressure    Pa
      varnamev(7)='CLDTOT'  !total cloud cover    fraction
      varnamev(8)='TREFHT'  !2m temperature    k
      varnamev(9)='PRECL'  !large scale precipitation m/s must become mm/hopur   in FLEXPARTECMWF ecmwf 142 but deaccumulated mm/hour /I must multiply by output time in s dvivide by hours and multipluy by 1000
      varnamev(10)='PRECC'  !convective precipitation  ms/ must become mm/hour   in FLEXPARTECMWF cmwf 143 but deaccumulated  mm/hour /I must multiply by output time in s dvivide by hours and multipluy by 1000
      varnamev(11)='SHFLX'  !sensible heat fluxes  W/m^2          in FLEXPARTECMWF ecmwf 146 but deaccumulated w/m^2 / this is already in Watt m^2 so nothing to do 
      varnamev(12)='TAUX'   !surface stress  east-west  N/m^2    in FLEXPARTECMWF ecmwf 180 but deaccumulated to N/m^2  / this is already in N/m^2 so nothign to do
      varnamev(13)='TAUY'   !surface stress  north-south N/m^2   in FLEXPARTECMWF ecmwf 181 but deaccumulated  to N/m^2 this is already in N/m^2 so nothing to do
      varnamev(14)='U10'     !ten meters wind speed m/s
      varnamev(15)='QREFHT'  !speciifc humidity 
      varnamev(16)='SNOWHLND' !m water equivakent snow depth !double check if must be also summed SNOWHICE 
      varnamev(17)='FSDS'   !downwelling solar flux atsurface! double checck if this correct for stomata opeining parameterizations !note that in ecmwf SSR ecmwf 176 but deaccumulated
      
      !---------------- check and load file contents    
      do i=1,nvar_exp_in_meteo_field !variable_loop 
        varname=varnamev(i)         
        call check_variable(varname,fnamenc,maxdim,nf_float, &
        id_var,ndims,id_dim,ierr,ncid)     
        if (ierr.ne.0) goto 100
        lendim_exp=0
        indextime=0
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))  
          if (id_dim(j).eq.dimtimenum) indextime=j !this indextime select whcih is the dimension associated with time
        end do
        if (dimtimenum.eq.0.and.indextime.eq.0) then !this means that there is no time dimension in the file and therefore we create an extra time dimension of leght 1.
          !notimedim=1 
          ndims=ndims+1
          lendim_exp(ndims)=1
          indextime=ndims
        end if 
      
        if (indextime.ne.ndims) stop 'ERROR readwind the time is expected &
      to be always the last dimension in the stored variables'
      
        call allocatedumarray(ndims,lendim_exp,maxdim,nf_float) 
       
        do j=1, ndims
          istart(j) = 1 !ndims
          icount(j) = lendim_exp(j)
        end do 
      
        if (ndims.eq.1) then
          iret = &
          nf_get_vara_real( ncid, id_var, istart, icount, dumarray1D_real)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        else if (ndims.eq.2) then       
          iret = &
          nf_get_vara_real( ncid, id_var, istart, icount, dumarray2D_real)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        else if (ndims.eq.3) then      
          iret = &
          nf_get_vara_real( ncid, id_var, istart, icount, dumarray3D_real)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5
            if ((varname.eq.'TAUx').or.(varname.eq.'TAUY').or.(varname.eq.&
            'SHFLX')) then
              print *,'************************************'
              print *,'* no stress available **************'
              print *,'************************************'
              cycle ! skyp the rest and go to the next variable
            end if
            goto 100 
          end if
        else if (ndims.eq.4) then      
          iret = &
          nf_get_vara_real( ncid, id_var, istart, icount, dumarray4D_real)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        end if
      !write(75,*) varname
      
        do jy=0,nymin1
          do ix=0, nxfield-1 
            
            if (varname.eq.'U') then
              do kz=1,nuvz-1 !recall that there is an augmented surface layer where U10 will be stored  
                uuh(ix,jy,nuvz+1-kz)=dumarray4D_real(ix+1,jy+1,kz,itime)
              end do
            else if (varname.eq.'V') then
              do kz=1,nuvz-1 !recall taht there is an augmented surface layer where V10 will be stored  
                vvh(ix,jy,nuvz+1-kz)=dumarray4D_real(ix+1,jy+1,kz,itime)
              end do  
            else if (varname.eq.'OMEGA') then
              do kz=1,nuvz-1 !recall that there is an augmented surface layer where W must be defined  
                wwh(ix,jy,nuvz+1-kz)=dumarray4D_real(ix+1,jy+1,kz,itime)
              end do  
            else if (varname.eq.'T') then
              do kz=1,nuvz-1 !recall that there is an augmented surface layer where T must be defined   
                tth(ix,jy,nuvz+1-kz,n)=dumarray4D_real(ix+1,jy+1,kz,itime) 
              end do
            else if (varname.eq.'Q') then
              do kz=1,nuvz-1 !recall that there is an augmented surface layer where Q must be defined   
                qvh(ix,jy,nuvz+1-kz,n)=dumarray4D_real(ix+1,jy+1,kz,itime) 
                if (qvh(ix,jy,nuvz+1-kz,n).lt.0) qvh(ix,jy,nuvz+1-kz,n)=0.
              end do
            else if (varname.eq.'PS') then 
              ps(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime) 
            else if (varname.eq.'CLDTOT') then 
              tcc(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)  
            else if (varname.eq.'TREFHT') then
              tt2(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)          
            else if (varname.eq.'PRECL') then
              lsprec(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)* &
              3600.*1000. !convet from m/s to mm/hour this must be controlled and eventually replaced with cumulated values
            else if (varname.eq.'PRECC') then 
              convprec(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)* &
              3600.*1000. !convet from m/s to mm/hour
            else if (varname.eq.'SHFLX') then
              sshf(ix,jy,1,n)=-dumarray3D_real(ix+1,jy+1,itime) ! note the minus sign here because FLEXPART assume positive FLUX downward like ECMWF: see e.g.http://tigge.ecmwf.int/tigge/d/show_object/table=parameters/name=time_integrated_surface_sensible_heat_flux/levtype=sfc/
              hflswitch=.true.    ! Heat flux available
            else if (varname.eq.'TAUX') then
              ewss(ix,jy)=dumarray3D_real(ix+1,jy+1,itime)
              strswitch=.true.  ! put this to true if stress are avilable
            else if (varname.eq.'TAUY') then
              nsss(ix,jy)=dumarray3D_real(ix+1,jy+1,itime)
            else if (varname.eq.'U10') then ! in NORESM U10 contain sped not component so it is partitioned using ground stresses, see below
              if  (strswitch) then
                u10(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)* ewss(ix,jy)/ &
                sqrt(ewss(ix,jy)**2+nsss(ix,jy)**2)  
                v10(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)* nsss(ix,jy)/ &
                sqrt(ewss(ix,jy)**2+nsss(ix,jy)**2)  
              else
                u10(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime) ! if stress not available then see below
                v10(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)
              end if
            else if (varname.eq.'QREFHT') then !
              td2(ix,jy,1,n)=dewpoint(tt2(ix,jy,1,n), &
              dumarray3D_real(ix+1,jy+1,itime),ps(ix,jy,1,n))
              qv2(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime)
            !else if (varname.eq.'SOLARRADIATION') then !
            ! ssr(ix,jy,1,n)=
            else if (varname.eq.'SNOWHLND') then 
              sd(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime) !m of water equivalent 
            else if (varname.eq.'FSDS') then
              ssr(ix,jy,1,n)=dumarray3D_real(ix+1,jy+1,itime) !downwelling solar radition flux in w/m^2
            end if
         
          end do
        end do
      
!c----------- deallocation added my mc on 12-2013 --------
        if (allocated(dumarray4D)) deallocate(dumarray4D)
        if (allocated(dumarray3D)) deallocate(dumarray3D)
        if (allocated(dumarray2D)) deallocate(dumarray2D)
        if (allocated(dumarray1D)) deallocate(dumarray1D)
        if (allocated(dumarray4D_real)) deallocate(dumarray4D_real)
        if (allocated(dumarray3D_real)) deallocate(dumarray3D_real)
        if (allocated(dumarray2D_real)) deallocate(dumarray2D_real)
        if (allocated(dumarray1D_real)) deallocate(dumarray1D_real)
        if (allocated(dumarray4D_int)) deallocate(dumarray4D_int)
        if (allocated(dumarray3D_int)) deallocate(dumarray3D_int)
        if (allocated(dumarray2D_int)) deallocate(dumarray2D_int)
        if (allocated(dumarray1D_int)) deallocate(dumarray1D_int)
  
!c-------------------------------------------------------------------------    
!c------------- read number of vertical level  
!c------------- note taht U, V, T, Q, OMEGA are colocated in CAM3.0/CAM4.0 see user's guide to thE NCAR CAM 3.0 page 38-
!c------------- so nwz and nuvz are teh same here while tadot levels form surface to top are nwz+1=27 for CAm 4.+0!
        if (varname.eq.'OMEGA') then
          if (nwz.ne.lendim_exp(3)) stop 'READWIND nwz NOT CONSISTENT' 
        else if (varname.eq.'U') then
          if (nuvz-1.ne.lendim_exp(3)) stop 'READWIND nwz NOT CONSITENT' 
        end if
!c----------         
      
      end do !variable_loop !
      
 !c---------------------------          
      iret = nf_close( ncid ) !closee netcdf files
      
  !c----------------------- more diagnostic for vertical velcoity ----------    
  !    do jy=0,nymin1
  !     do ix=0, nxfield-1 
  !       do kz=1,26
  !       write (76,'(3f12.4,4E15.6)')0+dx*ix,ylat0+dy*jy,kz*1.,wwh(ix,jy,kz)
  !       end do
  !     end do
  !    end do
      
!c----------------- here we have to put the call to the routine that will transform OMEGA in etadot * delta_p/delta_eta
      !deltat=wftime(indj+1)-wftime(indj-1) !time interval around indj
      if (method_w.eq.1) then ! <<<<---------------------if on method_w
        call transform_omega_etadot(uuh,vvh,wwh,n)
      end if
!C-----C more dignostic for vertical velocity ---------------   
      !do jy=0,nymin1
      ! do ix=0, nxfield-1 
      !   do kz=1,26
      !   write (77,'(3f12.4,4E15.6)')0+dx*ix,ylat0+dy*jy,kz*1.,wwh(ix,jy,kz)
      !   end do
      ! end do
      !end do
 !c-------------------------------------------------------------     

!c------------------------------------------------------------------------------------------------------------------
       
    
      
  !C----  If desired, shift all grids by nxshift grid cells ------
  
      if (xglobal) then
        call shift_field_0(ewss,nxfield,ny)
        call shift_field_0(nsss,nxfield,ny)
        !call shift_field_0(oro,nxfield,ny)
        !call shift_field_0(excessoro,nxfield,ny)
        !call shift_field_0(lsm,nxfield,ny)
        call shift_field(ps,nxfield,ny,1,1,2,n)
        call shift_field(sd,nxfield,ny,1,1,2,n)
        call shift_field(msl,nxfield,ny,1,1,2,n)
        call shift_field(tcc,nxfield,ny,1,1,2,n)
        call shift_field(u10,nxfield,ny,1,1,2,n)
        call shift_field(v10,nxfield,ny,1,1,2,n)
        call shift_field(tt2,nxfield,ny,1,1,2,n)
        call shift_field(td2,nxfield,ny,1,1,2,n)
        call shift_field(lsprec,nxfield,ny,1,1,2,n)
        call shift_field(convprec,nxfield,ny,1,1,2,n)
        call shift_field(sshf,nxfield,ny,1,1,2,n)
        call shift_field(ssr,nxfield,ny,1,1,2,n)
        call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
        call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
        call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)   
        call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
        call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
      endif
!C**********************************************

      do i=0,nxmin1
        do j=0,nymin1
          surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
        end do
      end do
      
 !c***************************************************** 
      if ((.not.hflswitch).or.(.not.strswitch)) then
        write(*,*) 'WARNING: No flux data contained in wind file ', &
        wfname(indj)

  ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
  ! As ECMWF has increased the model resolution, such that now the first model
  ! level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
  ! (3rd model level in FLEXPART) for the profile method
  !***************************************************************************

        do i=0,nxmin1
          do j=0,nymin1
            plev1=akz(3)+bkz(3)*ps(i,j,1,n)
            pmean=0.5*(ps(i,j,1,n)+plev1)
            tv=tth(i,j,3,n)*(1.+0.61*qvh(i,j,3,n))
            fu=-r_air*tv/ga/pmean
            hlev1=fu*(plev1-ps(i,j,1,n))   ! HEIGTH OF FIRST MODEL LAYER
            ff10m= u10(i,j,1,n) !note that  in NORESM, u10 and v10, if strswitch.eq.false, contain the speed not the components  
            fflev1=sqrt(uuh(i,j,3)**2+vvh(i,j,3)**2)
            call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
               tt2(i,j,1,n),tth(i,j,3,n),ff10m,fflev1, &
               surfstr(i,j,1,n),sshf(i,j,1,n))
            if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
            if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
          end do
        end do
        do i=0,nxmin1
          do j=0,nymin1
            u10(i,j,1,n)=u10(i,j,1,n)*uuh(i,j,2)/ & ! **NOTE:since we had no taux and tauy we used  speed at 10 meters to determine total stress see above 
            sqrt(uuh(i,j,2)**2+vvh(i,j,2)**2)         ! **therefore to assgin the wind component at 10 meters we have  to do a further assumption 
            v10(i,j,1,n)=v10(i,j,1,n)*vvh(i,j,2)/ & ! **this assumtpion is to take direction equal to the one of the first computed wind level 
            sqrt(uuh(i,j,2)**2+vvh(i,j,2)**2)         ! **which in general is probably not that correct..
          end do  
        end do    
       
      endif ! on surface stresses 
     
  ! Assign 10 m wind to model level at eta=1.0 to have one additional model
  ! level at the ground
  ! Specific humidity is taken the same as at one level above
  ! Temperature is taken as 2 m temperature
  !**************************************************************************

      do i=0,nxmin1
        do j=0,nymin1
           wwh(i,j,1)=0. !NORESM in ECMWf the OMEGA (hybrid system velocity in Pa/s) value at the ground is used  alternatively we shoudl use etadot (1/s) at the ground to obtain wwh at the ground
           uuh(i,j,1)=u10(i,j,1,n)
           vvh(i,j,1)=v10(i,j,1,n)
           !qvh(i,j,1,n)=qvh(i,j,2,n) !ECMWF version
           qvh(i,j,1,n)=qv2(i,j,1,n) ! this replace the above in NORESM version
           tth(i,j,1,n)=tt2(i,j,1,n)
        end do
      end do

       
      nlev_ec=nuvz-1 !in ECMWF version nelev_ec is the counter of ECMWF number of vertical levels for U (e.g. 91 with 92 levels)

      return
      
      
      
100   continue
      print *,'error reading',varname,'erros code',ierr
      
      stop
   
      
      return
      end

     
