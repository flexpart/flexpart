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
     

      subroutine read_delta_ps_intime(indj,index)

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!*                                                                             *
!*     Read pressure fileds at time t-dt and t+dt                              *
!*     note: part of netcdf reading strucutre adapted from the routines by     *
!*     Fast and Easter in FLEXPART-WRF                                         *
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
    
!C-------------- extra function to calculate the dew point   ------
      real(kind=4) :: dewpoint ! this is a function

!c---------------varibles used to calculate the time interval --------- 
      integer :: date_aid(maxtime),idumb,itime ,dimtimenum,indextime, &
      time_interval
      real(kind=dp) :: jul,juldate,dumb  !variable jul used for date in days, juldate is a function
      
!C-------------- some declaration for netcdf reading ------ 
      real(kind=4) :: duma
      real(kind=4), allocatable, dimension(:) :: duma_alloc
      integer :: ierr !error code message
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
!C--------------------------------------------------------------------------------
      integer :: index !index of the time locaction with respect to wftime(indj) in readwind: 1 for indj-1 & 2 for indj+1      
         
      real(kind=4) :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      logical :: hflswitch,strswitch
      integer :: indj,n       

!***************************************************************************************************


9100  format( / '** read_deltap_intime -- ', a /&
      'file = ', a )
9110  format( / '** read_deltap_intime -- ', a, 1x, i8 /&
      'file = ', a )
9120  format( / '** read_deltap_intime -- ', a, 2(1x,i8) /&
      'file = ', a )
9030  format( a, 2i6, 2(2x,a) )
9031  format( 1i4,1a20,2i4,1a30)
   
      idiagaa=0 !set it to zero supress some diagnostic in files 73,74 see opening below!
      !if (idiagaa.eq.1) then 
      !open(71,file='..\options\data_type.txt') !  for testing: mc
      !open(72,file='..\options\list_windfield.txt') !  for testing: mc
      !open(73,file='..\options\list_global_att.txt') !  for testing: mc
      !open(74,file='..\options\list_variable_att.txt') !  for testing: mc
      !open(75,file='..\options\seq_diagnostict.txt') !  for testing: mc
      !continue
      !end if
      
      gotgrid=0
      ierr=0

!c-------------- open 4D (3D plus dtime dimension) wind meteo file 
      if (indj.le.1) goto 100
      fnamenc=path(3)(1:length(3))&
      //trim(wfname(indj))
       
      iret = nf_open(fnamenc,NF_NOWRITE,ncid)
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error doing open', fnamenc
        ierr = -1
        !stop
      end if
!c-------------- get information on dimensions ---------
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
     
!C------------inquiring about global attributes -------------

      if (idiagaa .gt. 0) then
        write(73,*)
        write(73,*) 'attribute #, name, type, value'
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
          if (idiagaa .gt. 0) write(73,91010) &
          iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
        else if (ivtype .eq. 4) then
          iret = nf_get_att_int( ncid, nf_global, attname, idum )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) write(73,91020) &
          iatt, attname(1:40), ivtype, lenatt, idum
        else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
          iret = nf_get_att_real( ncid, nf_global, attname, duma )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) write(73,91030) &
          iatt, attname(1:40), ivtype, lenatt, duma
        else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
          allocate(duma_alloc(lenatt))
          iret = nf_get_att_real( ncid, nf_global, attname, duma_alloc )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) then
            write(73,91010) iatt, attname(1:40), ivtype, lenatt
            write(73,91040) (duma_alloc(i), i=1,lenatt)
          end if        
          deallocate(duma_alloc)
        else
          if (idiagaa .gt. 0) write(73,'(i4,1x,a,2(1x,i6))') &
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
      write(73,9110) 'error inquiring global attribute', iatt, fnamenc
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
            if (idiagaa .gt. 0) write(74,91010) &
            iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
          else if (ivtype .eq. 4) then
            iret = nf_get_att_int( ncid, varid, attname, idum )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) write(74,91020) &
            iatt, attname(1:40), ivtype, lenatt, idum
          else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
            iret = nf_get_att_real( ncid, varid, attname, duma )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) write(74,91030) &
            iatt, attname(1:40), ivtype, lenatt, duma
          else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
            allocate( duma_alloc(lenatt) )
            iret = nf_get_att_real( ncid, varid, attname, duma_alloc )
            if (iret .ne. nf_noerr) goto 3602
            if (idiagaa .gt. 0) then
              write(74,91010) iatt, attname(1:40), ivtype, lenatt
              write(74,91040) (duma_alloc(i), i=1,lenatt)
            end if
            deallocate( duma_alloc )
          else
            if (idiagaa .gt. 0) write(74,'(i4,1x,a,2(1x,i6))') &
            iatt, attname(1:40), ivtype, lenatt
            exit
          endif
    
        end do
        goto 3902
3602    write(*,9110) 'error inquiring variable attribute', varname, &
        iatt, attname, fnamenc
        write(74,9110) 'error inquiring variable attribute', varname, &
        iatt, attname, fnamenc
        stop   
3902    continue
      
      !if (idiagaa.eq.1) then
      !write(75,*)'start',ncid, varid, varname, xtype, & !more diagnostic for testing mc
      !ndims, dimids, natts,'fine'
      !end if
  
        xtypev(ii)=xtype
        varnamev(ii)=varname
        ndimsv(ii)=ndims
        do j=1,maxdim
          dimidsm(ii,j)=dimids(j)
        end do
      
      !if (idiagaa.eq.1) then
      !write(71,9031)i,varnamev(ii),xtypev(ii),ndims,units(ii) !WRITE diagnostic to file BY M C
      !end if
      
      
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
      itime=i   !selected right file and record....
      time_interval=wftime(indj+1)-wftime(indj-1) !for calculating time derivative of pressure see below
      !c--------  finish operation of finding the time --------------
      print *,wfname(indj),wftime(indj),itime,n  !write diagnostic to file for test reason by, by mc
      
      !c---------------- Inquiring about varnames and real fields ---------------------
      nvar_exp_in_meteo_field=1
      varnamev(1)='PS'  !horizontal wind    m/s
       
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
      !!if (idiagaa.eq.1) then
      !!write(75,*) varname
      !!end if
      
        do jy=0,nymin1
          do ix=0, nxfield-1             
            if (varname.eq.'PS') then 
              ps_tplus1_and_min1(ix,jy,index)=dumarray3D_real(ix+1,jy+1,itime)/time_interval
            end if
          end do
        end do
     
      end do
      
      iret = nf_close( ncid )
     
      return
      
      
      
100   continue
      print *,'error reading for read_delta_ps_intime.f90',varname,'erros code',ierr,'or indj < 1:',indj
      print *,'delta_ps_in_time is set to be zero'
      do jy=0,nymin1
        do ix=0,nxfield-1
           ps_tplus1_and_min1(ix,jy,index)=0.
        end do
      end do
      
      
   
      
      return
      end subroutine read_delta_ps_intime

     
