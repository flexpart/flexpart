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
   
      Subroutine gridcheck()
      
      use noresm_variables  
      use par_mod
      use com_mod
      use conv_mod
      use cmapf_mod, only: stlmbr,stcm2p
      
      implicit none
      
      include 'netcdf.inc'

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!*                                                                             *
!*     Read grid strucutre from netcdf files                                   *
!*     set grid for FLEXPART-NorESM                                            *
!*     check if variables are there in the netcdf files                        *
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

      
      !integer maxdim,maxvar,maxtime !moved in par_mod.f90
      real(kind=4) eps
      
      !parameter(maxdim=9,maxvar=70,maxtime=12,  !this must be moved in par_mod when finalized
      parameter(eps=0.0001)      
      
      integer :: gotGrid     
      real(kind=4) :: sizesouth,sizenorth,xauxa,pint
    
!c---------------------------------------------
      real(kind=4) :: xaux1,xaux2,yaux1,yaux2
      real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
      integer :: nyfieldin,nxfieldin
      
!C-------------- declaration for netcdf reading
     
      real(kind=4) :: duma
      real(kind=4), allocatable, dimension(:) :: duma_alloc

      integer :: ierr !error code message
      integer :: idiagaa !flag
      integer :: id_var,id_dim(maxdim)
      integer :: nvar_exp_in_grid_atm_nc,nvar_exp_in_meteo_field    
   
      integer :: i,j,iatt,idimid_unlim,idum,iret,ivtype
      integer :: ix,jy,ifn
      integer :: lenatt,lendim(maxdim)
      integer :: natts_tot,ncid,ndims_tot,nvars_tot
      integer :: sizetype
      character*110 :: fnamenc
      character*80 :: dimname(maxdim)
      character*80 :: attname
      character*160 :: varname,vartype
      character*160 :: varnamev(maxvar)
      character*160 :: vartypev(maxvar)

      integer :: ndimsv(maxvar)
      integer ::  dimidsv(maxvar,maxdim)
      character*1000 :: dumch1000

      !integer istart(maxdim),icount(maxdim)
      integer :: xtype,xtypev(maxvar)
      integer :: ndims
      integer :: dimids(maxdim),dimidsm(maxvar,maxdim)
      integer :: LENDIM_EXP(maxdim),LENDIM_MAX(maxdim)
      integer :: natts
      integer :: varid
      
!*************************************************************************************************      
      if(ideltas.gt.0) then
      ifn=1
      else
      ifn=numbwf
      endif

!*************************************************************************************************
9100  format( / '** read_noresmout_gridinfo -- ', a / &
        'file = ', a )
9110  format( / '** read_noresmout_gridinfo -- ', a, 1x, i8 / &
         'file = ', a )
9120  format( / '** read_noresmout_gridinfo -- ', a, 2(1x,i8) / &
         'file = ', a )

9030  format( a, 2i6, 2(2x,a) )
9031  format( 1i4,1a20,2i4)
   
 
      idiagaa=0  !diagnostic on reading and opening files if 1 write more info
      if (idiagaa.eq.1) then 
       !!open(71,file='..\options\data_type.txt') !  for testing: mc     
       !open(73,file='..\options\list_global_att.txt') !  for testing: mc
       open(unitdiagnostic2,file='list_global_att.txt') !  for testing: mc
       !!open(75,file='..\options\seq_diagnostict.txt') !  for testing: mc
      end if
      
      gotgrid=0
      ierr=0
!c-------------- open grid structure file
      
      fnamenc=path(5)(1:length(5)) !'..\..\..\data_set\grid_atm.nc' !grid_atm.nc' 
      ncid = 10
      iret = nf_open( fnamenc, NF_NOWRITE, ncid )
      if (iret .ne. nf_noerr) then
         write(*,9100) 'error doing open', fnamenc
         ierr = -1
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
      do i = 1, min(ndims_tot,maxdim)
         iret = nf_inq_dim( ncid, i, dimname(i), lendim(i) )
         if (iret .ne. nf_noerr) then
           write(*,9110) 'error inquiring dimensions for dim#', i, fnamenc
           ierr = -2
           stop
         end if
       end do
 
!c------------inquiring about global attributes ---------

       if (idiagaa .gt. 0) then
         write(73,*)
         write(73,*) 'grid strucuture: attribute #, name, type, value'
       end if
       do 3400 iatt = 1, natts_tot
         iret = nf_inq_attname( ncid, nf_global, iatt, attname)
         if (iret .ne. nf_noerr) goto 3600
         iret = nf_inq_att( ncid, nf_global, attname, ivtype, lenatt )
         if (iret .ne. nf_noerr) goto 3600
         if (ivtype .eq. 2) then
           iret = nf_get_att_text( ncid, nf_global, attname, dumch1000 )
           if (iret .ne. nf_noerr) goto 3600
           i = max(1,min(1000,lenatt))
           if (idiagaa .gt. 0) write(73,91010) &
               iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
         else if (ivtype .eq. 4) then
             iret = nf_get_att_int( ncid, nf_global, attname, idum )
             if (iret .ne. nf_noerr) goto 3600
             if (idiagaa .gt. 0) write(73,91020) &
               iatt, attname(1:40), ivtype, lenatt, idum
         else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
             iret = nf_get_att_real( ncid, nf_global, attname, duma )
             if (iret .ne. nf_noerr) goto 3600
             if (idiagaa .gt. 0) write(73,91030) &
               iatt, attname(1:40), ivtype, lenatt, duma
         else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
             allocate( duma_alloc(lenatt) )
             iret = nf_get_att_real( ncid, nf_global, attname, duma_alloc )
             if (iret .ne. nf_noerr) goto 3600
             if (idiagaa .gt. 0) then
               write(73,91010) iatt, attname(1:40), ivtype, lenatt
               write(73,91040) (duma_alloc(i), i=1,lenatt)
             end if          
             deallocate( duma_alloc )
         else
             if (idiagaa .gt. 0) write(73,'(i4,1x,a,2(1x,i6))')  &
                iatt, attname(1:40), ivtype, lenatt
             goto 3400
         endif
           
3400  continue
91010 format( i4, 1x, a, 2(1x,i6), 1x, a )
91020 format( i4, 1x, a, 2(1x,i6), 1x, i10 )
91030 format( i4, 1x, a, 2(1x,i6), 1x, 1pe12.4 )
91040 format(( 12x, 5(1pe12.4) ))

91050 format(a,1x,2(1x,i5), 12(1x,e15.8))

      goto 3900

3600  write(*,9110) 'error inquiring attribute', iatt, fnamenc
      write(73,9110) 'error inquiring attribute', iatt, fnamenc
      stop

3900  continue


!c---------------- Inquiring about varnames and real fields ---------------------      
      do i=1, nvars_tot
        varid=i
        iret = nf_inq_var ( ncid, varid, varname, xtype, &
        ndims, dimids, natts)
        if (iret .ne. nf_noerr) then
          write (*,*) 'error in nf_inq_var'
          stop
        end if
        !if (idiagaa.eq.1) then 
        !!write(75,*)'grid structure inquire',ncid, varid, varname, xtype, &
        !!ndims, dimids, natts,'fine'
        !end if
    
        xtypev(i)=xtype
        varnamev(i)=varname
        ndimsv(i)=ndims
        do j=1,maxdim
          dimidsm(i,j)=dimids(j)
        end do
        !if (idiagaa.eq.1) then 
        !!write(71,9031)i,varnamev(i),xtypev(i),ndims !WRITE diagnostic to file
        !end if
      end do      
! stop

!c---------- number of varibles to recover from netcdf file 
      nvar_exp_in_grid_atm_nc=6
!c--------- 1D fields       
      varnamev(1)='lat'       
      varnamev(2)='lon'  
!C--------- 2D fields       
      varnamev(3)='PHIS'     
      varnamev(4)='LANDFRAC'      
      varnamev(5)='SGH'
      varnamev(6)='xv'       
     
!c--------- read fields & set grid --------------------------------
      do i=1,nvar_exp_in_grid_atm_nc       
        varname=varnamev(i)         
        call check_variable(varname,fnamenc,maxdim,nf_double, &
          id_var,ndims,id_dim,ierr,ncid)     
        if (ierr.lt.0) goto 100
        lendim_exp=0
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))            
        end do
       
        call allocatedumarray(ndims,lendim_exp,maxdim,nf_double)  
        do j=1, ndims
          istart(j) = 1 !ndims
          icount(j) = lendim_exp(j)
        end do  
        if (ndims.eq.1) then
         iret = &
         nf_get_vara_double( ncid, id_var, istart, icount, dumarray1D)
         if (iret .ne. nf_noerr) then
         write(*,9100) 'error inquiring var value ' // varname, fnamenc
           ierr = -5 
           goto 100 
         end if
        else if (ndims.eq.2) then       
          iret = &
          nf_get_vara_double( ncid, id_var, istart, icount, dumarray2D)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        else if (ndims.eq.3) then       
          iret = &
          nf_get_vara_double( ncid, id_var, istart, icount, dumarray3D)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        end if
      
        !if (idiagaa.eq.1) then 
        !!write(75,*)'gridcheck',varname
        !end if
        
        !here load the arrays
        if (varname.eq.'lat') then
          nyfieldin=lendim_exp(1)
          yaux2in=dumarray1D(lendim_exp(1))  !last point 
          yaux1in=dumarray1D(1)              !first point
          ny=nyfieldin
        else if (varname.eq.'lon') then
          nxfieldin=lendim_exp(1)
          xaux2in=dumarray1D(lendim_exp(1))   !last point 
          xaux1in=dumarray1D(1)               !first point
          nxfield=nxfieldin
       
!c--------------    now set grid structure --------------
          xaux1=xaux1in
          xaux2=xaux2in
          yaux1=yaux1in
          yaux2=yaux2in
      
       
          if (xaux1.gt.180.) xaux1=xaux1-360.0
          if (xaux2.gt.180.) xaux2=xaux2-360.0
          if (xaux1.lt.-180.) xaux1=xaux1+360.0
          if (xaux2.lt.-180.) xaux2=xaux2+360.0
          if (xaux2.lt.xaux1) xaux2=xaux2+360.0
          xlon0=xaux1
          ylat0=yaux1
          !ny=nyfieldin
          dx=(xaux2-xaux1)/real(nxfield-1)
          dy=(yaux2-yaux1)/real(ny-1)
          dxconst=180./(dx*r_earth*pi)
          dyconst=180./(dy*r_earth*pi) 
          gotGrid=1

  ! Check whether fields are global
  ! If they contain the poles, specify polar stereographic map
  ! projections using the stlmbr- and stcm2p-calls
  !***********************************************************

          xauxa=abs(xaux2+dx-360.-xaux1)
          if (xauxa.lt.0.001) then
            nx=nxfield+1                 ! field is cyclic
            xglobal=.true.
            if (abs(nxshift).ge.nx) &
            stop 'nxshift in file par_mod is too large'
            xlon0=xlon0+real(nxshift)*dx
          else
            nx=nxfield
            xglobal=.false.
            if (xglobal.eqv..false.) &
            stop 'Noresm/CAM simualation supposed to have global wind field and &
      must use nxshift different from zero'
          endif
          nxmin1=nx-1
          nymin1=ny-1
          if (xlon0.gt.180.) xlon0=xlon0-360.
          xauxa=abs(yaux1+90.)
          if (xglobal.and.xauxa.lt.0.001) then
            sglobal=.true.               ! field contains south pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
            sizesouth=6.*(switchsouth+90.)/dy
            call stlmbr(southpolemap,-90.,0.)
            call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth, &
            sizesouth,switchsouth,180.)
            switchsouthg=(switchsouth-ylat0)/dy
          else
            sglobal=.false.
            switchsouthg=999999.
          endif
          xauxa=abs(yaux2-90.)
          if (xglobal.and.xauxa.lt.0.001) then
            nglobal=.true.               ! field contains north pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
            sizenorth=6.*(90.-switchnorth)/dy
            call stlmbr(northpolemap,90.,0.)
            call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth, &
             sizenorth,switchnorth,180.)
            switchnorthg=(switchnorth-ylat0)/dy
          else
            nglobal=.false.
            switchnorthg=999999.
          endif
          if (nxshift.lt.0) &
          stop 'nxshift (par_mod) must not be negative'
          if (nxshift.ge.nxfield) stop 'nxshift (par_mod) too large'

!c--------------- end of horizontal grid structure definitions -------- 
       
       
        else if (varname.eq.'PHIS') then
          do jy=0,ny-1
          do ix=0,nxfield-1
            oro(ix,jy)=sngl(dumarray2D(ix+1,jy+1))/ga
          end do
          end do
        else if (varname.eq.'SGH') then
          do jy=0,ny-1
          do ix=0,nxfield-1
            excessoro(ix,jy)=sngl(dumarray2D(ix+1,jy+1))
          end do
          end do
        else if (varname.eq.'LANDFRAC') then
          do jy=0,ny-1
          do ix=0,nxfield-1
            lsm(ix,jy)=sngl(dumarray2D(ix+1,jy+1))
           end do
           end do
        end if
       
      end do
!c------------------ set the grid ---------------------  

      print *,'gotgrid',gotgrid
      !error message if no fields found with correct first longitude in it
      if (gotGrid.eq.0) then
      print*,'***ERROR: horizontal grid strucutre not defined'
      stop
      endif
    
      ! goto 200
      iret = nf_close( ncid )
!c-------------- open 4D (3D plus time dimension) wind meteo file 


      fnamenc=path(3)(1:length(3))//trim(wfname(ifn)) 
                             
      iret = nf_open( fnamenc, NF_NOWRITE, ncid )
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error doing open', fnamenc
        ierr = -1
       !stop
      end if
!c-------------- get information on dimension
      iret = nf_inq( ncid, &
               ndims_tot, nvars_tot, natts_tot, idimid_unlim )
      if (iret .ne. nf_noerr) then
        write(*,9100) 'error inquiring dimensions', fnamenc
        ierr = -2
        stop
      end if
      idiagaa=1
!c-------------inquiring abouot dimensions name -------
      do i = 1, min(ndims_tot,maxdim)
        iret = nf_inq_dim( ncid, i, dimname(i), lendim(i) )
        if (iret .ne. nf_noerr) then
          write(*,9110) 'error inquiring dimensions for dim#',i,fnamenc
          ierr = -2
          stop
        end if
      end do
!c-----------------------------------------------------------      
!C------------inquiring about global attributes ---------

      if (idiagaa .gt. 0) then
        write(*,*)
        write(73,*)'gridcheck-windfiled: attribute #, name, type, value'
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
          allocate( duma_alloc(lenatt) )
          iret = nf_get_att_real( ncid, nf_global, attname, duma_alloc )
          if (iret .ne. nf_noerr) goto 3601
          if (idiagaa .gt. 0) then
            write(73,91010) iatt, attname(1:40), ivtype, lenatt
            write(73,91040) (duma_alloc(i), i=1,lenatt)
          end if
          
          deallocate( duma_alloc )
        else
          if (idiagaa .gt. 0) write(*,'(i4,1x,a,2(1x,i6))') &
              iatt, attname(1:40), ivtype, lenatt
          goto 3401
        endif
    
3401  continue

      goto 3901

3601  write(*,9110) 'gridcheck-windfiled:error inquiring attribute',iatt,fnamenc
      write(73,9110) 'gridcheck-windfiled:error inquiring attribute',iatt,fnamenc
      stop

3901  continue

    
  !c--------------------------------------------------    
      do i=1, nvars_tot
        varid=i
        iret = nf_inq_var ( ncid, varid, varname, xtype, &
        ndims, dimids, natts)
        if (iret .ne. nf_noerr) then
          write (*,*) 'error in nf_inq_var'
          stop
        end if
        !if (idiagaa.eq.1) then 
        !!write(75,*)'wind field inquire',ncid, varid, varname, xtype, &
        !!ndims, dimids, natts,'fine'
        !end if
  
        xtypev(i)=xtype
        varnamev(i)=varname
        ndimsv(i)=ndims
        do j=1,maxdim
          dimidsm(i,j)=dimids(j)
        end do
        !if (idiagaa.eq.1) then 
        !!write(71,9031)i,varnamev(i),xtypev(i),ndims !WRITE diagnostic to file
        !end if
      end do
      
      
      
      
      !c---------------- Inquiring about varnames and real fields ---------------------
      nvar_exp_in_meteo_field=17
      varnamev(1)='U'  !horizontal wind    
      varnamev(2)='V'   !horizontal wind    
      varnamev(3)='OMEGA'  ! vertical wind   pa/s
      varnamev(4)='T'   !temperature   
      varnamev(5)='Q'   !     Specific humidity 
      varnamev(6)='PS'  !surface pressure   
      varnamev(7)='CLDTOT'  !total cloud cover    
      varnamev(8)='U10'     !ten meters wind
      varnamev(9)='TREFHT'  !2m temperature    
      varnamev(10)='PRECL'  !large scale precipitation      
      varnamev(11)='PRECC'  !convective precipitation
      varnamev(12)='SHFLX'  !sensible heat fluxes    
      varnamev(13)='TAUX'   !surface stress  east-west
      varnamev(14)='TAUY'   !surface stress  north-south
      varnamev(15)='QREFHT'  !speciifc humidity 
      varnamev(16)='SNOWHLND' !water equivakent snow depth !double check if must be also summed SNOWHICE 
      varnamev(17)='FSDS'   !downwelling solar flux atsurface! double checck if this correct for stomata opeining parameterizations
      !---------------- check and load file contents    
      do i=1,nvar_exp_in_meteo_field 
        varname=varnamev(i)         
        call check_variable(varname,fnamenc,maxdim,nf_float, &
             id_var,ndims,id_dim,ierr,ncid)     
        if (ierr.lt.0) goto 100
        lendim_exp=0
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))            
        end do        
     
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
      !!write(75,*)'gridcheck read variable', varname
!c------------- read number of vertical level  
!c------------- seems that U, V, T, Q, OMEGA are co-located in CAM3.0/CAM4.0, see user's guide to thE NCAR CAM 3.0 page 38-
!c------------- so nwz and nuvz are the same here while etadot levels from surface to top are nwz+1=27 for CAM 4.+0!
        if (varname.eq.'OMEGA') then
          nwz=lendim_exp(3) 
        else if (varname.eq.'U') then
          nuvz=lendim_exp(3)
        end if
      
!c----------         
      
      end do
     
      
!200   continue

  !C If desired, shift all grids by nxshift grid cells
  !*************************************************

      if (xglobal) then
        call shift_field_0(oro,nxfield,ny)
        call shift_field_0(lsm,nxfield,ny)
        call shift_field_0(excessoro,nxfield,ny)
      endif

  ! Output of grid info
  !********************

      write(*,*)
      write(*,*)
      write(*,'(a,2i7)') '# of vertical levels in NORESM data: ', & 
      nuvz+1,nwz
      write(*,*)
      write(*,'(a)') 'Mother domain:'
      write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ', & 
      xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
      write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range: ', &
      ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy
      write(*,*)
      
!C-------------- load vertical grid strucutre --------------
      nvar_exp_in_meteo_field=6
      varnamev(1)='P0'   !omega levels 
      varnamev(2)='hyai'   !etadot levels
      varnamev(3)='hybi'   !etadot levels
      varnamev(4)='hyam'   !omega levels
      varnamev(5)='hybm'   !omega levels 
      varnamev(6)='time'   !omega levels 
      do i=1,nvar_exp_in_meteo_field 
        varname=varnamev(i)         
        call check_variable(varname,fnamenc,maxdim,nf_double, &
        id_var,ndims,id_dim,ierr,ncid)     
        if (ierr.lt.0) goto 100
        lendim_exp=0
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))            
        end do
      ! if (ndims.eq.0) then  !this alow teh use of a 1D array of 1 element for reading scalar variable (0-D)
      !  ndims=1
      !  lendim_exp(1)=1
      ! end if
        call allocatedumarray(ndims,lendim_exp,maxdim,nf_double) 
        do j=1, ndims
          istart(j) = 1 !ndims
          icount(j) = lendim_exp(j)
        end do 
            
        if (ndims.eq.0) then
          iret = &
          nf_get_vara_double( ncid, id_var, 1, 1, dumvar)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5
            goto 100 
          end if  
           
        else if (ndims.eq.1) then
          iret = &
          nf_get_vara_double( ncid, id_var, istart, icount, dumarray1D)
          if (iret .ne. nf_noerr) then
            write(*,9100) 'error inquiring var value ' // varname, fnamenc
            ierr = -5 
            goto 100 
          end if
        end if
       
        if (varname.eq.'hyai') then !etadot level 27 levels
          netadot=nwz+1 
          do j=1, nwz+1 
            akm(nwz+2-j)=sngl(dumarray1D(j)*dumvar) ! dumvar contains P0
          end do
        else if (varname.eq.'hybi') then  !etadot level 27 levels
          do j=1, nwz+1 
            bkm(nwz+2-j)=sngl(dumarray1D(j))
          end do
        else if (varname.eq.'hyam') then !omega and U,V  levels 26 levels, use 27 because added an "artificial" layer at the ground
          do j=1, nuvz 
            akz(nuvz+2-j)=sngl(dumarray1D(j)*dumvar) !dumvar contains P0
          end do
           akz(1)=0.
        else if (varname.eq.'hybm') then !omega and U,V  levels 26 levels, use 27 because added an "artificial" layer at the ground
          do j=1, nuvz
            bkz(nuvz+2-j)=sngl(dumarray1D(j)) 
          end do
          bkz(1)=1.0
        end if
      end do
      
!*************** Output of grid info!********************      

      write(*,*)
      write(*,*)
      write(*,'(a,2i7)') '# of vertical levels in NORESM data: ', &
       nuvz+1,nwz
      write(*,*)
      write(*,'(a)') 'Mother domain:'
      write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ', &
      xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
      write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range: ', &
      ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy
      write(*,*)
      nuvz=nuvz+1
      nz=nuvz
      
!c-- assign vertical discretization ------------------
      
      if (nz.gt.nzmax) stop 'nzmax too small'
      do i=1,nuvz
        aknew(i)=akz(i)
        bknew(i)=bkz(i)
      end do
   
!   THE DOUBLING IS NOT active in NORESM VERSION : Switch on following lines to use doubled vertical resolution
!   this is not active in FLEXPART-NorESM !!!! 
!  *************************************************************
!  nz=nuvz+nwz-1
!  if (nz.gt.nzmax) stop 'nzmax too small'
!  do 100 i=1,nwz
!    aknew(2*(i-1)+1)=akm(i)
!  100     bknew(2*(i-1)+1)=bkm(i)
!  do 110 i=2,nuvz
!    aknew(2*(i-1))=akz(i)
!  10     bknew(2*(i-1))=bkz(i)
!  !End doubled vertical resolution
      
!C-------------- load date --------------
      nvar_exp_in_meteo_field=1
      varnamev(1)='datesec'   !date
      
       
      do i=1,nvar_exp_in_meteo_field 
        varname=varnamev(i)         
        call check_variable(varname,fnamenc,maxdim,nf_int, &
        id_var,ndims,id_dim,ierr,ncid)     
        if (ierr.lt.0) goto 100
        lendim_exp=0
        do j = 1, ndims
          lendim_exp(j) =  lendim(id_dim(j))            
        end do
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
      end do
       
       iret = nf_close( ncid ) 

   !c----------- deallocation to free memory 12-2013 --------
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


  ! Determine the uppermost level for which the convection scheme shall be applied
  ! by assuming that there is no convection above 50 hPa (for standard SLP)
  !***************************************************************************** 
     
      do i=1,nuvz-2
        pint=akz(i)+bkz(i)*101325.
        if (pint.lt.5000.) goto 96
      end do
   96 nconvlev=i
      if (nconvlev.gt.nconvlevmax-1) then
        nconvlev=nconvlevmax-1
        write(*,*) 'Attention, convection only calculated up to ', &
        akz(nconvlev)+bkz(nconvlev)*101325,' Pa'
      endif

      return     
     
100   continue
      print *,'error reading',varname,'erros code',ierr
      stop
   
      return
      end    
