!**********************************************************************
! Copyright 2016                                                      *
! Andreas Stohl, Massimo Cassiani, Petra Seibert, A. Frank,           *
! Gerhard Wotawa,  Caroline Forster, Sabine Eckhardt, John Burkhart,  *
! Harald Sodemann                                                     *
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

subroutine allocatedumarray(ndims,lendim_exp,maxdim,vartype)

  !*****************************************************************************
  !                                                                            *
  !     This routine allocate temporary variables for reading netcdf input     *
  !     files from NorESM                                                      *
  !                                                                            *
  !     Author:                                                                *
  !     M. Cassiani  2016                                                      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! vartype  type of variable                                                  *
  ! ndims number of dimensions 0,1,2,3                                         *
  ! lendim_exp  expeceted elemnts in the dimension                             *
  ! vartype                                                                    *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

     use noresm_variables
     implicit none
     include 'netcdf.inc'
     integer ndims,lendim_exp(maxdim),maxdim,vartype
     
     if (vartype.eq.nf_double) then
      if (ndims.eq.0) then
       dumvar=dble(0.0)
    
      else if (ndims.eq.1) then
       call allocatedumarray1D(ndims,lendim_exp,maxdim)    
      else if (ndims.eq.2) then
       call allocatedumarray2D(ndims,lendim_exp,maxdim)      
      else if (ndims.eq.3) then
       call allocatedumarray3D(ndims,lendim_exp,maxdim)        
      else if (ndims.eq.4) then
       call allocatedumarray4D(ndims,lendim_exp,maxdim)        
      end if 
     
     else if (vartype.eq.nf_float) then
      if (ndims.eq.0) then
       dumvar_real=0.0
      else if (ndims.eq.1) then
       call allocatedumarray1D_real(ndims,lendim_exp,maxdim)      
      else if (ndims.eq.2) then
       call allocatedumarray2D_real(ndims,lendim_exp,maxdim)    
      else if (ndims.eq.3) then
       call allocatedumarray3D_real(ndims,lendim_exp,maxdim)    
      else if (ndims.eq.4) then
       call allocatedumarray4D_real(ndims,lendim_exp,maxdim)       
      end if    
     
     else if (vartype.eq.nf_int) then
       if (ndims.eq.0) then
        dumvar_int=0.0
       else if (ndims.eq.1) then
        call allocatedumarray1D_int(ndims,lendim_exp,maxdim)      
      else if (ndims.eq.2) then
        call allocatedumarray2D_int(ndims,lendim_exp,maxdim)    
      else if (ndims.eq.3) then
        call allocatedumarray3D_int(ndims,lendim_exp,maxdim)    
     !else if (ndims.eq.4) then
     !      call allocatedumarray4D_int(ndims,lendim_exp,maxdim)       
      end if 
     
     else if (vartype.eq.nf_char) then
      if (ndims.eq.1) then
       call allocatedumarray1D_char(ndims,lendim_exp,maxdim)      
      else if (ndims.eq.2) then
       call allocatedumarray2D_char(ndims,lendim_exp,maxdim)    
      end if   
     end if
     return
     end
