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

      subroutine check_variable(varname,fnamenc,maxdim,typetocheck, &
      id_var,ndims,id_dim,ierr,ncid)
      implicit none
      include 'netcdf.inc'

!*****************************************************************************
!                                                                            *
!     Thies routines execture some consistency checks on the variable        * 
!     varname in a netcdf file  fnamec                                       *
!                                                                            *
!     Author:                                                                *
!     M. Cassiani  2016                                                      *
!                                                                            *
!                                                                            *
!*****************************************************************************

      
      character*160 :: varname,vartype
      integer :: maxdim
      integer :: id_var,itype_var,id_dim(maxdim),sizetype
      integer :: typetocheck 
      integer :: natts_tot, ncid, ndims
      integer :: ierr,iret !error code message
      character*110 :: fnamenc
!c--------- check  --------------------------------
      
          
      iret = nf_inq_varid( ncid, varname, id_var )
      if (iret .ne. nf_noerr) then
      write(*,9100) 'error inquiring var id for ' // varname, fnamenc
      ierr = -2 
      else 
      iret = nf_inq_var( ncid, id_var, &
      varname, itype_var, ndims, id_dim, natts_tot )
      if (iret .ne. nf_noerr) then
      write(*,9100) 'error inquiring var info for ' // varname, fnamenc
      ierr = -3
      else
      iret= nf_inq_type(ncid, itype_var, vartype, sizetype)
!c   check variable type
      if (itype_var .ne. typetocheck) then
      write(*,9110) 'var type wrong for ' // varname, &
      itype_var, fnamenc
      ierr = -4          
      end if
      end if
      end if
      
9100  format( / '*** read NorESM gridinfo -- ', a / &
      'file = ', a )
9110  format( / '*** read NorESM gridinfo -- ', a, 1x, i8 / &
      'file = ', a )
9120  format( / '*** read NorESM gridinfo -- ', a, 2(1x,i8) / &
      'file = ', a )
      
      return
      end
