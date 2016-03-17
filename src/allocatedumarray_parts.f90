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
   

!*****************************************************************************
!                                                                            *
!     These routines allocates temporary variables for reading netcdf input  *
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
!                                                                            *
! Constants:                                                                 *
!                                                                            *
!*****************************************************************************
   
    
    
    subroutine  allocatedumarray1D(ndims,lendim_exp,maxdim)
      use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray1D)) deallocate(dumarray1D)
       allocate(dumarray1D(lendim_exp(1)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray2D(ndims,lendim_exp,maxdim)
      use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray2D)) deallocate(dumarray2D)
       allocate(dumarray2D(lendim_exp(1),lendim_exp(2)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray3D(ndims,lendim_exp,maxdim)
      use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray3D)) deallocate(dumarray3D)
       allocate(dumarray3D(lendim_exp(1),lendim_exp(2),lendim_exp(3)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
      
    subroutine  allocatedumarray4D(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray4D)) deallocate(dumarray4D)
       allocate(dumarray4D(lendim_exp(1),lendim_exp(2),lendim_exp(3),lendim_exp(4)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
      
      
      
        
    subroutine  allocatedumarray1D_real(ndims,lendim_exp,maxdim)
    use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
      if (allocated(dumarray1D_real)) deallocate(dumarray1D_real)
      allocate(dumarray1D_real(lendim_exp(1)))       
      if (allocated(istart)) deallocate(istart)
      allocate(istart(ndims)) 
      if (allocated(icount)) deallocate(icount)
      allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray2D_real(ndims,lendim_exp,maxdim)
    use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
      if (allocated(dumarray2D_real)) deallocate(dumarray2D_real)
      allocate(dumarray2D_real(lendim_exp(1),lendim_exp(2)))       
      if (allocated(istart)) deallocate(istart)
      allocate(istart(ndims)) 
      if (allocated(icount)) deallocate(icount)
      allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray3D_real(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray3D_real)) deallocate(dumarray3D_real)
       allocate(dumarray3D_real(lendim_exp(1),lendim_exp(2),lendim_exp(3)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
      
    subroutine  allocatedumarray4D_real(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray4D_real)) deallocate(dumarray4D_real)
       allocate(dumarray4D_real(lendim_exp(1),lendim_exp(2),lendim_exp(3),lendim_exp(4)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
      
    
      
    subroutine  allocatedumarray1D_int(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray1D_int)) deallocate(dumarray1D_int)
       allocate(dumarray1D_int(lendim_exp(1)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray2D_int(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray2D_int)) deallocate(dumarray2D_int)
       allocate(dumarray2D_int(lendim_exp(1),lendim_exp(2)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    
    subroutine  allocatedumarray3D_int(ndims,lendim_exp,maxdim)
    use noresm_variables
      integer maxdim
      integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray3D_int)) deallocate(dumarray3D_int)
       allocate(dumarray3D_int(lendim_exp(1),lendim_exp(2),lendim_exp(3)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    
    
    subroutine  allocatedumarray1D_char(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray1D_char)) deallocate(dumarray1D_char)
       allocate(dumarray1D_char(lendim_exp(1)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
    
    subroutine  allocatedumarray2D_char(ndims,lendim_exp,maxdim)
    use noresm_variables
       integer maxdim
       integer ndims,lendim_exp(maxdim)
      
       if (allocated(dumarray2D_char)) deallocate(dumarray2D_char)
       allocate(dumarray2D_char(lendim_exp(1),lendim_exp(2)))       
       if (allocated(istart)) deallocate(istart)
       allocate(istart(ndims)) 
       if (allocated(icount)) deallocate(icount)
       allocate(icount(ndims)) 
       
    return
    end
