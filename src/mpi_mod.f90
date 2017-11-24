!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

module mpi_mod

!*****************************************************************************
!                                                                            *
!  DESCRIPTION                                                               *
!    This module contains subroutines and common variables used for the      *
!    MPI parallelization of FLEXPART.                                        *
!                                                                            *
!  NOTE                                                                      *
!    Depending on the MPI library installed on your system (e.g. mpich2,     *
!    OpenMPI) you may need to choose below in this file between              *
!      use mpi                                                               *
!    (if the MPI library comes with the file 'mpi.mod'); or                  *
!      include 'mpif.h'                                                      *
!                                                                            *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
!                                                                            *
! mp_ierr                 MPI error code                                     *
! mp_np                   Number of MPI processes                            *
! mp_pid                  Process ID of each MPI process                     *
! mp_seed                 Parameter for random number seed                   *
! read_grp_min            Minimum number of processes at which one will be   *
!                         used as reader                                     *
! numpart_mpi,            Number of particles per node                       *
! maxpart_mpi                                                                *
! mp_partid               MPI process ID for particle calculation            *
! mp_partgroup_           Refers to the subset of processors performing      *
!                         loops over particles. Will be all processes        *
!                         unless a dedicated process runs getfields/readwind *
! lmp_sync                If .false., use asynchronous MPI                   *
! mp_cp                   Real precision to use for deposition fields        *
!                                                                            *
!                                                                            *
!                                                                            *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Modification by DJM, 2017-05-09 - added #ifdef USE_MPIINPLACE cpp          *
! directive to mpif_tm_reduce_grid() to insure that MPI_IN_PLACE is          *
! used in the MPI_Reduce() only if specifically compiled with that           *
! directive.                                                                 *
!                                                                            *
!*****************************************************************************

  use mpi 
  use par_mod, only: dp,sp
  use com_mod, only: lroot

  implicit none

!  include 'mpif.h'

  public

! Set aside a process for reading windfields if using at least these many processes
!==================================================
  integer, parameter, private :: read_grp_min=4
!==================================================

! Variables for each MPI process in the world group
  integer :: mp_ierr, mp_np, mp_pid, mp_partid
  integer, private :: world_group_id

! Variables for MPI processes in the 'particle' group
  integer, allocatable, dimension(:) :: mp_partgroup_rank
  integer :: mp_partgroup_comm, mp_partgroup_pid, mp_partgroup_np

  integer :: mp_seed=0
  integer, parameter :: mp_sp=MPI_REAL4, mp_dp=MPI_REAL8
  integer :: mp_cp
  integer, parameter :: id_root=0 ! master process

! MPI tags/requests for send/receive operation
  integer :: tm1
  integer, parameter :: nvar_async=26
!integer, dimension(:), allocatable :: tags
  integer, dimension(:), allocatable :: reqs

! Status array used for certain MPI operations (MPI_RECV)
  integer, dimension(MPI_STATUS_SIZE) :: mp_status


  integer :: id_read   ! readwind/getfield process
  integer :: numpart_mpi,maxpart_mpi ! number of particles per node
  integer :: tot_numpart=0
  integer :: mp_comm_used ! global or subgroup communicator

  logical :: lmpreader=.false. ! is set to true for reading process(es) only.
  logical :: lmp_use_reader=.false. ! true if separate readwind process is used

! .true. if only using synchronous MPI send/recv (default)
! If setting this to .false., numwfmem must be set to 3
!===============================================================================
  logical :: lmp_sync=.true. 
!===============================================================================

! mp_dbg_mode       Used for debugging MPI.
! mp_dev_mode       various checks related to debugging the parallel code
! mp_dbg_out        write some arrays to data file for debugging
! mp_measure_time   Measure cpu/wall time, write out at end of run
! mp_time_barrier   Measure MPI barrier time
! mp_exact_numpart  Use an extra MPI communication to give the exact number of particles
!                   to standard output (this does *not* otherwise affect the simulation) 
  logical, parameter :: mp_dbg_mode = .false.
  logical, parameter :: mp_dev_mode = .false.
  logical, parameter :: mp_dbg_out = .false.
  logical, parameter :: mp_time_barrier=.true.
  logical, parameter :: mp_measure_time=.false.
  logical, parameter :: mp_exact_numpart=.true.

! for measuring CPU/Wall time
  real(sp),private :: mp_comm_time_beg, mp_comm_time_end, mp_comm_time_total=0.
  real(dp),private :: mp_comm_wtime_beg, mp_comm_wtime_end, mp_comm_wtime_total=0.
  real(sp),private :: mp_root_time_beg, mp_root_time_end, mp_root_time_total=0.
  real(dp),private :: mp_root_wtime_beg, mp_root_wtime_end, mp_root_wtime_total=0.
  real(sp),private :: mp_barrier_time_beg, mp_barrier_time_end, mp_barrier_time_total=0.
  real(dp),private :: mp_barrier_wtime_beg, mp_barrier_wtime_end, mp_barrier_wtime_total=0.
  real(sp),private :: tm_nploop_beg, tm_nploop_end, tm_nploop_total=0.
  real(sp),private :: tm_tot_beg, tm_tot_end, tm_tot_total=0.
  real(dp),private :: mp_getfields_wtime_beg, mp_getfields_wtime_end, mp_getfields_wtime_total=0.
  real(sp),private :: mp_getfields_time_beg, mp_getfields_time_end, mp_getfields_time_total=0.
  real(dp),private :: mp_readwind_wtime_beg, mp_readwind_wtime_end, mp_readwind_wtime_total=0.
  real(sp),private :: mp_readwind_time_beg, mp_readwind_time_end, mp_readwind_time_total=0.
  real(dp),private :: mp_io_wtime_beg, mp_io_wtime_end, mp_io_wtime_total=0.
  real(sp),private :: mp_io_time_beg, mp_io_time_end, mp_io_time_total=0.
  real(dp),private :: mp_wetdepo_wtime_beg, mp_wetdepo_wtime_end, mp_wetdepo_wtime_total=0.
  real(sp),private :: mp_wetdepo_time_beg, mp_wetdepo_time_end, mp_wetdepo_time_total=0.
  real(dp),private :: mp_advance_wtime_beg, mp_advance_wtime_end, mp_advance_wtime_total=0.
  real(dp),private :: mp_conccalc_time_beg, mp_conccalc_time_end, mp_conccalc_time_total=0.
  real(dp),private :: mp_total_wtime_beg, mp_total_wtime_end, mp_total_wtime_total=0.
  real(dp),private :: mp_vt_wtime_beg, mp_vt_wtime_end, mp_vt_wtime_total
  real(sp),private :: mp_vt_time_beg, mp_vt_time_end, mp_vt_time_total

! dat_lun           logical unit number for i/o
  integer, private :: dat_lun 

! Temporary arrays for particles (allocated and deallocated as needed)
  integer, allocatable, dimension(:) :: nclass_tmp, npoint_tmp, itra1_tmp, idt_tmp, &
       & itramem_tmp, itrasplit_tmp
  real(kind=dp), allocatable, dimension(:) :: xtra1_tmp, ytra1_tmp
  real, allocatable, dimension(:) :: ztra1_tmp 
  real, allocatable, dimension(:,:) :: xmass1_tmp

! mp_redist_fract        Exchange particles between processes if relative difference
!                        is larger. A good value is between 0.0 and 0.5
! mp_maxpart_factor      Allocate more memory per process than strictly needed
!                        (safety factor). Recommended value between 1.5 and 2.5
! mp_min_redist          Do not redistribute particles if below this limit
  real, parameter :: mp_redist_fract=0.2, mp_maxpart_factor=1.5
  integer,parameter :: mp_min_redist=100000


contains

  subroutine mpif_init
!***********************************************************************
! DESCRIPTION
!   Initialize MPI.
!   
!   Create the global communicator MPI_COMM_WORLD
!   If using a separate MPI process for getfields/readwind, a subgroup
!   is created for the other processes.
!
! VARIABLES
!   mpi_mode    default 0, set to 2/3 if running MPI version
!   mp_np       number of running processes, decided at run-time
!***********************************************************************
    use par_mod, only: maxpart, numwfmem, dep_prec
    use com_mod, only: mpi_mode, verbosity

    implicit none

    integer :: i,j,s,addmaxpart=0

! Each process gets an ID (mp_pid) in the range 0,..,mp_np-1
    call MPI_INIT(mp_ierr)
    if (mp_ierr /= 0) goto 100
    call MPI_COMM_RANK(MPI_COMM_WORLD, mp_pid, mp_ierr)
    if (mp_ierr /= 0) goto 100
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mp_np, mp_ierr)
    if (mp_ierr /= 0) goto 100


! Variable mpi_mode is used to handle subroutines common to parallel/serial version
    if (lmp_sync) then
      mpi_mode=2 ! hold 2 windfields in memory
    else
      mpi_mode=3 ! hold 3 windfields in memory
    end if

    if (mp_pid.ne.0) then
      lroot = .false.
    end if

! Set MPI precision to use for transferring deposition fields
!************************************************************
    if (dep_prec==dp) then
      mp_cp = MPI_REAL8
! TODO: write info message for serial version as well 
      if (lroot.and.verbosity>0) write(*,*) 'Using double precision for deposition fields'
    else if (dep_prec==sp) then
      mp_cp = MPI_REAL4
      if (lroot.and.verbosity>0) write(*,*) 'Using single precision for deposition fields'
    else
      write(*,*) 'ERROR: something went wrong setting MPI real precision'
      stop
    end if

! Check for sensible combination of parameters
!*********************************************

    if (.not.lmp_sync.and.numwfmem.ne.3) then
      if (lroot) then
        write(*,FMT='(80("#"))')
        write(*,*) '#### mpi_mod::mpif_init> ERROR: ', &
             & 'numwfmem must be set to 3 for asyncronous reading ####'
        write(*,FMT='(80("#"))')
      end if
      call MPI_FINALIZE(mp_ierr)
      stop
    else if (lmp_sync.and.numwfmem.ne.2.and.lroot) then
      write(*,FMT='(80("#"))')
      write(*,*) '#### mpi_mod::mpif_init> WARNING: ', &
           & 'numwfmem should be set to 2 for syncronous'
      write(*,*) ' reading. Results will still be valid, but unneccesary memory &
           &is allocated.'
      write(*,FMT='(80("#"))')
! Force "syncronized" version if all processes will call getfields
    else if ((.not.lmp_sync.and.mp_np.lt.read_grp_min).or.(mp_np.eq.1)) then
      if (lroot) then
        write(*,FMT='(80("#"))')
        write(*,*) '#### mpi_mod::mpif_init> WARNING: ', &
             & 'all procs call getfields. Setting lmp_sync=.true.'
        write(*,FMT='(80("#"))')
      end if
      lmp_sync=.true.
    end if

! TODO: Add more warnings for unimplemented flexpart features

! Set ID of process that calls getfield/readwind. 
! Using the last process in the group ensures statistical identical results
! as running with one process less but not using separate read process
!**********************************************************************

!    id_read = min(mp_np-1, 1)
    id_read = mp_np-1

    if (mp_pid.eq.id_read) lmpreader=.true.

    call MPI_Comm_group (MPI_COMM_WORLD, world_group_id, mp_ierr)

! Create a MPI group/communicator that will calculate trajectories. 
! Skip this step if program is run with only a few processes
!************************************************************************
    allocate(mp_partgroup_rank(0:mp_np-2))

! This allows checking for allocation of arrays when no subgroup is used
    mp_partgroup_pid=0

    if (read_grp_min.lt.2) then
      write(*,*) '#### mpi_mod::mpif_init> ERROR ####', &
           & 'read_grp_min should be at least 2. Exiting'
      stop
    end if

    if (mp_np.ge.read_grp_min) then
      lmp_use_reader = .true.

! Map the subgroup IDs= 0,1,2,...,mp_np-2, skipping reader process
      j=-1 
      do i=0, mp_np-2 !loop over all (subgroup) IDs
        j=j+1
        if (i.eq.id_read) then
          j=j+1
        end if
        mp_partgroup_rank(i) = j
      end do

      call MPI_Group_incl (world_group_id, mp_np-1, mp_partgroup_rank, &
           &mp_partgroup_pid, mp_ierr)
      if (mp_ierr /= 0) goto 100
      call MPI_Comm_create (MPI_COMM_WORLD, mp_partgroup_pid, mp_partgroup_comm, mp_ierr)
      if (mp_ierr /= 0) goto 100

      if (mp_pid.ne.id_read) then
        call MPI_Comm_rank (mp_partgroup_comm, mp_partgroup_pid, mp_ierr)
        if (mp_ierr /= 0) goto 100

! Get size of new group (=mp_np-1)
        call MPI_COMM_SIZE(mp_partgroup_comm, mp_partgroup_np, mp_ierr)
        if (mp_ierr /= 0) goto 100
        if (mp_partgroup_np.ne.mp_np-1) then
          write(*,*)  '#### mpi_mod:: mpif_init> ERROR ####&
               & mp_partgroup_np.ne.mp_np-1'
          stop
        endif

      else
        mp_partgroup_pid = -1
      end if
    end if

    if (lmp_use_reader) then
      mp_comm_used = mp_partgroup_comm
      mp_partid = mp_partgroup_pid
      mp_partgroup_np=mp_np-1
    else
      mp_comm_used = MPI_COMM_WORLD
      mp_partgroup_np = mp_np
      mp_partid = mp_pid
    end if

! Set maxpart per process
! eso 08/2016: Increase maxpart per process, in case of unbalanced distribution
    maxpart_mpi=int(mp_maxpart_factor*real(maxpart)/real(mp_partgroup_np))
    if (mp_np == 1) maxpart_mpi = maxpart

! Do not allocate particle data arrays for readwind process
    if (lmpreader.and.lmp_use_reader) then
      maxpart_mpi=0
    end if

! Set random seed for each non-root process
    if (mp_pid.gt.0) then
      s = 244
      mp_seed = -abs(mod((s*181)*((mp_pid-83)*359), 104729))
    end if
    if (mp_dbg_mode) then
      mp_seed=0
      if (lroot) write(*,*) 'MPI: setting seed=0'
    end if

! Allocate request array for asynchronous MPI
    if (.not.lmp_sync) then
      allocate(reqs(0:nvar_async*mp_np-1))
      reqs(:)=MPI_REQUEST_NULL
    else
      allocate(reqs(0:1))
      reqs(:)=MPI_REQUEST_NULL
    end if

    goto 101

100 write(*,*) '#### mpi_mod::mpif_init> ERROR ####', mp_ierr
    stop

101 end subroutine mpif_init


  subroutine mpif_mpi_barrier
!***********************************************************************
! Set MPI barrier (all processes are synchronized here), with the
! possible exception of a separate 'readwind' process.
! Optionally measure cpu/wall time.
!
!***********************************************************************
    implicit none

    if (mp_time_barrier) then
      call cpu_time(mp_barrier_time_beg)
      mp_barrier_wtime_beg = mpi_wtime()
    endif

    call MPI_BARRIER(mp_comm_used, mp_ierr)
    if (mp_ierr /= 0) goto 100

    if (mp_time_barrier) then
      call cpu_time(mp_barrier_time_end)
      mp_barrier_wtime_end = mpi_wtime()

      mp_barrier_time_total=mp_barrier_time_total+&
           &(mp_barrier_time_end-mp_barrier_time_beg)
      mp_barrier_wtime_total=mp_barrier_wtime_total+&
           &(mp_barrier_wtime_end-mp_barrier_wtime_beg)
    end if

    goto 101

100 write(*,*) '#### mpi_mod::mpif_mpi_barrier> ERROR ####', mp_ierr
    stop

101 end subroutine mpif_mpi_barrier


  subroutine mpif_com_mod_allocate
!*******************************************************************************    
! Dynamic allocation of arrays (in serial code these have size maxpart)
!
! Not in use anymore, moved to com_mod for interoperability with serial version
! 
! TODO: remove 
!*******************************************************************************
    use com_mod

    implicit none 

    integer :: array_size

    array_size = maxpart_mpi

    allocate(itra1(array_size),npoint(array_size),&
         & nclass(array_size),idt(array_size),&
         & itramem(array_size),itrasplit(array_size),&
         & xtra1(array_size),ytra1(array_size),&
         & ztra1(array_size),xmass1(array_size, maxspec))

    allocate(uap(array_size),ucp(array_size),&
         & uzp(array_size),us(array_size),&
         & vs(array_size),&
         & ws(array_size),cbt(array_size))


  end subroutine mpif_com_mod_allocate


  subroutine mpif_tm_send_dims
!***********************************************************************
! Distribute array dimensions from pid0 to all processes.
! numpart_mpi/numpart is sent to allow dynamic allocation
!
! Currently not in use.
!
! Variables of similar type (integer, real) are placed in an array
! and sent collectively, to avoid the overhead associated with individual
! MPI send/receives
!
!
!***********************************************************************
    use com_mod

    implicit none

    integer :: add_part=0

    call MPI_Bcast(numpart, 1, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)

! MPI subgroup does the particle-loop 
    if (lmp_use_reader) then
      if (mod(numpart,mp_partgroup_np).ne.0) add_part=1
      numpart_mpi=int(numpart/mp_partgroup_np)+add_part
    else
      if (mod(numpart,mp_np).ne.0) add_part=1
      numpart_mpi=int(numpart/mp_np)+add_part
    end if


! redefine numpart as 'numpart per process' throughout the code
!**************************************************************
    numpart = numpart_mpi

  end subroutine mpif_tm_send_dims


  subroutine mpif_send_part_properties(num_part)
!***********************************************************************
! Distribute particle properties from root to all processes.
!  
! Used by init_domainfill_mpi
!
! Variables:
! 
! num_part        input, number of particles per process (rounded up)
!
!***********************************************************************
    use com_mod

    implicit none

    integer,intent(in) :: num_part
    integer :: i,jj, addone

! Exit if too many particles
    if (num_part.gt.maxpart_mpi) then
      write(*,*) '#####################################################'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES REQUIRED  ####'
      write(*,*) '#### EXCEEDS THE MAXIMUM ALLOWED NUMBER. REDUCE  ####'
      write(*,*) '#### EITHER NUMBER OF PARTICLES PER RELEASE POINT####'
      write(*,*) '#### OR INCREASE MAXPART.                        ####'
      write(*,*) '#####################################################'
!      call MPI_FINALIZE(mp_ierr)
      stop
    end if


! Time for MPI communications
!****************************
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Distribute variables, send from pid 0 to other processes (including itself)
!****************************************************************************

    call MPI_SCATTER(nclass_tmp,num_part,MPI_INTEGER,nclass,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(npoint_tmp,num_part,MPI_INTEGER,npoint,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(itra1_tmp,num_part,MPI_INTEGER,itra1,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(idt_tmp,num_part,MPI_INTEGER,idt,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(itramem_tmp,num_part,MPI_INTEGER,itramem,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(itrasplit_tmp,num_part,MPI_INTEGER,itrasplit,&
         &num_part,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(xtra1_tmp,num_part,mp_dp,xtra1,&
         &num_part,mp_dp,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(ytra1_tmp,num_part,mp_dp,ytra1,&
         &num_part,mp_dp,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_SCATTER(ztra1_tmp,num_part,mp_sp,ztra1,&
         &num_part,mp_sp,id_root,mp_comm_used,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    do i=1,nspec
      call MPI_SCATTER(xmass1_tmp(:,i),num_part,mp_sp,xmass1(:,i),&
           &num_part,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
    end do

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop

! After the transfer of particles it it possible that the value of
! "numpart" is larger than the actual used, so we reduce it if there are
! invalid particles at the end of the arrays

601 do i=num_part, 1, -1
      if (itra1(i).eq.-999999999) then
        numpart=numpart-1
      else
        exit
      end if
    end do


  end subroutine mpif_send_part_properties


  subroutine mpif_calculate_part_redist(itime)
!***********************************************************************
! Calculate number of particles to redistribute between processes. This routine 
! can be called at regular time intervals to keep a level number of
! particles on each process.
!
! First, all processes report their local "numpart" to each other, which is
! stored in array "numpart_mpi(np)". This array is sorted from low to
! high values, along with a corresponding process ID array "idx_arr(np)".
! If the relative difference between the highest and lowest "numpart_mpi" value
! is above a threshold, particles are transferred from process idx_arr(np-1)
! to process idx_arr(0). Repeat for processes idx_arr(np-i) and idx_arr(i)
! for all valid i. 
! Note: If np is an odd number, the process with the median
! number of particles is left unchanged
!
! VARIABLES
! itime        input, current time
!***********************************************************************
    use com_mod

    implicit none
    
    integer, intent(in) :: itime
    real :: pmin,z
    integer :: i,jj,nn, num_part=1,m,imin, num_trans
    logical :: first_iter
    integer,allocatable,dimension(:) :: numparticles_mpi, idx_arr
    real,allocatable,dimension(:) :: sorted ! TODO: we don't really need this

! Exit if running with only 1 process
!************************************************************************
    if (mp_np.eq.1) return

! All processes exchange information on number of particles
!****************************************************************************
    allocate(numparticles_mpi(0:mp_partgroup_np-1), &
         &idx_arr(0:mp_partgroup_np-1), sorted(0:mp_partgroup_np-1))

    call MPI_Allgather(numpart, 1, MPI_INTEGER, numparticles_mpi, &
         & 1, MPI_INTEGER, mp_comm_used, mp_ierr)


! Sort from lowest to highest
! Initial guess: correct order
    sorted(:) = numparticles_mpi(:)
    do i=0, mp_partgroup_np-1
      idx_arr(i) = i
    end do

! For each successive element in index array, see if a lower value exists
    do i=0, mp_partgroup_np-2
      pmin=sorted(i)
      imin=idx_arr(i)
      m=i+1
      do jj=m, mp_partgroup_np-1
        if (pmin.le.sorted(jj)) cycle
        z=pmin
        pmin=sorted(jj)
        sorted(jj)=z

        nn=imin
        imin=idx_arr(jj)
        idx_arr(jj)=nn
      end do
      sorted(i)=pmin
      idx_arr(i)=imin
    end do

! For each pair of processes, transfer particles if the difference is above a
! limit, and numpart at sending process large enough

    m=mp_partgroup_np-1 ! index for last sorted process (most particles)
    do i=0,mp_partgroup_np/2-1
      num_trans = numparticles_mpi(idx_arr(m)) - numparticles_mpi(idx_arr(i))
      if (mp_partid.eq.idx_arr(m).or.mp_partid.eq.idx_arr(i)) then
        if ( numparticles_mpi(idx_arr(m)).gt.mp_min_redist.and.&
             & real(num_trans)/real(numparticles_mpi(idx_arr(m))).gt.mp_redist_fract) then
! DBG
          ! write(*,*) 'mp_partid, idx_arr(m), idx_arr(i), mp_min_redist, num_trans, numparticles_mpi', &
          !      &mp_partid, idx_arr(m), idx_arr(i), mp_min_redist, num_trans, numparticles_mpi
! DBG
          call mpif_redist_part(itime, idx_arr(m), idx_arr(i), num_trans/2)
        end if
      end if
      m=m-1
    end do

    deallocate(numparticles_mpi, idx_arr, sorted)

  end subroutine mpif_calculate_part_redist


  subroutine mpif_redist_part(itime, src_proc, dest_proc, num_trans)
!***********************************************************************
! Transfer particle properties between two arbitrary processes.
! 
! VARIABLES
!
! itime           input, current time
! src_proc        input, ID of source process
! dest_proc       input, ID of destination process
! num_trans       input, number of particles to transfer
!
!************************************************************************
    use com_mod !TODO:  ,only: nclass etc

    implicit none 

    integer, intent(in) :: itime, src_proc, dest_proc, num_trans
    integer :: ll, ul ! lower and upper indices in arrays
    integer :: arr_sz ! size of temporary arrays
    integer :: mtag   ! MPI message tag
    integer :: i, j, minpart, ipart, maxnumpart
 
! Check for invalid input arguments
!**********************************
 if (src_proc.eq.dest_proc) then
   write(*,*) '<mpi_mod::mpif_redist_part>: Error: &
        &src_proc.eq.dest_proc' 
   stop
 end if

! Measure time for MPI communications
!************************************
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Send particles
!***************
    if (mp_partid.eq.src_proc) then
      mtag = 2000

! Array indices for the transferred particles
      ll=numpart-num_trans+1
      ul=numpart

      if (mp_dev_mode) then
        write(*,FMT='(72("#"))')
        write(*,*) "Sending ", num_trans, "particles (from/to)", src_proc, dest_proc
        write(*,*) "numpart before: ", numpart
      end if

      call MPI_SEND(nclass(ll:ul),num_trans,&
           & MPI_INTEGER,dest_proc,mtag+1,mp_comm_used,mp_ierr)

      call MPI_SEND(npoint(ll:ul),num_trans,&
           & MPI_INTEGER,dest_proc,mtag+2,mp_comm_used,mp_ierr)

      call MPI_SEND(itra1(ll:ul),num_trans, &
           & MPI_INTEGER,dest_proc,mtag+3,mp_comm_used,mp_ierr)

      call MPI_SEND(idt(ll:ul),num_trans, &
           & MPI_INTEGER,dest_proc,mtag+4,mp_comm_used,mp_ierr)

      call MPI_SEND(itramem(ll:ul),num_trans, &
           & MPI_INTEGER,dest_proc,mtag+5,mp_comm_used,mp_ierr)

      call MPI_SEND(itrasplit(ll:ul),num_trans,&
           & MPI_INTEGER,dest_proc,mtag+6,mp_comm_used,mp_ierr)

      call MPI_SEND(xtra1(ll:ul),num_trans, &
           & mp_dp,dest_proc,mtag+7,mp_comm_used,mp_ierr)

      call MPI_SEND(ytra1(ll:ul),num_trans,&
           & mp_dp,dest_proc,mtag+8,mp_comm_used,mp_ierr)

      call MPI_SEND(ztra1(ll:ul),num_trans,&
           & mp_sp,dest_proc,mtag+9,mp_comm_used,mp_ierr)

      do j=1,nspec
        call MPI_SEND(xmass1(ll:ul,j),num_trans,&
             & mp_sp,dest_proc,mtag+(9+j),mp_comm_used,mp_ierr)
      end do

! Terminate transferred particles, update numpart
      itra1(ll:ul) = -999999999

      numpart = numpart-num_trans

      if (mp_dev_mode) then
        write(*,*) "numpart after: ", numpart
        write(*,FMT='(72("#"))')
      end if

    else if (mp_partid.eq.dest_proc) then

! Receive particles
!******************
      mtag = 2000
      ! Array indices for the transferred particles
      ll=numpart+1
      ul=numpart+num_trans

      if (mp_dev_mode) then
        write(*,FMT='(72("#"))')
        write(*,*) "Receiving ", num_trans, "particles (from/to)", src_proc, dest_proc
        write(*,*) "numpart before: ", numpart
      end if

! We could receive the data directly at nclass(ll:ul) etc., but this leaves
! vacant spaces at lower indices. Using temporary arrays instead.
      arr_sz = ul-ll+1
      allocate(itra1_tmp(arr_sz),npoint_tmp(arr_sz),nclass_tmp(arr_sz),&
           & idt_tmp(arr_sz),itramem_tmp(arr_sz),itrasplit_tmp(arr_sz),&
           & xtra1_tmp(arr_sz),ytra1_tmp(arr_sz),ztra1_tmp(arr_sz),&
           & xmass1_tmp(arr_sz, maxspec))
      
      call MPI_RECV(nclass_tmp,num_trans,&
           & MPI_INTEGER,src_proc,mtag+1,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(npoint_tmp,num_trans,&
           & MPI_INTEGER,src_proc,mtag+2,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itra1_tmp,num_trans, &
           & MPI_INTEGER,src_proc,mtag+3,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(idt_tmp,num_trans, &
           & MPI_INTEGER,src_proc,mtag+4,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itramem_tmp,num_trans, &
           & MPI_INTEGER,src_proc,mtag+5,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(itrasplit_tmp,num_trans,&
           & MPI_INTEGER,src_proc,mtag+6,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(xtra1_tmp,num_trans, &
           & mp_dp,src_proc,mtag+7,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(ytra1_tmp,num_trans,&
           & mp_dp,src_proc,mtag+8,mp_comm_used,mp_status,mp_ierr)

      call MPI_RECV(ztra1_tmp,num_trans,&
           & mp_sp,src_proc,mtag+9,mp_comm_used,mp_status,mp_ierr)

      do j=1,nspec
        call MPI_RECV(xmass1_tmp(:,j),num_trans,&
             & mp_sp,src_proc,mtag+(9+j),mp_comm_used,mp_status,mp_ierr)
      end do

! This is the maximum value numpart can possibly have after the transfer
      maxnumpart=numpart+num_trans

! Search for vacant space and copy from temporary storage
!********************************************************
      minpart=1
      do i=1, num_trans
! Take into acount that we may have transferred invalid particles
        if (itra1_tmp(i).ne.itime) cycle
        do ipart=minpart,maxnumpart
          if (itra1(ipart).ne.itime) then
            itra1(ipart) = itra1_tmp(i)
            npoint(ipart) = npoint_tmp(i)
            nclass(ipart) = nclass_tmp(i)
            idt(ipart) = idt_tmp(i)
            itramem(ipart) = itramem_tmp(i)
            itrasplit(ipart) = itrasplit_tmp(i)
            xtra1(ipart) = xtra1_tmp(i)
            ytra1(ipart) = ytra1_tmp(i)
            ztra1(ipart) = ztra1_tmp(i)
            xmass1(ipart,:) = xmass1_tmp(i,:)
            goto 200 ! Storage space has been found, stop searching
          end if
! :TODO: add check for if too many particles requiried
        end do
200     minpart=ipart+1
      end do
! Increase numpart, if necessary
      numpart=max(numpart,ipart)

      deallocate(itra1_tmp,npoint_tmp,nclass_tmp,idt_tmp,itramem_tmp,&
           &itrasplit_tmp,xtra1_tmp,ytra1_tmp,ztra1_tmp,xmass1_tmp)

      if (mp_dev_mode) then
        write(*,*) "numpart after: ", numpart
        write(*,FMT='(72("#"))')
      end if

    else
! This routine should only be called by the two participating processes
      write(*,*) "ERROR: wrong process has entered mpi_mod::mpif_redist_part"
      stop
      return
    end if

! Measure time for MPI communications
!************************************
    if (mp_measure_time) call mpif_mtime('commtime',1)

  end subroutine mpif_redist_part


  subroutine mpif_tm_send_vars
!***********************************************************************
! Distribute particle variables from pid0 to all processes.
! Called from timemanager
! *NOT IN USE* at the moment, but can be useful for debugging
!
!***********************************************************************
    use com_mod

    implicit none

    integer :: i

! Time for MPI communications
!****************************
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Distribute variables, send from pid 0 to other processes
!**********************************************************************

! integers
    if (lroot) then
      call MPI_SCATTER(npoint,numpart_mpi,MPI_INTEGER,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(idt,numpart_mpi,MPI_INTEGER,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(itra1,numpart_mpi,MPI_INTEGER,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(nclass,numpart_mpi,MPI_INTEGER,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(itramem,numpart_mpi,MPI_INTEGER,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 

! int2
      call MPI_SCATTER(cbt,numpart_mpi,MPI_INTEGER2,MPI_IN_PLACE,&
           & numpart_mpi,MPI_INTEGER2,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

! real
      call MPI_SCATTER(uap,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(ucp,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(uzp,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

      call MPI_SCATTER(us,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(vs,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(ws,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

      call MPI_SCATTER(xtra1,numpart_mpi,mp_dp,MPI_IN_PLACE,&
           & numpart_mpi,mp_dp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(ytra1,numpart_mpi,mp_dp,MPI_IN_PLACE,&
           & numpart_mpi,mp_dp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(ztra1,numpart_mpi,mp_sp,MPI_IN_PLACE,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 

      do i=1,nspec
        call MPI_SCATTER(xmass1(:,i),numpart_mpi,mp_sp,MPI_IN_PLACE,&
             & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
        if (mp_ierr /= 0) goto 600 
      end do

    else ! (mp_pid >= 1)
! integers
      call MPI_SCATTER(npoint,numpart_mpi,MPI_INTEGER,npoint,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(idt,numpart_mpi,MPI_INTEGER,idt,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(itra1,numpart_mpi,MPI_INTEGER,itra1,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(nclass,numpart_mpi,MPI_INTEGER,nclass,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(itramem,numpart_mpi,MPI_INTEGER,itramem,&
           & numpart_mpi,MPI_INTEGER,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 

! int2
      call MPI_SCATTER(cbt,numpart_mpi,MPI_INTEGER2,cbt,&
           & numpart_mpi,MPI_INTEGER2,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

! reals
      call MPI_SCATTER(uap,numpart_mpi,mp_sp,uap,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(ucp,numpart_mpi,mp_sp,ucp,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(uzp,numpart_mpi,mp_sp,uzp,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

      call MPI_SCATTER(us,numpart_mpi,mp_sp,us,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(vs,numpart_mpi,mp_sp,vs,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_SCATTER(ws,numpart_mpi,mp_sp,ws,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600

      call MPI_SCATTER(xtra1,numpart_mpi,mp_dp,xtra1,&
           & numpart_mpi,mp_dp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(ytra1,numpart_mpi,mp_dp,ytra1,&
           & numpart_mpi,mp_dp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_SCATTER(ztra1,numpart_mpi,mp_sp,ztra1,&
           & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600 

      do i=1,nspec
        call MPI_SCATTER(xmass1(:,i),numpart_mpi,mp_sp,xmass1(:,i),&
             & numpart_mpi,mp_sp,id_root,mp_comm_used,mp_ierr)
        if (mp_ierr /= 0) goto 600 
      end do

    end if !(lroot)

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop
601 end subroutine mpif_tm_send_vars


  subroutine mpif_tm_collect_vars
!*******************************************************************************
! Collect particle data to PID 0 from all processes.                           *
!                                                                              *
! This can be called from timemanager_mpi, before caclulating                  *
! concentrations/writing output grids.                                         *
!                                                                              *
! Currently *not in use*, as each process calculates concentrations            *
! separately, but can be useful for debugging                                  *
!                                                                              *
! To use this routine (together with mpif_tm_send_vars) to transfer data       *
! to the MPI root process (useful for debugging), insert code like this        *
! (in timemanager_mpi.f90)                                                     *
!                                                                              *
!      if (lroot) tot_numpart=numpart ! root stores total numpart              *
!      call mpif_tm_send_dims                                                  *
!      if (numpart>1) then                                                     *
!        call mpif_tm_send_vars                                                *
!      end if                                                                  *
!                                                                              *
!    To collect data from all processes to root process after a parallel       *
!    region:                                                                   *
!                                                                              *
!      if (numpart>0) then                                                     *
!          if (lroot) numpart = tot_numpart                                    *
!       call mpif_tm_collect_vars                                              *
!      end if                                                                  *
!                                                                              *
! Then a section that begins with "if (lroot) ..." will behave like            *
! serial flexpart                                                              *
!                                                                              *
!*******************************************************************************
    use com_mod !, only: numpart, nspec, itra1, npoint, nclass

    implicit none

    integer :: i

    logical :: use_mp_vars = .false.! use mp_... type allocated vars
    logical :: use_in_place = .true.! use MPI_IN_PLACE to collect vars.


! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Distribute variables, send from pid 0 to other processes
!**********************************************************************

    if (.not. use_mp_vars.and..not.use_in_place) then

! Integers:
      call MPI_GATHER(npoint, numpart_mpi, MPI_INTEGER, npoint, &
           &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(idt, numpart_mpi, MPI_INTEGER, idt, &
           &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(itra1, numpart_mpi, MPI_INTEGER, itra1, &
           &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(nclass, numpart_mpi, MPI_INTEGER, nclass, &
           &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(itramem, numpart_mpi, MPI_INTEGER, itramem, &
           &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 

! int2
      call MPI_GATHER(cbt, numpart_mpi, MPI_INTEGER2, cbt, &
           &numpart_mpi, MPI_INTEGER2, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 

! Reals:
      call MPI_GATHER(uap, numpart_mpi, mp_sp, uap, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(ucp, numpart_mpi, mp_sp, ucp, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(uzp, numpart_mpi, mp_sp, uzp, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 

      call MPI_GATHER(us, numpart_mpi, mp_sp, us, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_GATHER(vs, numpart_mpi, mp_sp, vs, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(ws, numpart_mpi, mp_sp, ws, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600


      call MPI_GATHER(xtra1, numpart_mpi, mp_dp, xtra1, &
           &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(ytra1, numpart_mpi, mp_dp, ytra1, &
           &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_GATHER(ztra1, numpart_mpi, mp_sp, ztra1, &
           &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600 

      do i=1, nspec
        call MPI_GATHER(xmass1(:,i), numpart_mpi, mp_sp, xmass1(:,i), &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

      end do

! Use variables npoint etc directly for communications
!***********************************************************************
    else if (use_in_place.and..not.use_mp_vars) then
      if (lroot) then
! Integers:
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER, npoint, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER, idt, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER, itra1, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER, nclass, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER, itramem, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

! int2
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, MPI_INTEGER2, cbt, &
             &numpart_mpi, MPI_INTEGER2, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

! Reals:
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, uap, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, ucp, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, uzp, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, us, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, vs, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, ws, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600


        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_dp, xtra1, &
             &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_dp, ytra1, &
             &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, ztra1, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

        do i=1, nspec
          call MPI_GATHER(MPI_IN_PLACE, numpart_mpi, mp_sp, xmass1(:,i), &
               &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
          if (mp_ierr /= 0) goto 600 
        end do

      else ! (for mp_pid >= 1)

! Integers:
        call MPI_GATHER(npoint, numpart_mpi, MPI_INTEGER, npoint, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(idt, numpart_mpi, MPI_INTEGER, idt, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(itra1, numpart_mpi, MPI_INTEGER, itra1, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(nclass, numpart_mpi, MPI_INTEGER, nclass, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(itramem, numpart_mpi, MPI_INTEGER, itramem, &
             &numpart_mpi, MPI_INTEGER, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

! int2
        call MPI_GATHER(cbt, numpart_mpi, MPI_INTEGER2, cbt, &
             &numpart_mpi, MPI_INTEGER2, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

! Reals:
        call MPI_GATHER(uap, numpart_mpi, mp_sp, uap, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(ucp, numpart_mpi, mp_sp, ucp, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(uzp, numpart_mpi, mp_sp, uzp, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

        call MPI_GATHER(us, numpart_mpi, mp_sp, us, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600
        call MPI_GATHER(vs, numpart_mpi, mp_sp, vs, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(ws, numpart_mpi, mp_sp, ws, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600


        call MPI_GATHER(xtra1, numpart_mpi, mp_dp, xtra1, &
             &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(ytra1, numpart_mpi, mp_dp, ytra1, &
             &numpart_mpi, mp_dp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 
        call MPI_GATHER(ztra1, numpart_mpi, mp_sp, ztra1, &
             &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
        if (mp_ierr /= 0) goto 600 

        do i=1, nspec
          call MPI_GATHER(xmass1(:,i), numpart_mpi, mp_sp, xmass1(:,i), &
               &numpart_mpi, mp_sp, id_root, mp_comm_used, mp_ierr)
          if (mp_ierr /= 0) goto 600 
        end do
      end if ! (lroot)
    end if !  (use_in_place)

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop
601 end subroutine mpif_tm_collect_vars


  subroutine mpif_gf_send_vars(memstat)
!*******************************************************************************
! DESCRIPTION
!   Distribute 'getfield' variables from reader process
! 
!   Called from timemanager
!
! NOTE
!   This subroutine distributes windfields read from the reader process to
!   all other processes. Usually one set of fields are transfered, but at
!   first timestep when there are no fields in memory, both are transfered.
!   MPI_Bcast is used, so implicitly all processes are synchronized at this
!   step
!
!
!*******************************************************************************
    use com_mod
    use par_mod,only: numwfmem

    implicit none

    integer, intent(in) :: memstat

! Common array sizes used for communications
    integer, parameter :: d3_size1 = nxmax*nymax*nzmax
    integer, parameter :: d3_size2 = nxmax*nymax*nuvzmax
    integer, parameter :: d2_size1 = nxmax*nymax
    integer, parameter :: d2_size2 = nxmax*nymax*maxspec
    integer, parameter :: d2_size3 = nxmax*nymax
    integer, parameter :: d1_size1 = maxwf

    integer :: d3s1,d3s2,d2s1,d2s2

! li,ui    lower,upper indices of fields to be transfered (1-3, ui>=li) 
    integer :: li,ui

! First time routine is called the unchangeable fields will be transfered    
    logical :: first_call=.true.

!**********************************************************************

! Sizes of arrays transfered
    d3s1=d3_size1
    d3s2=d3_size2
    d2s1=d2_size1 
    d2s2=d2_size2

! Decide which field(s) need to be transfered
    if (memstat.ge.32) then ! distribute 2 fields, to 1st/2nd indices 
      li=1; ui=2
      d3s1=2*d3_size1
      d3s2=2*d3_size2
      d2s1=2*d2_size1 
      d2s2=2*d2_size2
    else if (memstat.eq.17) then ! distribute 1 field, place on 1st index
      li=1; ui=1
    else if (memstat.eq.18) then ! distribute 1 field, place on 2nd index
      li=2; ui=2
    else if (memstat.eq.19) then ! distribute 1 field, place on 3nd index
      li=3; ui=3
    else
      write(*,*) '#### mpi_mod::mpif_gf_send_vars> ERROR: ', &
           & 'wrong value of memstat, exiting ####', memstat
      stop
    end if


! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Send variables from getfield process (id_read) to other processes
!**********************************************************************

! The non-reader processes need to know if cloud water was read.
    call MPI_Bcast(readclouds,1,MPI_LOGICAL,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 

! Static fields/variables sent only at startup
    if (first_call) then
      call MPI_Bcast(oro(:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(excessoro(:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(lsm(:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(xlanduse(:,:,:),d2_size3*numclass,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(wftime,d1_size1,MPI_INTEGER,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(numbwf,1,MPI_INTEGER,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(nmixz,1,MPI_INTEGER,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(height,nzmax,MPI_INTEGER,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600

      first_call=.false.
    endif

    call MPI_Bcast(uu(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(vv(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(uupol(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(vvpol(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(ww(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(tt(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(rho(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(drhodz(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(tth(:,:,:,li:ui),d3s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(qvh(:,:,:,li:ui),d3s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(qv(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(pv(:,:,:,li:ui),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(clouds(:,:,:,li:ui),d3s1,MPI_INTEGER1,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600

! cloud water/ice:
    if (readclouds) then
      call MPI_Bcast(ctwc(:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if

! 2D fields
    call MPI_Bcast(cloudsh(:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(vdep(:,:,:,li:ui),d2s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(ps(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(sd(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(tcc(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(tt2(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(td2(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(lsprec(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(convprec(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(ustar(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(wstar(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(hmix(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Bcast(tropopause(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600
    call MPI_Bcast(oli(:,:,:,li:ui),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop

601 end subroutine mpif_gf_send_vars


  subroutine mpif_gf_send_vars_nest(memstat)
!***********************************************************************
! DESCRIPTION
!   Distribute 'getfield' variables from reader process to all processes.
!   For nested fields
! 
!   Called from timemanager
!
! NOTE
!   This subroutine distributes nested windfields from the reader process to
!   all other processes. Usually one set of fields is transfered, but at
!   first timestep when there are no fields in memory, both are transfered.
!   MPI_Bcast is used, so implicitly all processes are synchronized at this
!   step
!
!
!***********************************************************************
    use com_mod
    use par_mod,only: numwfmem

    implicit none

    integer, intent(in) :: memstat
    integer :: i
! li,ui    lower,upper indices of fields to be transfered (1-3, ui>=li) 
    integer :: li,ui

! Common array sizes used for communications
    integer :: d3_size1 = nxmaxn*nymaxn*nzmax
    integer :: d3_size2 = nxmaxn*nymaxn*nuvzmax
    integer :: d2_size1 = nxmaxn*nymaxn
    integer :: d2_size2 = nxmaxn*nymaxn*maxspec
    integer :: d2_size3 = nxmaxn*nymaxn

    integer :: d3s1,d3s2,d2s1,d2s2

! First time routine is called the unchangeable fields will be transfered    
    logical :: first_call=.true.

!**********************************************************************

! Sizes of arrays transfered
    d3s1=d3_size1
    d3s2=d3_size2
    d2s1=d2_size1 
    d2s2=d2_size2

! Decide which field(s) need to be transfered
    if (memstat.ge.32) then ! distribute 2 fields 
      li=1; ui=2
      d3s1=2*d3_size1
      d3s2=2*d3_size2
      d2s1=2*d2_size1 
      d2s2=2*d2_size2
    else if (memstat.eq.17) then ! distribute 1 field, on 1st index
      li=1; ui=1
    else if (memstat.eq.18) then ! distribute 1 field, on 2nd index
      li=2; ui=2
    else if (memstat.eq.19) then ! distribute 1 field, on 3nd index
      li=3; ui=3
    else
      write(*,*) '#### mpi_mod::mpif_gf_send_vars_nest> ERROR: ', &
           & 'wrong value of memstat, exiting ####', memstat
      stop
    end if


! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Distribute variables, send from getfield process (id_read) to other 
! processes
!**********************************************************************

! The non-reader processes need to know if cloud water was read.
    call MPI_Bcast(readclouds_nest,maxnests,MPI_LOGICAL,id_read,MPI_COMM_WORLD,mp_ierr)
    if (mp_ierr /= 0) goto 600 

! Static fields/variables sent only at startup
    if (first_call) then
      call MPI_Bcast(oron(:,:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(excessoron(:,:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(lsmn(:,:,:),d2_size3,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(xlandusen(:,:,:,:),d2_size3*numclass,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      first_call=.false.
    end if

! MPI prefers contiguous arrays for sending (else a buffer is created),
! hence the loop over nests
!**********************************************************************
    do i=1, numbnests 
! 3D fields
      call MPI_Bcast(uun(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(vvn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(wwn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(ttn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(rhon(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(drhodzn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(tthn(:,:,:,li:ui,i),d3s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(qvhn(:,:,:,li:ui,i),d3s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(qvn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(pvn(:,:,:,li:ui,i),d3s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(cloudsn(:,:,:,li:ui,i),d3s1,MPI_INTEGER1,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600

! cloud water/ice:
      if (readclouds_nest(i)) then
! call MPI_Bcast(icloud_stats(:,:,:,li:ui),d2s1*5,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
! if (mp_ierr /= 0) goto 600
        call MPI_Bcast(ctwcn(:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
        if (mp_ierr /= 0) goto 600
      end if

! 2D fields
      call MPI_Bcast(cloudshn(:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(vdepn(:,:,:,li:ui,i),d2s2,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(psn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(sdn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(tccn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600
      call MPI_Bcast(tt2n(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(td2n(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(lsprecn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(convprecn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(ustarn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(wstarn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(olin(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(hmixn(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Bcast(tropopausen(:,:,:,li:ui,i),d2s1,mp_sp,id_read,MPI_COMM_WORLD,mp_ierr)
      if (mp_ierr /= 0) goto 600 
    end do

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0",mp_ierr
    stop

601 end subroutine mpif_gf_send_vars_nest


  subroutine mpif_gf_send_vars_async(memstat)
!*******************************************************************************
! DESCRIPTION
!   Distribute 'getfield' variables from reader process to all processes.
!   Called from timemanager by reader process only
!
! NOTE
!   This version uses asynchronious sends. The newest fields are sent in the
!   background, so calculations can continue while
!   MPI communications are performed.
!
!   The way the MPI tags/requests are sequenced means that this routine must
!   carefully match routine 'mpif_gf_recv_vars_async'
!
! VARIABLES
!   memstat -- input, for resolving pointer to windfield index being read
!   mind    -- index where to place new fields
!
!
!*******************************************************************************
    use com_mod

    implicit none

    integer, intent(in) :: memstat
    integer :: mind
    integer :: dest,i

! Common array sizes used for communications
    integer :: d3s1 = nxmax*nymax*nzmax
    integer :: d3s2 = nxmax*nymax*nuvzmax
    integer :: d2s1 = nxmax*nymax
    integer :: d2s2 = nxmax*nymax*maxspec
!    integer :: d1s1 = maxwf

!*******************************************************************************

! At the time the send is posted, the reader process is one step further
! in the permutation of memstat compared with the receiving processes

    if (memstat.ge.32) then
! last read was synchronous, to indices 1 and 2, 3 is free
      write(*,*) "#### mpi_mod::mpif_gf_send_vars_async> ERROR: &
           & memstat>=32 should never happen here."
      stop
    else if (memstat.eq.17) then
! old fields on 1,2, send 3
      mind=3
    else if (memstat.eq.18) then
! old fields on 2,3, send 1
      mind=1
    else if (memstat.eq.19) then
! old fields on 3,1, send 2
      mind=2
    else
      write(*,*) "#### mpi_mod::mpif_gf_send_vars_async> ERROR: &
           & invalid memstat"
      mind=-1
      stop
    end if

    if (mp_dev_mode) then
      if (mp_pid.ne.id_read) then
        write(*,*) 'MPI_DEV_MODE: error calling mpif_gf_send_vars_async'
      end if
    end if

    if (mp_dev_mode) write(*,*) '## in mpif_gf_send_vars_async, memstat:', memstat 

! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Loop over receiving processes, initiate data sending
!*****************************************************

    do dest=0,mp_np-1 ! mp_np-2 will also work if last proc reserved for reading
! TODO: use mp_partgroup_np here
      if (dest.eq.id_read) cycle
      i=dest*nvar_async
      call MPI_Isend(uu(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(vv(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(uupol(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(vvpol(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(ww(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(tt(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(rho(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(drhodz(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(tth(:,:,:,mind),d3s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(qvh(:,:,:,mind),d3s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(qv(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(pv(:,:,:,mind),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(clouds(:,:,:,mind),d3s1,MPI_INTEGER1,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      i=i+1
      if (mp_ierr /= 0) goto 600 
      call MPI_Isend(cloudsh(:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(vdep(:,:,:,mind),d2s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
! 15
      call MPI_Isend(ps(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(sd(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(tcc(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(tt2(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(td2(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(lsprec(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(convprec(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(ustar(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(wstar(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
      call MPI_Isend(hmix(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      i=i+1
! 25
      call MPI_Isend(tropopause(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600
      i=i+1
      call MPI_Isend(oli(:,:,:,mind),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
      if (mp_ierr /= 0) goto 600 

! Send cloud water if it exists. Increment counter always (as on receiving end)
      if (readclouds) then
        i=i+1
        call MPI_Isend(ctwc(:,:,mind),d2s1,mp_sp,dest,tm1,&
             &MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
      end if
    end do

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "#### mpi_mod::mpif_gf_send_vars_async> mp_ierr \= 0", mp_ierr
    stop

601 end subroutine mpif_gf_send_vars_async


  subroutine mpif_gf_recv_vars_async(memstat)
!*******************************************************************************
! DESCRIPTION
!   Receive 'getfield' variables from reader process.
!   Called from timemanager by all processes except reader 
!
! NOTE
!   This version uses asynchronious communications.
!
! VARIABLES
!   memstat -- input, used to resolve windfield index being received
!
!
!*******************************************************************************
    use com_mod

    implicit none

    integer, intent(in) :: memstat
    integer :: mind,j

! Common array sizes used for communications
    integer :: d3s1 = nxmax*nymax*nzmax
    integer :: d3s2 = nxmax*nymax*nuvzmax
    integer :: d2s1 = nxmax*nymax
    integer :: d2s2 = nxmax*nymax*maxspec
!integer :: d1_size1 = maxwf

!    integer :: d3s1,d3s2,d2s1,d2s2
!*******************************************************************************

! At the time this immediate receive is posted, memstat is the state of
! windfield indices at the previous/current time. From this, the future
! state is deduced.
    if (memstat.eq.32) then
! last read was synchronous to indices 1 and 2, 3 is free
      mind=3
    else if (memstat.eq.17) then
! last read was asynchronous to index 3, 1 is free
      mind=1
    else if (memstat.eq.18) then
! last read was asynchronous to index 1, 2 is free
      mind=2
    else if (memstat.eq.19) then
! last read was asynchronous to index 2, 3 is free
      mind=3
    else
! illegal state
      write(*,*) 'mpi_mod> FLEXPART ERROR: Illegal memstat value. Exiting.'
      stop
    end if

    if (mp_dev_mode) then
      if (mp_pid.eq.id_read) then
        write(*,*) 'MPI_DEV_MODE: error calling mpif_gf_recv_vars_async'
      end if
    end if

! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

    if (mp_dev_mode) write(*,*) '## in mpif_gf_send_vars_async, memstat:', memstat 

! Initiate receiving of data 
!***************************

! Get MPI tags/requests for communications
    j=mp_pid*nvar_async
    call MPI_Irecv(uu(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(vv(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(uupol(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(vvpol(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(ww(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(tt(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(rho(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(drhodz(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(tth(:,:,:,mind),d3s2,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(qvh(:,:,:,mind),d3s2,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(qv(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(pv(:,:,:,mind),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(clouds(:,:,:,mind),d3s1,MPI_INTEGER1,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)    
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(cloudsh(:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(vdep(:,:,:,mind),d2s2,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(ps(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(sd(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(tcc(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(tt2(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(td2(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(lsprec(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(convprec(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    call MPI_Irecv(ustar(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(wstar(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(hmix(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 
    j=j+1
    call MPI_Irecv(tropopause(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600
    j=j+1
    call MPI_Irecv(oli(:,:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
         &MPI_COMM_WORLD,reqs(j),mp_ierr)
    if (mp_ierr /= 0) goto 600 

! Post request for clwc. These data are possibly not sent, request must then be cancelled
! For now assume that data at all steps either have or do not have water 
    if (readclouds) then
      j=j+1
      call MPI_Irecv(ctwc(:,:,mind),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if


    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "#### mpi_mod::mpif_gf_recv_vars_async> MPI ERROR ", mp_ierr
    stop

601 end subroutine mpif_gf_recv_vars_async


  subroutine mpif_gf_send_vars_nest_async(memstat)
!*******************************************************************************
! DESCRIPTION
!   Distribute 'getfield' variables from reader process to all processes.
!   Called from timemanager by reader process only.
!   Version for nested wind fields
!
! NOTE
!   This version uses asynchronious sends. The newest fields are sent in the
!   background, so calculations can continue while
!   MPI communications are performed.
!
!   The way the MPI tags/requests are sequenced means that this routine must
!   carefully match routine 'mpif_gf_recv_vars_async'
!
! VARIABLES
!   memstat -- input, for resolving pointer to windfield index being read
!   mind    -- index where to place new fields
!
! TODO
!   Some unused arrays are currently sent (uupoln,..) 
!*******************************************************************************
    use com_mod

    implicit none

    integer, intent(in) :: memstat
    integer :: mind
    integer :: dest,i,k

! Common array sizes used for communications
    integer :: d3s1 = nxmaxn*nymaxn*nzmax
    integer :: d3s2 = nxmaxn*nymaxn*nuvzmax
    integer :: d2s1 = nxmaxn*nymaxn
    integer :: d2s2 = nxmaxn*nymaxn*maxspec

!*******************************************************************************

! At the time the send is posted, the reader process is one step further
! in the permutation of memstat compared with the receiving processes

    if (memstat.ge.32) then
! last read was synchronous, to indices 1 and 2, 3 is free
      write(*,*) "#### mpi_mod::mpif_gf_send_vars_nest_async> ERROR: &
           & memstat>=32 should never happen here."
      stop
    else if (memstat.eq.17) then
! old fields on 1,2, send 3
      mind=3
    else if (memstat.eq.18) then
! old fields on 2,3, send 1
      mind=1
    else if (memstat.eq.19) then
! old fields on 3,1, send 2
      mind=2
    else
      write(*,*) "#### mpi_mod::mpif_gf_send_vars_nest_async> ERROR: &
           & invalid memstat"
      mind=-1
      stop
    end if

    if (mp_dev_mode) then
      if (mp_pid.ne.id_read) then
        write(*,*) 'MPI_DEV_MODE: error calling mpif_gf_send_vars_nest_async'
      end if
    end if

    if (mp_dev_mode) write(*,*) '## in mpif_gf_send_vars_nest_async, memstat:', memstat 

! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Loop over receiving processes, initiate data sending
!*****************************************************

    do dest=0,mp_np-1 ! mp_np-2 will also work if last proc reserved for reading
! TODO: use mp_partgroup_np here
      if (dest.eq.id_read) cycle
      do k=1, numbnests 
        i=dest*nvar_async
        call MPI_Isend(uun(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(vvn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(wwn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(ttn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(rhon(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(drhodzn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(tthn(:,:,:,mind,k),d3s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(qvhn(:,:,:,mind,k),d3s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(qvn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(pvn(:,:,:,mind,k),d3s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(cloudsn(:,:,:,mind,k),d3s1,MPI_INTEGER1,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        i=i+1
        if (mp_ierr /= 0) goto 600 
        call MPI_Isend(cloudshn(:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(vdepn(:,:,:,mind,k),d2s2,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(psn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(sdn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
! 15
        call MPI_Isend(tccn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(tt2n(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(td2n(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(lsprecn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(convprecn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(ustarn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(wstarn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(hmixn(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600 
        i=i+1
        call MPI_Isend(tropopausen(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
        i=i+1
        call MPI_Isend(olin(:,:,:,mind,k),d2s1,mp_sp,dest,tm1,MPI_COMM_WORLD,reqs(i),mp_ierr)
        if (mp_ierr /= 0) goto 600
! 25

! Send cloud water if it exists. Increment counter always (as on receiving end)
        if (readclouds) then
          i=i+1
          call MPI_Isend(ctwcn(:,:,mind,k),d2s1,mp_sp,dest,tm1,&
               &MPI_COMM_WORLD,reqs(i),mp_ierr)
          if (mp_ierr /= 0) goto 600
        end if
      end do
    end do

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "#### mpi_mod::mpif_gf_send_vars_nest_async> mp_ierr \= 0", mp_ierr
    stop

601 end subroutine mpif_gf_send_vars_nest_async


  subroutine mpif_gf_recv_vars_nest_async(memstat)
!*******************************************************************************
! DESCRIPTION
!   Receive 'getfield' variables from reader process.
!   Called from timemanager by all processes except reader
!   Version for nested wind fields
!
! NOTE
!   This version uses asynchronious communications.
!
! VARIABLES
!   memstat -- input, used to resolve windfield index being received
!
!
!*******************************************************************************
    use com_mod

    implicit none

    integer, intent(in) :: memstat
    integer :: mind,j,k

! Common array sizes used for communications
    integer :: d3s1 = nxmaxn*nymaxn*nzmax
    integer :: d3s2 = nxmaxn*nymaxn*nuvzmax
    integer :: d2s1 = nxmaxn*nymaxn
    integer :: d2s2 = nxmaxn*nymaxn*maxspec

!*******************************************************************************

! At the time this immediate receive is posted, memstat is the state of
! windfield indices at the previous/current time. From this, the future
! state is deduced.
    if (memstat.eq.32) then
! last read was synchronous to indices 1 and 2, 3 is free
      mind=3
    else if (memstat.eq.17) then
! last read was asynchronous to index 3, 1 is free
      mind=1
    else if (memstat.eq.18) then
! last read was asynchronous to index 1, 2 is free
      mind=2
    else if (memstat.eq.19) then
! last read was asynchronous to index 2, 3 is free
      mind=3
    else
! illegal state
      write(*,*) 'mpi_mod> FLEXPART ERROR: Illegal memstat value. Exiting.'
      stop
    end if

    if (mp_dev_mode) then
      if (mp_pid.eq.id_read) then
        write(*,*) 'MPI_DEV_MODE: error calling mpif_gf_recv_vars_async'
      end if
    end if

! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

    if (mp_dev_mode) write(*,*) '## in mpif_gf_send_vars_async, memstat:', memstat 

! Initiate receiving of data 
!***************************

    do k=1, numbnests 
! Get MPI tags/requests for communications
      j=mp_pid*nvar_async
      call MPI_Irecv(uun(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(vvn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(wwn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(ttn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(rhon(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(drhodzn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(tthn(:,:,:,mind,k),d3s2,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(qvhn(:,:,:,mind,k),d3s2,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(qvn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(pvn(:,:,:,mind,k),d3s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(cloudsn(:,:,:,mind,k),d3s1,MPI_INTEGER1,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)    
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(cloudshn(:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(vdepn(:,:,:,mind,k),d2s2,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(psn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(sdn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(tccn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(tt2n(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(td2n(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(lsprecn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(convprecn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      call MPI_Irecv(ustarn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(wstarn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(hmixn(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 
      j=j+1
      call MPI_Irecv(tropopausen(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600
      j=j+1
      call MPI_Irecv(olin(:,:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
           &MPI_COMM_WORLD,reqs(j),mp_ierr)
      if (mp_ierr /= 0) goto 600 

! Post request for clwc. These data are possibly not sent, request must then be cancelled
! For now assume that data at all steps either have or do not have water 
      if (readclouds) then
        j=j+1
        call MPI_Irecv(ctwcn(:,:,mind,k),d2s1,mp_sp,id_read,MPI_ANY_TAG,&
             &MPI_COMM_WORLD,reqs(j),mp_ierr)
        if (mp_ierr /= 0) goto 600
      end if
    end do

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "#### mpi_mod::mpif_gf_recv_vars_nest_async> MPI ERROR ", mp_ierr
    stop

601 end subroutine mpif_gf_recv_vars_nest_async


  subroutine mpif_gf_request
!*******************************************************************************
! DESCRIPTION
!   Check for completion of MPI_Isend/MPI_Irecv data transfer.
!   
!   
! NOTE
!   implicit synchronisation between all processes takes place here
!
! TODO
!   NB: take into account nested wind fields by using a separate 
!   subroutine mpif_gf_request_nest (and separate request arrays for the
!   nested variables)
!
!*******************************************************************************
!    use com_mod,only: readclouds

    implicit none


    integer :: n_req !,i

!***********************************************************************

    n_req=nvar_async*mp_np

    if (mp_measure_time) call mpif_mtime('commtime',0)

!    call MPI_Wait(rm1,MPI_STATUS_IGNORE,mp_ierr)

! TODO: cancel recv request if readclouds=.false.
!    if (readclouds) then
    call MPI_Waitall(n_req,reqs,MPI_STATUSES_IGNORE,mp_ierr)
!    endif
! else
!   do i = 0, nvar_async*mp_np-1
!     if (mod(i,27).eq.0 .or. mod(i,28).eq.0) then
!       call MPI_Cancel(reqs(i),mp_ierr)
!       cycle
!     end if
!     call MPI_Wait(reqs(i),MPI_STATUS_IGNORE,mp_ierr)
!   end do
! end if

    if (mp_ierr /= 0) goto 600 

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "#### mpi_mod::mpif_gf_request> MPI ERROR ", mp_ierr
    stop

601 end subroutine mpif_gf_request


  subroutine mpif_tm_reduce_grid
!***********************************************************************
! Collect grid variable to PID 0, adding from all processes.
!
! NOTE
!   - Older MPI implementations may lack the MPI_IN_PLACE operation.
!     If so use 1) below and change unc_mod
!
!***********************************************************************
    use com_mod
    use unc_mod
    use par_mod

    implicit none

    integer :: grid_size2d,grid_size3d
    integer, parameter :: rcpt_size=maxreceptor*maxspec

!**********************************************************************
    grid_size2d=numxgrid*numygrid*maxspec* &
         & maxpointspec_act*nclassunc*maxageclass
    grid_size3d=numxgrid*numygrid*numzgrid*maxspec* &
         & maxpointspec_act*nclassunc*maxageclass


! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! 1) Using a separate grid (gridunc0) for received values
! call MPI_Reduce(gridunc, gridunc0, grid_size3d, mp_sp, MPI_SUM, id_root, &
!      & mp_comm_used, mp_ierr)
! if (mp_ierr /= 0) goto 600

! 2) Using in-place reduction

!!!--------------------------------------------------------------------
!!! DJM - 2017-05-09 change - MPI_IN_PLACE option for MPI_Reduce() causes
!!! severe numerical problems in some cases.  I'm guessing there are memory
!!! overrun issues in this complex code, but have so far failed to identify
!!! a specific problem.  And/or, when one searches the Internet for this 
!!! issue, there is "some" hint that the implementation may be buggy.  
!!! 
!!! At this point, with the following CPP USE_MPIINPLACE directive, the
!!! default behaviour will be to not use the MPI_IN_PLACE option.
!!! Users will have to compile with -DUSE_MPIINPLACE if they want that option.
!!! Introduction of the CPP directives also requires that the code be compiled
!!! with the "-x f95-cpp-input" option.
!!!
!!! Modification of this section requires the addition of array gridunc0, which
!!! requires an allocation in outgrid_init.f90 and initial declaration in
!!! unc_mod.f90.
!!!--------------------------------------------------------------------

#ifdef USE_MPIINPLACE

    if (lroot) then
      call MPI_Reduce(MPI_IN_PLACE, gridunc, grid_size3d, mp_sp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    else
      call MPI_Reduce(gridunc, 0, grid_size3d, mp_sp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
    end if

#else

      call MPI_Reduce(gridunc, gridunc0, grid_size3d, mp_sp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (lroot) gridunc = gridunc0

#endif

    if ((WETDEP).and.(ldirect.gt.0)) then
      call MPI_Reduce(wetgridunc, wetgridunc0, grid_size2d, mp_cp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if

    if ((DRYDEP).and.(ldirect.gt.0)) then
      call MPI_Reduce(drygridunc, drygridunc0, grid_size2d, mp_cp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if

! Receptor concentrations    
    if (lroot) then
      call MPI_Reduce(MPI_IN_PLACE,creceptor,rcpt_size,mp_sp,MPI_SUM,id_root, &
           & mp_comm_used,mp_ierr)
      if (mp_ierr /= 0) goto 600
    else
      call MPI_Reduce(creceptor,0,rcpt_size,mp_sp,MPI_SUM,id_root, &
           & mp_comm_used,mp_ierr)
    end if

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop

601 end subroutine mpif_tm_reduce_grid


  subroutine mpif_tm_reduce_grid_nest
!***********************************************************************
! Collect nested grids to PID 0, adding from all processes.
!
! NOTE
!   - Compiling with 'fcheck=all' (gfortran) will cause a run-time error
!     as wetgriduncn0 is only allocated for root process. Possibly valid
!     MPI code but invalid gfortran.
!
!***********************************************************************
    use com_mod
    use unc_mod
    use par_mod

    implicit none

    integer :: grid_size2d,grid_size3d

!**********************************************************************

    grid_size3d=numxgridn*numygridn*numzgrid*maxspec* &
         & maxpointspec_act*nclassunc*maxageclass
    grid_size2d=numxgridn*numygridn*maxspec* &
         & maxpointspec_act*nclassunc*maxageclass


! Time for MPI communications
    if (mp_measure_time) call mpif_mtime('commtime',0)

! Using a separate grid (gridunc0) for received values, for debugging
! call MPI_Reduce(griduncn, griduncn0, grid_size3d, mp_sp, MPI_SUM, id_root, &
!      & mp_comm_used, mp_ierr)
! if (mp_ierr /= 0) goto 600

! Using in-place reduction
    if (lroot) then
      call MPI_Reduce(MPI_IN_PLACE, griduncn, grid_size3d, mp_sp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    else
      call MPI_Reduce(griduncn, 0, grid_size3d, mp_sp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
    end if

    if ((WETDEP).and.(ldirect.gt.0)) then
      call MPI_Reduce(wetgriduncn, wetgriduncn0, grid_size2d, mp_cp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if

    if ((DRYDEP).and.(ldirect.gt.0)) then
      call MPI_Reduce(drygriduncn, drygriduncn0, grid_size2d, mp_cp, MPI_SUM, id_root, &
           & mp_comm_used, mp_ierr)
      if (mp_ierr /= 0) goto 600
    end if

    if (mp_measure_time) call mpif_mtime('commtime',1)

    goto 601

600 write(*,*) "mpi_mod> mp_ierr \= 0", mp_ierr
    stop

601 end subroutine mpif_tm_reduce_grid_nest


  subroutine mpif_mtime(ident,imode)
!***********************************************************************
! Measure CPU/Wall time in various parts of the code
!
! VARIABLES
!   ident        character, identifies routine to measure
!   imode        integer, 0:start clock(s)  1: stop clock(s) 
!
!***********************************************************************
    implicit none

    character(LEN=*), intent(in) :: ident
    integer, intent(in) :: imode

!***********************************************************************

    select case(ident)

    case ('timemanager')
      if (imode.eq.0) then
        call cpu_time(tm_tot_beg)
      else
        call cpu_time(tm_tot_end)
        tm_tot_total = tm_tot_total + (tm_tot_end-tm_tot_beg)
      end if

    case ('wetdepo')
      if (imode.eq.0) then
        mp_wetdepo_wtime_beg = mpi_wtime()
        call cpu_time(mp_wetdepo_time_beg)
      else
        mp_wetdepo_wtime_end = mpi_wtime()
        call cpu_time(mp_wetdepo_time_end)

        mp_wetdepo_wtime_total = mp_wetdepo_wtime_total + &
             &(mp_wetdepo_wtime_end - mp_wetdepo_wtime_beg)
        mp_wetdepo_time_total = mp_wetdepo_time_total + &
             &(mp_wetdepo_time_end - mp_wetdepo_time_beg)
      end if

    case ('advance')
      if (imode.eq.0) then
        mp_advance_wtime_beg = mpi_wtime()
      else
        mp_advance_wtime_end = mpi_wtime()

        mp_advance_wtime_total = mp_advance_wtime_total + &
             &(mp_advance_wtime_end - mp_advance_wtime_beg)
      end if

    case ('getfields')
      if (imode.eq.0) then
        mp_getfields_wtime_beg = mpi_wtime()
        call cpu_time(mp_getfields_time_beg)
      else
        mp_getfields_wtime_end = mpi_wtime()
        call cpu_time(mp_getfields_time_end)

        mp_getfields_wtime_total = mp_getfields_wtime_total + &
             &(mp_getfields_wtime_end - mp_getfields_wtime_beg)
        mp_getfields_time_total = mp_getfields_time_total + &
             &(mp_getfields_time_end - mp_getfields_time_beg)
      end if

    case ('partloop1')
      if (imode.eq.0) then
        call cpu_time(tm_nploop_beg)
      else
        call cpu_time(tm_nploop_end)
        tm_nploop_total = tm_nploop_total + (tm_nploop_end - tm_nploop_beg)
      end if

    case ('conccalc')
      if (imode.eq.0) then
        call cpu_time(mp_conccalc_time_beg)
      else
        call cpu_time(mp_conccalc_time_end)
        mp_conccalc_time_total = mp_conccalc_time_total + mp_conccalc_time_end - &
             &mp_conccalc_time_beg
      end if

    case ('rootonly')
      if (imode.eq.0) then
        call cpu_time(mp_root_time_beg)
        mp_root_wtime_beg = mpi_wtime()
      else
        call cpu_time(mp_root_time_end)
        mp_root_wtime_end = mpi_wtime()
        mp_root_time_total = mp_root_time_total + &
             &(mp_root_time_end - mp_root_time_beg)
        mp_root_wtime_total = mp_root_wtime_total + &
             &(mp_root_wtime_end - mp_root_wtime_beg)
      end if

    case ('iotime')
      if (imode.eq.0) then
        mp_io_wtime_beg = mpi_wtime()
        call cpu_time(mp_io_time_beg)
      else
        mp_io_wtime_end = mpi_wtime()
        call cpu_time(mp_io_time_end)

        mp_io_wtime_total = mp_io_wtime_total + (mp_io_wtime_end - &
             & mp_io_wtime_beg)
        mp_io_time_total = mp_io_time_total + (mp_io_time_end - &
             & mp_io_time_beg)
      end if

    case ('verttransform')
      if (imode.eq.0) then
        mp_vt_wtime_beg = mpi_wtime()
        call cpu_time(mp_vt_time_beg)
      else
        mp_vt_wtime_end = mpi_wtime()
        call cpu_time(mp_vt_time_end)

        mp_vt_wtime_total = mp_vt_wtime_total + (mp_vt_wtime_end - &
             & mp_vt_wtime_beg)
        mp_vt_time_total = mp_vt_time_total + (mp_vt_time_end - &
             & mp_vt_time_beg)
      end if

    case ('readwind')
      if (imode.eq.0) then
        call cpu_time(mp_readwind_time_beg)
        mp_readwind_wtime_beg = mpi_wtime()
      else
        call cpu_time(mp_readwind_time_end)
        mp_readwind_wtime_end = mpi_wtime()

        mp_readwind_time_total = mp_readwind_time_total + &
             &(mp_readwind_time_end - mp_readwind_time_beg)
        mp_readwind_wtime_total = mp_readwind_wtime_total + &
             &(mp_readwind_wtime_end - mp_readwind_wtime_beg)
      end if

    case ('commtime')
      if (imode.eq.0) then
        call cpu_time(mp_comm_time_beg)
        mp_comm_wtime_beg = mpi_wtime()
      else
        call cpu_time(mp_comm_time_end)
        mp_comm_wtime_end = mpi_wtime()
        mp_comm_time_total = mp_comm_time_total + &
             &(mp_comm_time_end - mp_comm_time_beg)
        mp_comm_wtime_total = mp_comm_wtime_total + &
             &(mp_comm_wtime_end - mp_comm_wtime_beg)
      end if

    case ('flexpart')
      if (imode.eq.0) then
        mp_total_wtime_beg=mpi_wtime()
      else
        mp_total_wtime_end=mpi_wtime()
        mp_total_wtime_total = mp_total_wtime_end-mp_total_wtime_beg
      end if

    case default
      write(*,*) 'mpi_mod::mpif_mtime> WARNING: unknown case identifier'

    end select


  end subroutine mpif_mtime


  subroutine mpif_finalize
!***********************************************************************
! Finalize MPI
! Optionally print detailed time diagnostics
!***********************************************************************
    implicit none

    integer :: ip !,j,r

!***********************************************************************

    IF (mp_measure_time) THEN
      do ip=0, mp_np-1
        call MPI_BARRIER(MPI_COMM_WORLD, mp_ierr)

        if (ip == mp_pid) then
          write(*,FMT='(72("#"))')
          write(*,FMT='(A60,I3)') 'STATISTICS FOR MPI PROCESS:', mp_pid 
          if (ip == 0) then
            write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR ROOT-ONLY (PID 0) &
                 &SECTIONS: ', mp_root_time_total
            write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR ROOT-ONLY (PID 0) &
                 &SECTIONS:', mp_root_wtime_total
          end if
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR FLEXPART: ' &
               &, mp_total_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR MPI COMMUNICATIONS:',&
               & mp_comm_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR MPI COMMUNICATIONS:',&
               & mp_comm_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR MPI BARRIERS:',&
               & mp_barrier_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR MPI BARRIERS:',&
               & mp_barrier_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR TIMEMANAGER:',&
               & tm_tot_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR TIMEMANAGER/NUMPART LOOP:',&
               & tm_nploop_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR ADVANCE:',&
               & mp_advance_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR GETFIELDS:',&
               & mp_getfields_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR GETFIELDS:',&
               & mp_getfields_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR READWIND:',&
               & mp_readwind_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR READWIND:',&
               & mp_readwind_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR FILE IO:',&
               & mp_io_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR FILE IO:',&
               & mp_io_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR WETDEPO:',&
               & mp_wetdepo_wtime_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR WETDEPO:',&
               & mp_wetdepo_time_total
          write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR CONCCALC:',&
               & mp_conccalc_time_total
! write(*,FMT='(A60,TR1,F9.2)') 'TOTAL WALL TIME FOR VERTTRANSFORM:',&
!      & mp_vt_wtime_total
! write(*,FMT='(A60,TR1,F9.2)') 'TOTAL CPU TIME FOR VERTTRANSFORM:',&
!      & mp_vt_time_total
! NB: the 'flush' function is possibly a gfortran-specific extension,
! comment out if it gives problems
!          call flush()
        end if
      end do
    end if

! This call to barrier is for correctly formatting output
    call MPI_BARRIER(MPI_COMM_WORLD, mp_ierr)

    if (lroot.and.mp_measure_time) then
      write(*,FMT='(72("#"))')
      WRITE(*,*) "To turn off output of time measurements, set "
      WRITE(*,*) "    mp_measure_time=.false."
      WRITE(*,*) "in file mpi_mod.f90"
      write(*,FMT='(72("#"))')
    end if

! j=mp_pid*nvar_async
! In the implementation with 3 fields, the processes may have posted
! MPI_Irecv requests that should be cancelled here
! if (.not.lmp_sync) then
!   r=mp_pid*nvar_async
!   do j=r,r+nvar_async-1
!     call MPI_Cancel(j,mp_ierr)
!     if (mp_ierr /= 0) write(*,*) '#### mpif_finalize::MPI_Cancel> ERROR ####'
!   end do
! end if

    call MPI_FINALIZE(mp_ierr)
    if (mp_ierr /= 0) then
      write(*,*) '#### mpif_finalize::MPI_FINALIZE> MPI ERROR ', mp_ierr, ' ####'
      stop
    end if


  end subroutine mpif_finalize


  subroutine get_lun(my_lun)
!***********************************************************************
! get_lun: 
!   Starting from 100, get next free logical unit number
!***********************************************************************

    implicit none

    integer, intent(inout) :: my_lun
    integer, save :: free_lun=100
    logical :: exists, iopen

!***********************************************************************

    loop1: do
      inquire(UNIT=free_lun, EXIST=exists, OPENED=iopen)
      if (exists .and. .not.iopen) exit loop1
      free_lun = free_lun+1
    end do loop1
    my_lun = free_lun

  end subroutine get_lun


  subroutine write_data_dbg(array_in, array_name, tstep, ident)
!***********************************************************************
! Write one-dimensional arrays to file (for debugging purposes)
!***********************************************************************
    implicit none 

    real, intent(in), dimension(:) :: array_in
    integer, intent(in) :: tstep
    integer :: lios
    character(LEN=*), intent(in) :: ident, array_name

    character(LEN=8) :: c_ts
    character(LEN=40) :: fn_1, fn_2

!***********************************************************************

    write(c_ts, FMT='(I8.8,BZ)') tstep
    fn_1='-'//trim(adjustl(c_ts))//'-'//trim(ident)
    write(c_ts, FMT='(I2.2,BZ)') mp_np
    fn_2= trim(adjustl(array_name))//trim(adjustl(fn_1))//'-np'//trim(adjustl(c_ts))//'.dat'

    call get_lun(dat_lun)
    open(UNIT=dat_lun, FILE=fn_2, IOSTAT=lios, ACTION='WRITE', &
         FORM='UNFORMATTED', STATUS='REPLACE')
    write(UNIT=dat_lun, IOSTAT=lios) array_in
    close(UNIT=dat_lun)

  end subroutine write_data_dbg


  subroutine set_fields_synthetic()
!*******************************************************************************
! DESCRIPTION
!   Set all meteorological fields to synthetic (usually constant/homogeneous)
!   values.
!   Used for validation and error-checking
!
! NOTE
!   This version uses asynchronious communications.
!
! VARIABLES
!   
!
!
!*******************************************************************************
    use com_mod

    implicit none

    integer :: li=1, ui=2 ! wfmem indices (i.e, operate on entire field)

    if (.not.lmp_sync) ui=3


! Variables transferred at initialization only
!*********************************************
!    readclouds=readclouds_
    oro(:,:)=0.0
    excessoro(:,:)=0.0
    lsm(:,:)=0.0
    xlanduse(:,:,:)=0.0
!    wftime
!    numbwf
!    nmixz
!    height

! Time-varying fields:
    uu(:,:,:,li:ui) = 10.0
    vv(:,:,:,li:ui) = 0.0
    uupol(:,:,:,li:ui) = 10.0
    vvpol(:,:,:,li:ui)=0.0
    ww(:,:,:,li:ui)=0.
    tt(:,:,:,li:ui)=300.
    rho(:,:,:,li:ui)=1.3
    drhodz(:,:,:,li:ui)=.0
    tth(:,:,:,li:ui)=0.0
    qvh(:,:,:,li:ui)=1.0
    qv(:,:,:,li:ui)=1.0

    pv(:,:,:,li:ui)=1.0
    clouds(:,:,:,li:ui)=0

    clwc(:,:,:,li:ui)=0.0
    ciwc(:,:,:,li:ui)=0.0

! 2D fields

    cloudsh(:,:,li:ui)=0
    vdep(:,:,:,li:ui)=0.0
    ps(:,:,:,li:ui)=1.0e5
    sd(:,:,:,li:ui)=0.0
    tcc(:,:,:,li:ui)=0.0
    tt2(:,:,:,li:ui)=300.
    td2(:,:,:,li:ui)=250.
    lsprec(:,:,:,li:ui)=0.0
    convprec(:,:,:,li:ui)=0.0
    ustar(:,:,:,li:ui)=1.0
    wstar(:,:,:,li:ui)=1.0
    hmix(:,:,:,li:ui)=10000.
    tropopause(:,:,:,li:ui)=10000.
    oli(:,:,:,li:ui)=0.01

  end subroutine set_fields_synthetic

end module mpi_mod
