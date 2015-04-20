
		FLEXPART VERSION 10.0 beta (MPI)

Description
-----------

  This branch contains both the standard (serial) FLEXPART, and a parallel 
  version (implemented with MPI). The latter is under developement, so not 
  every FLEXPART option is implemented yet.

  MPI related subroutines and variables are in file mpi_mod.f90.

  Most of the source files are identical/shared between the serial and 
  parallel versions. Those that depend on the MPI module have '_mpi' 
  apppended to their names, e.g. 'timemanager_mpi.f90'


Installation
------------

  A MPI library must be installed on the target platform, either as a 
  system library or compiled from source.

  So far, we have tested the following freely available implementations: 	   
  mpich2  -- versions 3.0.1, 3.0.4, 3.1, 3.1.3
  OpenMPI -- version 1.8.3

  Based on testing so far, OpenMPI is recommended.

  Compiling the parallel version (executable: FP_ecmwf_MPI) is done by

    'make [-j] ecmwf-mpi'

  The makefile has resolved dependencies, so 'make -j' will compile 
  and link in parallel. 

  The included makefile must be edited to match the target platform 
  (location of system libraries, compiler etc.).


Usage
-----

  Running the parallel version with MPI is done with the "mpirun" command
  (some MPI implementations may use a "mpiexec" command instead). The 
  simplest case is:

    'mpirun -n [number] ./FP_ecmwf_MPI'

  where 'number' is the number of processes to launch. Depending on the
  target platform, useful options regarding process-to-processor bindings
  can be specified (for performance reasons), e.g,

    'mpirun --bind-to l3cache -n [number] ./FP_ecmwf_MPI'


Implementation
--------------

  The current parallel model is based on distributing particles equally
  among the running processes. In the code, variables like 'maxpart' and 
  'numpart' are complemented by variables 'maxpart_mpi' and 'numpart_mpi'
  which are the run-time determined number of particles per process, i.e,
  maxpart_mpi = maxpart/np, where np are the number of processes. The variable 'numpart' 
  is still used in the code, but redefined to mean 'number of particles
  per MPI process'

  The root MPI process writes concentrations to file, following a MPI
  communication step where each process sends its contributions to root, 
  where the individual contributions are summed.

  In the parallel version one can choose to set aside a process dedicated
  to reading and distributing meteorological data ("windfields"). This process will
  thus not participate in the calculation of trajectories. This might not be
  the optimal choice when running with very few processes.
  As an example, running with a total number of processes np=4 and
  using one of these processes for reading windfields will normally
  be faster than running with np=3 and no dedicated 'reader' process. 
  But it is also possible that the
  program will run even faster if the 4th process is participating in 
  the calculation of particle trajectories instead. This will largely depend on
  the problem size (total number of particles in the simulation, resolution
  of grids etc) and hardware being used (disk speed/buffering, memory
  bandwidth etc).

  To control this
  behavior, edit the parameter 'read_grp_min' in file mpi_mod.f90. This 
  sets the minimum number of total processes at which one will be set 
  aside for reading the fields. Experimentation is required to find 
  the optimum value. At typical NILU machines (austre.nilu.no, 
  dmz-proc01.nilu.no) with 24-32 cores, a value of 6-8 seems to be a 
  good choice.

  An experimental feature, which is an extension of the functionality
  described above, is to hold 3 fields in memory instead of the usual 2.
  Here, the transfer of fields from the "reader" process to the "particle"
  processes is done on the vacant field index, simultaneously while the
  "particle" processes are calculating trajectories. To use this feature,
  set 'lmp_sync=.false'. in file mpi_mod.f90 and set numwfmem=3 in file
  par_mod.f90. At the moment, this method does not seem to produce faster
  running code (about the same as the "2-fields" version).
  

Performance efficency considerations
------------------------------------

  A couple of reference runs have been set up to measure performace of the
  MPI version (as well as checking for errors in the implementation).
  They are as follows:
   
  Reference run 1 (REF1):
    * Forward modelling (24h) of I2-131, variable number of particles
    * Two release locations 
    * 360x720 Global grid, no nested grid
    * Species file modified to include (not realistic) values for
        scavenging/deposition
 

  As the parallization is based on particles, it follows that if  
  FLEXPART-MPI is run with no (or just a few) particles, no performance 
  improvement is possible. In this case, most processing time is spent
  in the 'getfields'-routine (ECMWF).

  A) Running without dedicated reader process
  ----------------------------------------
  Running REF1 with 100M particles on 16 processes (NILU machine 'dmz-proc04'), 
  a speedup close to 8 is observed (~50% efficiency).

  Running REF1 with 10M particles on 8 processes (NILU machine 'dmz-proc04'), 
  a speedup close to 3 is observed (~40% efficiency). Running with 16
  processes gives only marginal improvements (speedup ~3.5) because of the 'getfields'
  bottleneck.
  
  Running REF1 with 1M particles: Here 'getfields' consumes ~70% of the CPU
  time. Running with 4 processes gives a speedup of ~1.5. Running with more
  processes does not help much here.

  B) Running with dedicated reader process
  ----------------------------------------

  Running REF1 with 40M particles on 16 processes (NILU machine 'dmz-proc04'), 
  a speedup above 10 is observed (~63% efficiency).

  :TODO: more to come...


Advice  
------
  From the tests referred to above, the following advice can be given:

    * Do not run with too many processes.
    * Do not use the parallel version when running with very few particles.
      

What is implemented in the MPI version
--------------------------------------

 -The following should work (have been through initial testing): 

    * Forward runs
    * OH fields
    * Radioactive decay
    * Particle splitting
    * Dumping particle positions to file
    * ECMWF data
    * Wet/dry deposition
    * Nested grid output
    * NetCDF output
    * Namelist input/output

 -Implemented but untested:
    * Domain-filling trajectory calculations
    * Nested wind fields

 -The following will most probably not work (untested/under developement): 

    * Backward runs

 -This will positively NOT work yet

    * Subroutine partoutput_short (MQUASILAG = 1) will not dump particles
      correctly at the moment
    * Reading particle positions from file (the tools to implement this
      are available in mpi_mod.f90 so it will be possible soon)


  Please keep in mind that running the serial version (FP_ecmwf_gfortran)
  should yield identical results as running the parallel version
  (FP_ecmwf_MPI) using only one process, i.e. "mpirun -n 1 FP_ecmwf_MPI".
  If not, this indicates a bug.
  
  When running with multiple processes, statistical differences are expected
  in the results.

Contact
-------

  If you have questions, or wish to work with the parallel version, please 
  contact Espen Sollum (eso@nilu.no). Please report any errors/anomalies!
