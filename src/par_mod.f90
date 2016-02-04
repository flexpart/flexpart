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

!*******************************************************************************
!   Include file for calculation of particle trajectories (Program FLEXPART)   *
!        This file contains the parameter statements used in FLEXPART          *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        1997                                                                  *
!                                                                              *
!        Last update 15 August 2013 IP                                         *
!                                                                              *
!*******************************************************************************

module par_mod

! wind_mod: is gfs_mod.f90 for target gfs, ecmwf_mod.f90 for target ecmwf
!************************************************************************
  use wind_mod 

  implicit none

  !****************************************************************
  ! Parameters defining KIND parameter for double/single precision
  !****************************************************************

  integer,parameter :: dp=selected_real_kind(P=15)
  integer,parameter :: sp=selected_real_kind(6)

  !****************************************************************
  ! dep_prec sets the precision for deposition calculations (sp or 
  ! dp). sp is default, dp can be used for increased precision.
  !****************************************************************

  integer,parameter :: dep_prec=sp

  !***********************************************************
  ! Number of directories/files used for FLEXPART input/output
  !***********************************************************

  integer,parameter :: numpath=4

  ! numpath                 Number of different pathnames for input/output files


  !*****************************
  ! Physical and other constants
  !*****************************

  real,parameter :: pi=3.14159265, r_earth=6.371e6, r_air=287.05, ga=9.81
  real,parameter :: cpa=1004.6, kappa=0.286, pi180=pi/180., vonkarman=0.4

  ! pi                      number "pi"
  ! pi180                   pi/180.
  ! r_earth                 radius of earth [m]
  ! r_air                   individual gas constant for dry air [J/kg/K]
  ! ga                      gravity acceleration of earth [m/s**2]
  ! cpa                     specific heat for dry air
  ! kappa                   exponent of formula for potential temperature
  ! vonkarman               von Karman constant

  real,parameter :: karman=0.40, href=15., convke=2.0
  real,parameter :: hmixmin=100., hmixmax=4500., turbmesoscale=0.16
  real,parameter :: d_trop=50., d_strat=0.1
  real,parameter :: rho_water=1000. !ZHG 2015 [kg/m3]
  ! karman                  Karman's constant
  ! href [m]                Reference height for dry deposition
  ! konvke                  Relative share of kinetic energy used for parcel lifting
  ! hmixmin,hmixmax         Minimum and maximum allowed PBL height
  ! turbmesoscale           the factor by which standard deviations of winds at grid
  !                    points surrounding the particle positions are scaled to
  !                    yield the scales for the mesoscale wind velocity fluctuations
  ! d_trop [m2/s]           Turbulent diffusivity for horizontal components in the troposphere
  ! d_strat [m2/s]          Turbulent diffusivity for vertical component in the stratosphere

  real,parameter :: xmwml=18.016/28.960

  ! xmwml   ratio of molar weights of water vapor and dry air
  !****************************************************
  ! Constants related to the stratospheric ozone tracer
  !****************************************************

  real,parameter :: ozonescale=60., pvcrit=2.0

  ! ozonescale              ppbv O3 per PV unit
  ! pvcrit                  PV level of the tropopause



  !********************
  ! Some time constants
  !********************

  integer,parameter :: idiffnorm=10800, idiffmax=2*idiffnorm, minstep=1

  ! idiffnorm [s]           normal time interval between two wind fields
  ! idiffmax [s]            maximum time interval between two wind fields
  ! minstep [s]             minimum time step to be used within FLEXPART


  !*****************************************************************
  ! Parameters for polar stereographic projection close to the poles
  !*****************************************************************

  real,parameter :: switchnorth=75., switchsouth=-75.

  ! switchnorth    use polar stereographic grid north of switchnorth
  ! switchsouth    use polar stereographic grid south of switchsouth


  !*********************************************
  ! Maximum dimensions of the input mother grids
  !*********************************************
  
  ! nxmax,nymax,nuvzmax,nwzmax,nzmax:
  ! Moved to ecmwf.f90 (for ECMWF) / gfs.f90 (GFS)
  

  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92,nzmax=92 !FNL XF
  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=152,nwzmax=152,nzmax=152 !ECMWF new 
  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92,nzmax=92 !ECMWF
  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=26,nwzmax=26,nzmax=26
  !integer,parameter :: nxmax=721,nymax=361,nuvzmax=64,nwzmax=64,nzmax=64
  !integer,parameter :: nxmax=1201,nymax=235,nuvzmax=58,nwzmax=58,nzmax=58

!  integer,parameter :: nxshift=359 ! for ECMWF
  !integer,parameter :: nxshift=0     ! for GFS or FNL

  integer,parameter :: nconvlevmax = nuvzmax-1
  integer,parameter :: na = nconvlevmax+1

  ! moved to gfs_mod.f90 / ecmwf_mod.f90
  !
  ! nxmax,nymax        maximum dimension of wind fields in x and y
  !                    direction, respectively
  ! nuvzmax,nwzmax     maximum dimension of (u,v) and (w) wind fields in z
  !                    direction (for fields on eta levels)
  ! nzmax              maximum dimension of wind fields in z direction
  !                    for the transformed Cartesian coordinates
  ! nxshift            for global grids (in x), the grid can be shifted by
  !                    nxshift grid points, in order to accomodate nested
  !                    grids, and output grids overlapping the domain "boundary"
  !                    nxshift must not be negative; "normal" setting would be 0

  ! ntracermax         maximum number of tracer species in convection
  ! nconvlevmax        maximum number of levels for convection
  ! na                 parameter used in Emanuel's convect subroutine


  !*********************************************
  ! Maximum dimensions of the nested input grids
  !*********************************************

  !integer,parameter :: maxnests=0, nxmaxn=0, nymaxn=0
  !integer,parameter :: maxnests=0,nxmaxn=351,nymaxn=351 !ECMWF

  !integer,parameter :: maxnests=1, nxmaxn=201, nymaxn=161 ! FNL XF
  ! maxnests                maximum number of nested grids
  ! nxmaxn,nymaxn           maximum dimension of nested wind fields in
  !                         x and y direction, respectively


  !*********************************
  ! Parmaters for GRIB file decoding
  !*********************************

  integer,parameter :: jpack=4*nxmax*nymax, jpunp=4*jpack

  ! jpack,jpunp             maximum dimensions needed for GRIB file decoding


  !**************************************
  ! Maximum dimensions of the output grid
  !**************************************

  !integer,parameter :: maxageclass=1,maxzgrid=10,nclassunc=1
  integer,parameter :: maxageclass=1,nclassunc=1

  ! nclassunc               number of classes used to calculate the uncertainty
  !                         of the output
  ! maxageclass             maximum number of age classes used for output

  ! Sabine Eckhardt, June, 2008
  ! the dimensions of the OUTGRID are now set dynamically during runtime
  ! maxxgrid,maxygrid,maxzgrid    maximum dimensions in x,y,z direction
  ! maxxgridn,maxygridn           maximum dimension of the nested grid
  !integer maxxgrid,maxygrid,maxzgrid,maxxgridn,maxygridn
  !integer,parameter :: maxxgrid=361,maxygrid=181,maxxgridn=0,maxygridn=0)

  integer,parameter :: maxreceptor=200

  ! maxreceptor             maximum number of receptor points


  !**************************************************
  ! Maximum number of particles, species, and similar
  !**************************************************

  integer,parameter :: maxpart=40000000
  integer,parameter :: maxspec=6
  integer,parameter :: minmass=0.0 !0.0001

  ! maxpart                 Maximum number of particles
  ! maxspec                 Maximum number of chemical species per release
  ! minmass                 Terminate particles carrying less mass

  ! maxpoint is also set dynamically during runtime
  ! maxpoint                Maximum number of release locations

  ! ---------
  ! Sabine Eckhardt: change of landuse inventary numclass=13
  ! ---------
  integer,parameter :: maxwf=50000, maxtable=1000, numclass=13, ni=11
  integer,parameter :: numwfmem=2 ! Serial version/MPI with 2 fields
  !integer,parameter :: numwfmem=3 ! MPI with 3 fields

  ! maxwf                   maximum number of wind fields to be used for simulation
  ! maxtable                Maximum number of chemical species that can be
  !                         tabulated for FLEXPART
  ! numclass                Number of landuse classes available to FLEXPART
  ! ni                      Number of diameter classes of particles
  ! numwfmem                Number of windfields kept in memory. 2 for serial
  !                         version, 2 or 3 for MPI version

  !**************************************************************************
  ! dimension of the OH field
  !**************************************************************************
  integer,parameter :: maxxOH=72, maxyOH=46, maxzOH=7

  !**************************************************************************
  ! Maximum number of particles to be released in a single atmospheric column
  ! for the domain-filling trajectories option
  !**************************************************************************

  integer,parameter :: maxcolumn=3000


  !*********************************
  ! Dimension of random number field
  !*********************************

!  integer,parameter :: maxrand=120000000
  integer,parameter :: maxrand=200000000
!
  ! maxrand                 number of random numbers used
  

  !*****************************************************
  ! Number of clusters to be used for plume trajectories
  !*****************************************************

  integer,parameter :: ncluster=5

  !************************************
  ! Unit numbers for input/output files
  !************************************

  integer,parameter :: unitpath=1, unitcommand=1, unitageclasses=1, unitgrid=1
  integer,parameter :: unitavailab=1, unitreleases=88, unitpartout=93
  integer,parameter :: unitpartin=93, unitflux=98, unitouttraj=96
  integer,parameter :: unitvert=1, unitoro=1, unitpoin=1, unitreceptor=1
  integer,parameter :: unitoutgrid=97, unitoutgridppt=99, unitoutinfo=1
  integer,parameter :: unitspecies=1, unitoutrecept=91, unitoutreceptppt=92
  integer,parameter :: unitlsm=1, unitsurfdata=1, unitland=1, unitwesely=1
  integer,parameter :: unitOH=1
  integer,parameter :: unitdates=94, unitheader=90,unitheader_txt=100, unitshortpart=95
  integer,parameter :: unitboundcond=89
  integer,parameter :: unittmp=101

!******************************************************
! integer code for missing values, used in wet scavenging (PS, 2012)
!******************************************************

  integer,parameter ::  icmv=-9999

! Parameters for testing
!*******************************************
!  integer :: verbosity=0

end module par_mod
