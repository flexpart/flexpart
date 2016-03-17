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
      

!*******************************************************************************
!   Include file for calculation of particle trajectories (Program FLEXPART)   *
!        This file contains the parameter statements used in FLEXPART          *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        1997                                                                  *
!                                                                              *
!        Last update 10 August 2000                                            *
!                                                                              * 
!        modified by m.Cassiani, 2016 to use NorESM meteo files                *
!                                                                              *
!*******************************************************************************

module par_mod

  implicit none
  
  integer, parameter :: yearorigin=-4800 !look in juldate for this: comment by mc

  !****************************************************************
  ! Parameter defining KIND parameter for "double precision"
  !****************************************************************

  integer,parameter :: dp=selected_real_kind(P=15)


  !***********************************************************
  ! Number of directories/files used for FLEXPART input/output
  !***********************************************************

  integer,parameter :: numpath=5 !one additional path for NorESM TO ACCOUNT FOR THE GRID FILE IN CAM4: comment by mc

  ! numpath                 Number of different pathnames for input/output files

  integer, parameter :: maxdim=10,maxvar=70,maxtime=12 ! parameter used in readwind and gridcheck for reading NORESM output file: commmnt added by mc
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
  integer,parameter :: nxmax=145,nymax=97,nuvzmax=27,nwzmax=27,nzmax=27  !for NorESM with 1.875x2.5 grid: added by mc
  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92,nzmax=92
  !integer,parameter :: nxmax=361,nymax=181,nuvzmax=61,nwzmax=61,nzmax=61
  !integer,parameter :: nxmax=721,nymax=361,nuvzmax=64,nwzmax=64,nzmax=64
  !integer,parameter :: nxshift=359  ! for ECMWF
  !nteger,parameter :: nxshift=0     ! for GFS
  integer,parameter :: nxshift=73     ! added by mc for CAM 4.0 NOTE THAT THE DATA OF CAM 4.0 HAVE LON(1)=0 with 144 nodes of 2.5 degree, therefore we move origin to -177.5.
  
  integer,parameter :: nconvlevmax = nuvzmax-1
  integer,parameter :: na = nconvlevmax+1


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

  integer,parameter :: maxnests=0, nxmaxn=0, nymaxn=0 ! <<<-----------------these must be always zero for NorESM: comment by mc
 

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
  integer,parameter :: maxageclass=9,nclassunc=1

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

  integer,parameter :: maxpart=20000000
  integer,parameter :: maxspec=4


  ! maxpart                 Maximum number of particles
  ! maxspec                 Maximum number of chemical species per release

  ! maxpoint is also set dynamically during runtime
  ! maxpoint                Maximum number of release locations

  ! ---------
  ! Sabine Eckhardt: change of landuse inventary numclass=13
  ! ---------
  integer,parameter :: maxwf=50000, maxtable=1000, numclass=13, ni=11

  ! maxwf                   maximum number of wind fields to be used for simulation
  ! maxtable                Maximum number of chemical species that can be
  !                         tabulated for FLEXPART
  ! numclass                Number of landuse classes available to FLEXPART
  ! ni                      Number of diameter classes of particles

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

  integer,parameter :: maxrand=2000000

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
  integer,parameter :: unitdates=94, unitheader=90, unitshortpart=95
  integer,parameter :: unitboundcond=89
  integer,parameter :: unitdiagnostic1=72,unitdiagnostic2=73,unitdiagnostic3=74 !number from 70 to 79 reserved for diagnostic

end module par_mod
