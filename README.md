# README #


### What is this repository for? ###

* This repository contains the development version of the Lagrangian model FLEXPART 

 ### How do I get set up? ###

* Configuration

  Edit the makefile with the paths to libraries and include files

* Dependencies

 * Jasper and grib_api or ECCodes
 * NetCDF (optional)

* Compilation

```
> cd src
> make 
```

* Deployment instructions 

   FLEXPART is a standalone executable   

### Contribution guidelines ###

* The version contributed should compile on a reference version of the system and compiler. The current reference is gfortran 5.4 on Ubuntu 16.04
* Code contribution including new features and bug fixes should be complemented with appropriate tests
   An essential test consists of a set of input files and directories that allow FLEXPART to run.
   A test can be accompanied by output files for verification
* Code review

[comment]: # "### Who do I talk to? ###"

ignacio.pisso@nilu.no
