# README #

This documentation shows how to use these python scripts to extract ECMWF ERA-Interim data and generate WINDFIELDS for running FLEXPART. 

### Overview ###

To run FLEXPART with ECMWF ERA-Interim data, you first need to retrieve ECMWF ERA-Interim GRIB fields and generate FLEXPART WINDFIELDS to run FLEXPART.

The two main programs are called respectively `getEIdata.py` and `prepareFLEXPART.py`. 

To get the usage of these programs, use `-h` option:

`getEIdata.py -h`

Optional arguments are mentionned in squared brackets. 

* `getEIdata.py`

This program allows you to download ECMWF ERA-Interim data from ECMWF using [ecmwfapi](https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home). 
It requires `ecmwfapi` python library (see Requirements below). Check with your local IT group as it may be already available.


{{{
Usage: getEIdata.py --start_date=YYYYMMDD [--end_date=YYYYMMDD] [--times=tt1/tt2/tt3] [--levels=nlevels]
                                                  [--area=north/west/south/east]  [--outputdir=output_directory]

Options:
  -h, --help            show this help message and exit
  --start_date=start_date
                        start date YYYYMMDD
  --end_date=end_date   end_date YYYYMMDD
  --times=times         times such as 00/12
  --levels=levels       number of vertical levels
  --area=area           area defined as north/west/south/east with default
                        90.0/-179.0/-90.0/180.0
  --outputdir=outputdir
                        root directory for storing output files
}}}

* `prepareFLEXPART.py`

This program allow you to generate FLEXPART WINDFIELDS (inputs for FLEXPART). It requires python interface to grib_api and `CONVERT2` program (located in `src` directory with instruction on how to compile it). You also need to provide a namelist for CONVERT2 (see test_1).

`
Usage: prepareFLEXPART.py --start_date=YYYYMMDD [--end_date=YYYYMMDD] [--namelist=namelist_for_convert] [--inputdir=input_root_directory] [--outputdir=output_directory]

Options:
  -h, --help            show this help message and exit
  --start_date=start_date
                        start date YYYYMMDD
  --end_date=end_date   end_date YYYYMMDD
  --namelist=namelist   namelist used for converting
  --inputdir=inputdir   root directory for reading input files
  --outputdir=outputdir
                        root directory for storing output files

`
### Requirements ###

| Python Support        |
|----------------------:| 
| python                | [http://www.python.org](http://www.python.org)  | We have use [python Anaconda](https://store.continuum.io/cshop/anaconda/) for our testings
| python-numpy          | [http://www.numpy.org/](http://www.numpy.org/)  | Not necessary if you have installed python Anaconda
| ecmwfapi              | [https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home](https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home) | You also need to install your API key (as explained in the documentation)
| Utilities             |
|----------------------:| 
|  grib-api             | [https://software.ecmwf.int/wiki/display/GRIB/Home](https://software.ecmwf.int/wiki/display/GRIB/Home) | Make sure you install GRIB-API with JPEG support and python GRIB-API.
|----------------------:| 
| FLEXPART programs     | to run these programs (prepareFLEXPART.py), you need to compile CONVERT2 program (located in preproc/src). See separate README file to get instructions on how to compile this code.
|----------------------:| 


### Installation ###

*  Environment

At UIO, Red Hat 6 Linux systems (64 bits) were used for testing. We use the [Module package](http://modules.sourceforge.net/) to set-up user environment.

* Getting the source code
In the directory of your choice:

`
git clone git@bitbucket.org:flexpart/flexpart.git
`
This command will create a subdirectory called flexpart: it contains the latest FLEXPART version.

Set then environment variable `FLEXPART_HOME`:

** Korn-shell or Bash users:
`
cd flexpart
export FLEXPART_HOME=$PWD
`

** C-shell users:

`
cd flexpart
setenv FLEXPART_HOME=$PWD
`

* Installation

Make sure your first generate `CONVERT2` program (see separate instructions in `preproc/src`).

Users need to be able to execute prepareFLEXPART.py and getEIdata.py so make sure they have the correct unix permissions:

`
cd preproc/python
chmod uog+rx getEIdata.py prepareFLEXPART.py
`

These two programs must be in the user PATH. At UIO this is done automatically when loading flexpart. If not, you would need to do the following:

** Korn-shell or Bash users:
`
export PATH=$FLEXPART_HOME/preproc/python:$PATH
`

** C-shell users:

`
setenv PATH $FLEXPART_HOME/preproc/python:$PATH
`

Where `$FLEXPART_HOME` is the directory where FLEXPART 

* Testing your installation

First check that grib-api python interface is correctly installed on your platform:
`
python
>>> from gribapi import *
>>>
`
Use `CTRL-D` to quit python.

Then check that `ecmwfapi` is properly installed:
`
python
>>> from ecmwfapi import *
>>>
`

If the two previous tests were successful, you can run `tests/preproc` (See separate instructions in `tests/preproc`). 

If any of these two tests fail, this probably means that either `ecmwfapi` or `grib-api` have not been installed properly. 

Please report any problems.

###  Installation FAQ ###

