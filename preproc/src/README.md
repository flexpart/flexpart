# README #

This documentation shows how to compile `CONVERT2` program

### Overview ###

### Requirements ###

| Utilities             |
|----------------------:| 
|  make                 | [http://www.gnu.org/software/make/](http://www.gnu.org/software/make/) |
| Third party libraries |
|----------------------:| 
|  EMOS                 | [https://software.ecmwf.int/wiki/display/EMOS/Emoslib](https://software.ecmwf.int/wiki/display/EMOS/Emoslib) | Make sure you install EMOS lib with 64 bits reals.
|  grib-api             | [https://software.ecmwf.int/wiki/display/GRIB/Home](https://software.ecmwf.int/wiki/display/GRIB/Home) | Make sure you install GRIB-API with JPEG support and python GRIB-API.


### Installation ###

*  Environment

At UIO, Red Hat 6 Linux systems (64 bits) were used for testing. We use the [Module package](http://modules.sourceforge.net/) to set-up user environment.

* Getting the source code

If you still haven't got the source from bitbucket (Please note that it needs to be done once only).

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

There is a makefile example called `makefile` which you can use as a template. It uses GNU gfortran version 4.8.2. The version of your `gfortran`compiler can be tested with the following command:
`
gfortran --version
`

For using our `makefile`, you need to define the following environment variables:

Please note that these environment variables may be automatically defined on your system. Please check with your local IT administrator. For instance at UIO, these environment variables are defined by loading the corresponding modulefiles (`module load grib_api emos`).

+ `GRIB_API_FFLAGS` 
+ `GRIB_API_LDFLAGS`
+ `EMOS_LDFLAGS`

For instance at UIO on sverdrup.uio.no:
`
> echo $GRIB_API_FFLAGS
-I/site/opt/grib_api/1.12.3_gnu/include

> echo $GRIB_API_LDFLAGS
-L/site/opt/grib_api/1.12.3_gnu/lib -lgrib_api_f90 -lgrib_api

> echo $EMOS_LDFLAGS
-L/site/opt/emos/000392_grib_api_1.12.3_gnu -lemos

`

If `EMOS` and `GRIB-API` have been installed properly and you have adjusted your makefile then you can compile `CONVERT2`:
`
cd preproc/src

gmake
`

If successful, a program called `CONVERT2` is created in the current directory. You now need to add this program in your `PATH`:
** Korn-shell or Bash users:
`
export PATH=$FLEXPART_HOME/preproc/src:$PATH
`

** C-shell users:

`
setenv PATH $FLEXPART_HOME/preproc/src:$PATH
`

Where `$FLEXPART_HOME` is the directory where FLEXPART

* Testing your installation

See separate instructions in `tests/preproc`.

Please report any problems.

###  Installation FAQ ###

