# README #

This documentation shows how to extract ECMWF ERA-Interim GRIB fields and generate FLEXPART WINDFIELDS. 


Make sure you have installed FLEXPART preproc properly (see corresponding installation instructions).

./run_retrieve

--> It will generate a directory called 2007 containing ECMWF ERA-Interim GRIB fields

--> And a directory called WIND_FIELDS containing FLEXPART WINDFIELDS to run FLEXPART.

After this test, you can remove ECMWF ERA-Interim GRIB fields as they are no longer needed:

rm -rf 2007

Then you can test your new generated WIND_FIELDS (see separate example).

