#!/usr/bin/env python
#
# (C) Copyright 2014 UIO.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# 
# Creation: July 2014 - Anne Fouilloux - University of Oslo
#
# Get ERA-Interim GRIB fields from ECMWF to prepare WINDFIELDS for running
# FLEXPART examples.
# 
# To run this example, you need to self-register at ECMWF data server and
# get an API key available from https://api.ecmwf.int/v1/key/
#
# Inputs parameters:
#  start_date: starting date for retrieving ECMWF data
#  end_date  : end date for retrieving ECMWF data; set to start_date if not specified
#  times     : forecast or analysis time for MARS retrievals
#  levels    : number of levels; set to 60 by default (for ERA-interim)
#  area      : you can define your area with --area=north/west/south/east (values are in degrees)
#  outputdir : location where you want to store ERA-Interim outputs
# 
from ecmwfapi import ECMWFDataServer
import calendar
import shutil
import datetime
import time
import os
from string import strip
from optparse import OptionParser
from FlexpartTools import MARSretrieval, EIFlexpart, daterange, mkdir_p, silentremove, months_between


def main():
    usage = """usage: %prog --start_date=YYYYMMDD [--end_date=YYYYMMDD] [--times=tt1/tt2/tt3] [--levels=nlevels]
                                                  [--area=north/west/south/east]  [--outputdir=output_directory] """
    parser = OptionParser(usage=usage)
    parser.add_option("--start_date", dest="start_date",
                      help="start date YYYYMMDD", metavar="start_date" )
    parser.add_option( "--end_date", dest="end_date",
                      help="end_date YYYYMMDD", metavar="end_date")
    parser.add_option("--times", dest="times",
                      default="00/03/06/09/12/15/18/21", help="times such as 00/12", metavar="times")
    parser.add_option("--levels", dest="levels",
                      default="60",help="number of vertical levels", metavar="levels")
    parser.add_option("--area", dest="area",
                      default="90.0/-179.0/-90.0/180.0",help="area defined as north/west/south/east with default 90.0/-179.0/-90.0/180.0", metavar="area")
    parser.add_option("--outputdir", dest="outputdir",
                      help="root directory for storing output files", metavar="outputdir")
    (options, args) = parser.parse_args()

    if not options.start_date:
        parser.error("start date must be specified!")
    else:
        start_date=options.start_date

    if not options.end_date:
        end_date=start_date
    else:
        end_date=options.end_date

    if not options.outputdir:
# if WORKDIR is defined, we will use it otherwise files 
# will be stored in the current directory
        outputdir=os.environ.get("WORKDIR",".")
    else:
        outputdir=options.outputdir


    print "start date %s "%(start_date)
    print "end date %s "%(end_date)

       
    server = ECMWFDataServer()

# Retrieve ERA interim data for running flexpart

    syear=int(start_date[:4])
    smonth=int(start_date[4:6])
    sday=int(start_date[6:])
    start = datetime.date( year = syear, month = smonth, day = sday )
    eyear=int(end_date[:4])
    emonth=int(end_date[4:6])
    eday=int(end_date[6:])

    end = datetime.date( year = eyear, month = emonth, day = eday )

    current_ym = ""
    ir_date = start
    retrieve="no"
    for date in daterange( start, end ):
# if new year & month then we create a new directory to store output files
        if date.strftime("%Y%m") != current_ym and current_ym != "":
               retrieve="yes"

        if date == end:
            retrieve="yes"

        if retrieve == "yes":
                # we need to retrieve MARS data for this period (maximum one month)
                flexpart = EIFlexpart()
                dates= ir_date.strftime("%Y%m%d") + "/to/" + er_date.strftime("%Y%m%d") 
                current_outputdir =  outputdir + "/"  + ir_date.strftime("%Y") + '/' + ir_date.strftime("%m") + '/' 
                mkdir_p(current_outputdir)
                print "retrieve " + dates + " in dir " + current_outputdir
                flexpart.retrieve(server, dates, options.times, options.area, options.levels, current_outputdir)
                ir_date = date
                retrieve="no"

        er_date = date

        current_ym =  date.strftime("%Y%m")



if __name__ == "__main__":
    main()
