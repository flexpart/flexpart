#!/usr/bin/env python
#
# (C) Copyright 2014 UIO.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# 
# Creation: October 2014 - Anne Fouilloux - University of Oslo
#
# Compute WINDFILEDS from ECMWF GRIB data
# 
import calendar
import shutil
import datetime
import time
import os
from UIOTools import UIOFiles
from string import strip
from optparse import OptionParser
from GribTools import GribTools
from FlexpartTools import EIFlexpart, daterange


def main():
    usage = """usage: %prog --start_date=YYYYMMDD [--end_date=YYYYMMDD] [--namelist=namelist_for_convert] [--inputdir=input_root_directory] [--outputdir=output_directory] """
    parser = OptionParser(usage=usage)
    parser.add_option("--start_date", dest="start_date",
                      help="start date YYYYMMDD", metavar="start_date" )
    parser.add_option( "--end_date", dest="end_date",
                      help="end_date YYYYMMDD", metavar="end_date")
    parser.add_option("--namelist", dest="namelist",
                      help="namelist used for converting", metavar="namelist")
    parser.add_option("--inputdir", dest="inputdir",
                      help="root directory for reading input files", metavar="inputdir")
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

    if not options.namelist:
        namelist='fort.4'
    else:
        namelist=options.namelist

    if not options.inputdir:
# if WORKDIR is defined, we will use it otherwise files 
# will be stored in the current directory
        inputdir=os.environ.get("WORKDIR",".")
    else:
        inputdir=options.inputdir

    if not options.outputdir:
# if FLEXPART_WINDS is defined, we will use it otherwise files 
# will be stored in the current directory
        outputdur=os.environ.get("FLEXPART_WINDS",".")
    else:
        outputdir=options.outputdir



    syear=int(start_date[:4])
    smonth=int(start_date[4:6])
    sday=int(start_date[6:])
    start = datetime.date( year = syear, month = smonth, day = sday )
    eyear=int(end_date[:4])
    emonth=int(end_date[4:6])
    eday=int(end_date[6:])

    end = datetime.date( year = eyear, month = emonth, day = eday )



    cyear = -1
    cmont = -1
    inputfiles=UIOFiles(['.grib', '.grb', '.grib1', '.grib2', '.grb1','.grb2'])
    if (not os.path.exists('fort.4')):
        fnamelist=open('fort.4','w')
        shutil.copyfileobj(open(namelist,'r'), fnamelist)
        fnamelist.close()
    for date in daterange( start, end ):
# data retrieved by year/month
           if cyear != date.year or cmonth != date.month:
             print 'year : ' + str(date.year) + ' month : ', date.month
             cyear = date.year
             cmonth = date.month

             # we will make the list of files from the root inputdir
             if cmonth < 10:
                 inputfiles.listFiles(inputdir + '/'+str(cyear)+'/0'+str(date.month))
             else:
                 inputfiles.listFiles(inputdir + '/'+str(cyear)+'/'+str(date.month))
                 
             flexpart = EIFlexpart()
             flexpart.create(inputfiles, outputdir)

if __name__ == "__main__":
    main()
