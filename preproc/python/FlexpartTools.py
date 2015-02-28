#
# (C) Copyright 2014 UIO.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# 
# Creation: October 2014 - Anne Fouilloux - University of Oslo
#
#
import subprocess
import shutil
import os, errno
import datetime
from dateutil.relativedelta import relativedelta
import re

from string import atoi
from numpy  import *
from ecmwfapi import ECMWFDataServer
from gribapi import *
from GribTools import GribTools

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

###############################################################
# utility to remove a file if it exists
# it does not fail if the file does not exist
###############################################################
def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured

###############################################################
# Utility to create directories
###############################################################
def mkdir_p(file):
    path = os.path.dirname(file)
    try:
        os.stat(path)
    except:
        os.makedirs(path)

##############################################################
# function to iterate over dates
##############################################################
def daterange( start_date, end_date ):
    if start_date <= end_date:
        for n in range( ( end_date - start_date ).days + 1 ):
            yield start_date + datetime.timedelta( n )
    else:
        for n in range( ( start_date - end_date ).days + 1 ):
            yield start_date - datetime.timedelta( n )

##############################################################
# function to iterate over years
##############################################################
def years_between(start_date,end_date):
    years = []
    cursor = start_date

    while cursor <= end_date:
        if cursor.year not in years:
            years.append(cursor.year)
        cursor += relativedelta(years=1)

    return years
##############################################################
# function to iterate over months
##############################################################
def months_between(start_date,end_date):
    months = []
    cursor = start_date

    while cursor <= end_date:
        if cursor.month not in months:
            months.append(cursor.month)
        cursor += datetime.timedelta(days=1)

    return months

##############################################################
# MARSretrieval class
##############################################################
class MARSretrieval:
    'class for MARS retrievals'

    def __init__(self,server, dataset="interim",marsclass="ei",type="",levtype="",levelist="", 
                 repres="", date="",resol="",stream="",area="",time="",step="",expver="1",
                 number="",accuracy="", grid="", gaussian="", target="",param=""):
        self.dataset=dataset
        self.marsclass=marsclass
        self.type=type
        self.levtype=levtype
        self.levelist=levelist
        self.repres=repres
        self.date=date
        self.resol=resol
        self.stream=stream
        self.area=area
        self.time=time
        self.step=step
        self.expver=expver
        self.target=target
        self.param=param
        self.number=number
        self.accuracy=accuracy
        self.grid=grid
        self.gaussian=gaussian
        self.server=server

    def displayInfo(self):
        if self.dataset:
            print "dataset: ", self.dataset
        if self.marsclass:
            print "class: ", self.marsclass
        if self.type:
             print "type: ", self.type
        if self.levtype:
            print "levtype: ", self.levtype
        if self.levelist:
            print "levelist: ", self.levelist
        if self.repres:
            print "repres: ", self.repres
        if self.date:
            print "date: ", self.date
        if self.resol:
            print "resol: ", self.resol
        if self.stream:
            print "stream: ", self.stream
        if self.area:
            print "area: ", self.area
        if self.time:
            print "time: ", self.time
        if self.step:
            print "step: ", self.step
        if self.expver:
            print "expver: ", self.expver
        if self.target:
            print "target: ", self.target
        if self.param:
            print "param: ", self.param
        if self.number:
            print "number: ", self.number
        if self.accuracy:
            print "accuracy: ", self.accuracy
        if self.grid:
            print "grid: ", self.grid
        if self.gaussian:
            print "gaussian: ", self.gaussian

    def dataRetrieve(self):
        dicolist = {}
        if self.dataset:
	        dicolist['dataset']="%s"%(self.dataset)
        if self.marsclass:
	        dicolist['class']="%s"%(self.marsclass)
        if self.date:
	        dicolist['date']="%s"%(self.date)
        if self.time:
	        dicolist['time']="%s"%(self.time)
        if self.expver:
	        dicolist['expver']="%s"%(self.expver)
        if self.param:
	        dicolist['param']="%s"%(self.param)
        if self.type:
	        dicolist['type']="%s"%(self.type)
        if self.levtype:
	        dicolist['levtype']="%s"%(self.levtype)
        if self.levelist:
	        dicolist['levelist']="%s"%(self.levelist)
        if self.repres:
	        dicolist['repres']="%s"%(self.repres)
        if self.resol:
	        dicolist['resol']="%s"%(self.resol)
        if self.stream:
	        dicolist['stream']="%s"%(self.stream)
        if self.area:
	        dicolist['area']="%s"%(self.area)
        if self.step:
	        dicolist['step']="%s"%(self.step)
        if self.target:
	        dicolist['target']="%s"%(self.target)
        if self.number:
	        dicolist['number']="%s"%(self.number)
        if self.accuracy:
	        dicolist['accuracy']="%s"%(self.accuracy)
        if self.grid:
	        dicolist['grid']="%s"%(self.grid)
        if self.gaussian:
	        dicolist['gaussian']="%s"%(self.gaussian)

        self.server.retrieve(dicolist)


##############################################################
class EIFlexpart:
    'class to retrieve Era Interim data for running FLEXPART'
##############################################################

    def __init__(self):
# different mars types for retrieving reanalysis data for flexpart
       self.types=["an","fc"]

       self.mars={}
# set type (an/fc) and times and steps for EI
#         (real) time    [type     time  step ] as in mars  (step for forecast only)     
       self.mars["00"] = ["an",    "00"]
       self.mars["01"] = ["fc",    "00", "01"]
       self.mars["02"] = ["fc",    "00", "02"]
       self.mars["03"] = ["fc",    "00", "03"]
       self.mars["04"] = ["fc",    "00", "04"]
       self.mars["05"] = ["fc",    "00", "05"]
       self.mars["06"] = ["an",    "06"]
       self.mars["07"] = ["fc",    "00", "07"]
       self.mars["08"] = ["fc",    "00", "08"]
       self.mars["09"] = ["fc",    "00", "09"]
       self.mars["10"] = ["fc",    "00", "10"]
       self.mars["11"] = ["fc",    "00", "11"]
       self.mars["12"] = ["an",    "12"]
       self.mars["13"] = ["fc",    "12", "01"]
       self.mars["14"] = ["fc",    "12", "02"]
       self.mars["15"] = ["fc",    "12", "03"]
       self.mars["16"] = ["fc",    "12", "04"]
       self.mars["17"] = ["fc",    "12", "05"]
       self.mars["18"] = ["an",    "18"]
       self.mars["19"] = ["fc",    "12", "07"]
       self.mars["20"] = ["fc",    "12", "08"]
       self.mars["21"] = ["fc",    "12", "09"]
       self.mars["22"] = ["fc",    "12", "10"]
       self.mars["23"] = ["fc",    "12", "11"]
       
    def retrieve(self, server, dates,times, area, levels, outputdir):
        self.dates=dates
        self.times=times
        self.server=server
        self.area=area
        self.levels=levels
        if outputdir=="":
	  self.outputdir='.'
        else:
          self.outputdir=outputdir
        # surface data
        dataset="interim"
        marsclass="ei"
        expver="0001"
        levels="1/to/"+ self.levels
# Retrieve Q not for using Q but as a template for a reduced gaussian grid one date and time is enough
# Take analysis at 00
        qdate=self.dates
        idx=qdate.find("/")
        if (idx >0):
          qdate=self.dates[:idx]
        
        Q= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type="an", levtype="ML", levelist="1",
                         gaussian="reduced",grid="80", resol="159",accuracy="24",target=self.outputdir+"/"+"Q.grb",
                         date=qdate, time="00",expver=expver, param="133.128")
        Q.displayInfo()
        Q.dataRetrieve()

        for type in self.types:
            stime=set()
            sstep=set()
            for time in self.times.split('/'):
# not so nice... to review later!
                if self.mars[time][0] == type:
                  stime.add(self.mars[time][1])
                  if (len(self.mars[time]) > 2):
                     sstep.add(self.mars[time][2])
          
            rtime="/".join(stime)
            if len(sstep) > 0:
              rstep="/".join(sstep)
            else:
              rstep=""
            print type
            print rtime
            if rstep!="":
               print rstep
            print "MARS retrieve start... "
            paramUVD="131.128/132.128/155.128"
            paramT="130.128"
            paramLNSP="152.128"
            param2D2T="167.128/168.128"
            if rstep!="":
                UVD= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist=levels,
                         resol="159", accuracy="24",grid="OFF",target=self.outputdir+"/"+type+"UVD.grb",
                         area=self.area, date=self.dates, time=rtime,step=rstep, expver=expver, param=paramUVD)
                T= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist=levels,
                         resol="159", accuracy="24",grid="1.0/1.0",target=self.outputdir+"/"+type+"T.grb",
                         area=self.area, date=self.dates, time=rtime,step=rstep, expver=expver, param=paramT)
                LNSP= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist="1",
                         accuracy="24",target=self.outputdir+"/"+type + "LNSP.grb",
                         area=self.area, date=self.dates, time=rtime,step=rstep, expver=expver, param=paramLNSP)
                DT= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="SFC", levelist="OFF",
                         resol="159",accuracy="24",grid="1.0/1.0",target=self.outputdir+"/"+type + "2D2T.grb",
                         area=self.area, date=self.dates, time=rtime,step=rstep, expver=expver, param=param2D2T)
                print "do not retrieve Q for forecast"
            else:
                UVD= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist=levels,
                         resol="159", accuracy="24",grid="OFF",target=self.outputdir+"/"+type + "UVD.grb",
                         area=self.area, date=self.dates, time=rtime,expver=expver, param=paramUVD)
                T= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist=levels,
                         resol="159", accuracy="24",grid="1.0/1.0",target=self.outputdir+"/"+type + "T.grb",
                         area=self.area, date=self.dates, time=rtime,expver=expver, param=paramT)
                LNSP= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="ML", levelist="1",
                         accuracy="24",target=self.outputdir+"/"+type + "LNSP.grb",
                         area=self.area, date=self.dates, time=rtime,expver=expver, param=paramLNSP)
                DT= MARSretrieval(self.server, dataset=dataset, marsclass=marsclass, stream="oper", type=type, levtype="SFC", levelist="OFF",
                         resol="159",accuracy="24",grid="1.0/1.0",target=self.outputdir+"/"+type + "2D2T.grb",
                         area=self.area, date=self.dates, time=rtime,expver=expver, param=param2D2T)

            UVD.displayInfo()
            UVD.dataRetrieve()
            T.dataRetrieve()
            LNSP.dataRetrieve()
            DT.dataRetrieve()
        print "MARS retrieve done... "
        
    def getFlexpartTime(self, type,step, time):
        if int(step) < 10:
            cstep = '0' +  str(step)
        else:
            cstep = str(step)
        if int(time/100) < 10:
            ctime = '0' +  str(int(time/100))
        else:
            ctime = str(int(time/100))

        ctype = str(type)
        if ctype == 'fc':
            myinfo = [ctype,ctime, cstep]
        else:
            myinfo = [ctype,ctime]
        cflextime = None
        for t, marsinfo in self.mars.items():
            if myinfo == marsinfo:
                cflextime=t
        return cflextime

    def create(self, inputfiles, outputdir):
        index_keys=["date","time","stepRange"]
        indexfile="date_time_stepRange.idx"
        silentremove(indexfile)
        grib=GribTools(inputfiles.files)
        iid=grib.index(index_keys=index_keys, index_file = indexfile)

        print 'index done...' 
        silentremove("fort.10")
        silentremove("fort.11")
        silentremove("fort.12")
        silentremove("fort.13")
        silentremove("fort.18")
        index_vals = []
        for key in index_keys:
            key_vals = grib_index_get(iid,key)

            index_vals.append(key_vals)


        for prod in product(*index_vals):
            for i in range(len(index_keys)):
                grib_index_select(iid,index_keys[i],prod[i])

            f10 = open('fort.10','w')
            f11 = open('fort.11','w')
            f12 = open('fort.12','w')
            f13 = open('fort.13','w')
            f16 = open('fort.16','w')
            gid = grib_new_from_index(iid)
            hid = gid
            cflextime = None
            if gid is not None: 
                cdate = str(grib_get(gid, 'date'))
                cyear = cdate[:4]
                cmonth = cdate[4:6]
                cday = cdate[6:8]
                time = grib_get(gid, 'time')
                type = grib_get(gid, 'type')
                step = grib_get(gid, 'stepRange')
                cflextime = self.getFlexpartTime(type,step, time)
#                print 'cyear '+cyear+'/'+cmonth+'/'+'/EI'+cyear[2:4]+cmonth+cday+cflextime
            while 1:  
                if gid is None: break
                paramId = grib_get(gid, 'paramId')
                if paramId == 133:
# Relative humidity (Q.grb) is used as a template only so we need the first we "meet"
                    fout=open('fort.18','w')
                    grib_write(gid,fout)
                    fout.close()
                if paramId == 131 or paramId == 132:
                    grib_write(gid, f10)
                if paramId == 130:
                    grib_write(gid,f11)
                if paramId == 152:
                    grib_write(gid,f12)
                if paramId == 155:
                    grib_write(gid,f13)
                if paramId == 167 or paramId == 168:
                    grib_write(gid, f16)                    
                grib_release(gid)
                gid = grib_new_from_index(iid)
# call for CONVERT2
            
            f10.close()
            f11.close()
            f12.close()
            f13.close()
            f16.close()
            if hid is not None:
                p=subprocess.check_call(['CONVERT2'])
# create the corresponding output file  fort.15 (generated by CONVERT2) + fort.16 (paramId 167 and paramId 168)  
                mkdir_p(outputdir+'/'+cyear+'/'+cmonth+'/')
                print "outputdir = " + outputdir+'/'+cyear+'/'+cmonth+'/'+'/EI'+cyear[2:4]+cmonth+cday+cflextime
                fout = open(outputdir+'/'+cyear+'/'+cmonth+'/EI'+cyear[2:4]+cmonth+cday+cflextime,'wb')    
                shutil.copyfileobj(open('fort.15','rb'), fout)     
                shutil.copyfileobj(open('fort.16','rb'), fout)
                fout.close()

        grib_index_release(iid)

    def __del__(self):
        print "clean"
        silentremove("fort.10")
        silentremove("fort.11")
        silentremove("fort.12")
        silentremove("fort.13")
        silentremove("fort.15")
        silentremove("fort.16")
        silentremove("fort.18")
        silentremove("VERTICAL.EC")
        silentremove("date_time_stepRange.idx")

