#
# (C) Copyright 2014 UIO.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# 
# 
# Creation: July 2014 - Anne Fouilloux - University of Oslo
#
#

from gribapi import *
import traceback
import sys,os


##############################################################
#  GRIB utilities
##############################################################
class GribTools:
    'class for GRIB API with new methods'
    def __init__(self,filename):
        self.filename=filename
    
# get keyvalues for a given list of keynames
# a where statment can be given (list of key and list of values)

    def getkeys(self,keynames,wherekeynames=[],wherekeyvalues=[]):
        fileid=open(self.filename,'r')

        return_list=[]
        
        while 1:
            gid_in = grib_new_from_file(fileid)  
            if gid_in is None: break
            select=True
            i=0
            if len(wherekeynames) != len(wherekeyvalues): raise Exception("Give a value for each keyname!")

            for wherekey in wherekeynames:
                if not grib_is_defined(gid_in, wherekey): raise Exception("where Key was not defined")
                select=select and (str(wherekeyvalues[i])==str(grib_get(gid_in, wherekey)))
                i=i+1
            if select:
                llist = []
                for key in keynames:
                    llist.extend([str(grib_get(gid_in, key))])
                return_list.append(llist)
            grib_release(gid_in)
        fileid.close()
        return return_list
      
# set keyvalues for a given list of keynames
# a where statment can be given (list of key and list of values)
# an input file must be given as an input for reading grib messages
# note that by default all messages are written out
# if you want to get only those meeting the where statement, use
# strict=true
    def setkeys(self,fromfile,keynames,keyvalues, wherekeynames=[],wherekeyvalues=[], strict=False, filemode='w'):
        fout=open(self.filename,filemode)
        fin=open(fromfile)
        
        while 1:
            gid_in = grib_new_from_file(fin)  
            if gid_in is None: break
            
            select=True
            i=0
            if len(wherekeynames) != len(wherekeyvalues): raise Exception("Give a value for each keyname!")

            for wherekey in wherekeynames:
                if not grib_is_defined(gid_in, wherekey): raise Exception("where Key was not defined")
                select=select and (str(wherekeyvalues[i])==str(grib_get(gid_in, wherekey)))
                i=i+1
            if select:
                i=0
                for key in keynames:
                    grib_set(gid_in, key, keyvalues[i])
                    i=i+1
            if strict:
                if select:
                    grib_write(gid_in,fout)
            else:
                grib_write(gid_in,fout)
            grib_release(gid_in)
        fin.close()
        fout.close()

# Add the content of a grib file but only messages
# corresponding to keys/values        
# if selectWhere is False select fields that are different from keynames/keyvalues
    def copy(self,filename_in, selectWhere=True, keynames=[], keyvalues=[],filemode='w'):
        fin=open(filename_in)
        fout=open(self.filename,filemode)
        
        while 1:
            gid_in = grib_new_from_file(fin)  
            if gid_in is None: break
            
            select=True
            i=0
            if len(keynames) != len(keyvalues): raise Exception("Give a value for each keyname!")

            for key in keynames:
                if not grib_is_defined(gid_in, key): raise Exception("Key was not defined")
                
                if selectWhere:
                    select=select and (str(keyvalues[i])==str(grib_get(gid_in, key)))
                else:
                    select=select and (str(keyvalues[i])!=str(grib_get(gid_in, key)))
                i=i+1
            if select:
                grib_write(gid_in,fout)
            grib_release(gid_in)
        fin.close()
        fout.close()

# Create index from a list of files if it does not exist or read it
    def index(self,index_keys=["mars"], index_file = "my.idx"):
        self.iid = None
 
        if (os.path.exists(index_file)):
            self.iid = grib_index_read(index_file)
            print "Use existing index file: %s "%(index_file)
        else:
            for file in self.filename:
                print "Inputfile: %s "%(file)
                if self.iid == None:
                    self.iid = grib_index_new_from_file(file,index_keys)
                else:
                     grib_index_add_file(self.iid,file)

            if self.iid != None:
              grib_index_write(self.iid,index_file)
        return self.iid 


       

        
