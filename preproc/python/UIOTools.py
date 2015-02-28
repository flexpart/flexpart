#!/usr/bin/env python                                                                                                                                                                      
#                                                                                                                                                                                          
# (C) Copyright 2014 UIO.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# Creation: October 2014 - Anne Fouilloux - University of Oslo
#

import os

class UIOFiles:
    'class to manipulate files'
    def __init__(self,suffix):
# type of files to manipulate such as ['.grib', 'grb', 'grib1', 'grib2', 'grb1','grb2']
        self.suffix=suffix

    def listFiles(self,pathname):
        ''' list files (suffix previously given) within this directory. '''
    # Get the absolute path of the pathname parameter
        pathname = os.path.abspath(pathname)
 
    # Get a list of files in pathname
        filesInCurDir = os.listdir(pathname)
        self.counter = 0
        self.files = []
    # Traverse through all files
        for file in filesInCurDir:
            curFile = os.path.join(pathname, file)
 
        # Check if it's a normal file or directory
            if os.path.isfile(curFile):
                # Get the file extension
                fileNoExt,curFileExtension = os.path.splitext(curFile)
                # Check if the file has an extension of typical video files
                if curFileExtension in self.suffix:
                    # We have got a file file! Increment the counter
                    self.counter += 1
                    # add this filename in the list
                    self.files.append(curFile)
                    
            else:
                # We got a directory, enter into it for further processing
                self.listFiles(curFile)
