#!/usr/bin/env python

''' This is a version of a script that executes the Data Reduction Pipeline
    from a folder one level above the images. The final goal is to have a script 
    connected to Yerkes' File Manager server that will automatically execute 
    whenever new images are uploaded, or on a set time interval (whichever is 
    more practical).

    This version of the pipeline is optimized for use with i/r/g-band images.
    It will still work if they are not there, but the order of images run through
    the pipeline may be incorrect.  If a different combination of filters are
    being used, the code can be edited appropriately (this version contains 'dummy'
    code that is currently commented out that can be easily edited/uncommented to
    work with any set of desired filters).
'''

import os
import sys
import logging

# Set logging format
logging.basicConfig(level = logging.INFO)

# Change directory & import the pipeline settings
sys.path.append('/data/scripts/DataReduction/source/')
from drp.pipeline import PipeLine
import datetime

today = datetime.date.today()
year = str(today.year)
month = str(today.month)
day = str(today.day)
if len(month) == 1:
    month = '0'+month
if len(day) == 1:
    day = '0'+day

date = year + '-' + month + '-' + day

datefilepath = '/data/images/StoneEdge/0.5meter/'+year+'/'+date

def execute():
    # Call the pipeline configuration
    pipe = PipeLine(config = '/data/scripts/DataReduction/pipeconf_master_edit_auto.txt')
    # This version only needs to be executed from a terminal. A specific image folder
    # (like the ones on the stars base) is specified for the pipeline.  The pipeline
    # will look in the folder and find any of the sub-folders that contain the FITS images.
    # It will then automatically take the files it finds and run them through the pipeline.
    print sys.argv
    if len(sys.argv) > 1:
        topdirectory = str(sys.argv[1])      # This is the specified directory -- entered as the second argument in the terminal command
    else:
        topdirectory = datefilepath
    rawlist = os.listdir(topdirectory)   # This is a list of everything in the specified directory.
    objectlist = []
    for Object in rawlist:
        if not '.' in Object:            # This line makes sure to exlude any stray files
            objectlist.append(Object)    # objectlist is simply a cleaned version of rawlist -- guaranteed to only contain the actual object folders from topdirectory
    print objectlist
    for entry in objectlist:                         # Run this loop for each object folder ('entry') found in objectlist
        imagelist = []                               # The empty list is created here so that it keeps getting rewritten for the pipeline (only 3 files at a time)
        fullentry = os.path.join(topdirectory,entry) # This line must be included for the pipeline to find the files
        for image in os.listdir(fullentry):          # Run this loop for each file ('image') found in the object folder ('fullentry')
	    if not '.fits' in image:                 # Makes sure the images collected are FITS images
                continue
            if 'dark' in image or 'flat' in image or 'bias' in image:    # Makes sure the images are regular filters
                continue
            imagelist.append(os.path.join(fullentry,image))   # Adds the correct images to imagelist
        print imagelist
        ilist = []
        rlist = []
        glist = []
        for File in imagelist:
            if 'i-band' in File or 'iband' in File or '%si' % entry.lower() in File:
                ilist.append(File)
            if 'r-band' in File or 'rband' in File or '%sr' % entry.lower() in File:
                rlist.append(File)
            if 'g-band' in File or 'gband' in File or '%sg' % entry.lower() in File:
                glist.append(File)
            else:
                continue
        print ilist
        print rlist
        print glist
        '''
        redlist = []                          # This is the section of code that can be edited/uncommented to optimize the pipeline for different filters
        greenlist = []                        # These lines create empty lists that organize the full imagelist into separate lists for each filter (only edit the list names if you want)
        bluelist = []
        for File in imagelist:                # These lines look for specified 'keywords' in all the files in imagelist to decide if they should be put in one of the filter lists
            if 'redkeyword' in File:          # The 'keyword' is what should be edited; in addition, several keywords can be used by adding "or 'keyword2' in File" etc.
                redlist.append(File)          # The names of the lists here must match the three empty lists that were just created.
            if 'greenkeyword' in File:
                greenlist.append(File)
            if 'bluekeyword' in File:
                bluelist.append(File)
            else:
                continue
        print redlist                         # These lines just print off all the contents of the separate filter lists; not actually necessary for operation
        print greenlist
        print bluelist
        '''
        # Now the program will run the files placed in (ilist, rlist, and glist) or imagelist through the pipeline.
        # It will do this for every entry in objectlist (i.e. for every object)
        # Allow for use with only one or two FITS files (won't result in normal color).
        pipe.reset()                   # Resets pipeline for running through multiple folders
        try:
            '''
            if len(redlist) >= 1 and len(greenlist) >= 1 and len(bluelist) >= 1:             # This is section of code is used to optimize pipeline for other filter combinations
                result = pipe([redlist[0], greenlist[0], bluelist[0]])                       # red/green/bluelist variables should be the same as those in the previous section of edited code
            '''
            if len(ilist) >= 1 and len(rlist) >= 1 and len(glist) >= 1:                      # This occurs if there is at least one i-, r-, and g-band filter found in imagelist (optimal circumstance)
                result = pipe([ilist[0], rlist[0], glist[0]])                                # In this case, the first image from each filter list will be put through the pipeline in the correct order.
            elif len(imagelist) == 1:
                result = pipe([imagelist[0], imagelist[0], imagelist[0]])                    # Cases where there is only one image in imagelist
            elif len(imagelist) == 2:
                result = pipe([imagelist[0], imagelist[1], imagelist[1]])                    # Cases where there are only two images in imagelist
            elif len(imagelist) >= 3:
                result = pipe([imagelist[0], imagelist[1], imagelist[2]])                    # Cases where there are three or more images in imagelist, but not one of each filter.
            else:
                print "***ERROR: Invalid entries. No files found in folder: %s***" %entry    # No images found in imagelist
        except:
            print "***ERROR: Issue occured while reducing files in folder: %s***" %entry     # Other miscellaneous error, written so that the pipeline won't be killed by such errors.

execute()
''' HISTORY '''

''' 2014/07/29: Original script created by Neil Stilin. Only used in local folder.
    2014/08/07: Added code to allow for use with only one or two FITS files -NS
    2015/01/05: Edited code to automatically find and process all image files from a specified multi-level directory  -- Neil S.
    2015/01/07: This version of the automatic pipeline should be used to reduce all images taken on a specific day.  -- Neil S.
    2015/02/23: This version can be executed without system arguments to process files from the current day -- Ben Mahon
'''
