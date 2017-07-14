#!/usr/bin/env python

''' This is a version of a script that executes the Data Reduction Pipeline
    from a folder two levels above the images. The final goal is to have a script 
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
import traceback

# Set system variables
logfile = '/data/scripts/DataReduction/PipeLineLog.txt'

# Set logging format
logging.basicConfig(filename = logfile, level = logging.DEBUG, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
log = logging.getLogger('pipe.ExecuteAutoYear')
log.info('Starting up')


# Change directory & import the pipeline settings
sys.path.append('/data/scripts/DataReduction/source/')
from drp.pipeline import PipeLine

def execute():
    # Call the pipeline configuration
    pipe = PipeLine(config = '/data/scripts/DataReduction/pipeconf_stonedge_auto.txt')
    # This version only needs to be executed from a terminal. A specific image folder
    # (like the ones on the stars base) is specified for the pipeline.  The pipeline
    # will look in the folder and find any of the sub-folders that contain the FITS images.
    # It will then automatically take the files it finds and run them through the pipeline.
    topdirectory = str(sys.argv[1])    # This is the specified directory  -- determined by the second argument in the command line
    rawlist = os.listdir(topdirectory) # This is a list of everything in the specified directory.
    datelist = []
    for date in rawlist:
        if not '.' in date:            # This line makes sure to exlude any stray files
            datelist.append(date)      # datelist is simply a cleaned version of rawlist -- and only contains the actual date folders from topdirectory
    objectlist = []
    for day in datelist:                                         # Run this loop for every folder ('day') found in datelist
        fullday = os.path.join(topdirectory,day)                 # Full directory path is required for the os.listdir to function properly
        for Object in os.listdir(fullday):                       # Run loop for each object folder found in each of the days in datelist
            fullobject = os.path.join(topdirectory,day,Object)   # Full directory path is needed (again)
            if not '.' in Object:                                # This line makes sure to exlude any stray files
                objectlist.append(fullobject)       # Adds all object folders found in each day folder to objectlist
    log.info('Object list = %s' %repr(objectlist))
    for entry in objectlist:                        # Run this loop for each object folder ('entry') found in objectlist
        imagelist = []                              # The empty list is created here so that it keeps getting rewritten for the pipeline
        fullentry = os.path.join(fullobject,entry)  # This line must be included for the pipeline to find the files
        for image in os.listdir(fullentry):         # Run this loop for each file ('image') found in the object folder ('fullentry')
	    num = image[-14:-5]
	    if not 'seo.fits' in image[-8:]:                # Makes sure the images collected are FITS images, not KEYS or WCS
                if not 'seo%s.fits' % num in image[-17:]:
                    continue
            if 'dark' in image or 'flat' in image or 'bias' in image:    # Ignore dark, flat and bias files 
                continue
            imagelist.append(os.path.join(fullentry,image))   # Adds the correct images to imagelist
        log.info('Object = %s Image list = %s' % (entry, repr(imagelist)))
        if len(imagelist) == 0 :
            log.warning('Image List is Empty, skipping object = %s' % entry)
            continue
        # Now the program will run the files placed in imagelist through the pipeline.
        # It will do this for every entry in objectlist (i.e. for every object)
        pipe.reset()
        # Run the pipeline (return with error message)
        try:
            result = pipe(imagelist)
        except:
            log.warning("Pipeline for object = %s returned Error" % entry)

# Run the setup code in an error with reporting traceback
try:
    execute()
except Exception, e:
    log.error('Found Error = %s' % repr(e))
    trb = traceback.format_exc().split('\n').reverse()
    for tr in trb:
        log.error(tr)
    raise e

def confirmation():
    response = raw_input("Are you sure you wish to reduce an entire year's worth of images? yes/no: ")
    while True:
        if response == "Yes" or response == 'yes':
            print "Alright, here we go. You may want to go for a run, this is going to take a while."
            execute()
            break
        elif response == "No" or response == "no":
            print "Aborting reduction. Thanks for saving me a ton of work :)"
            break
        else:
            response = raw_input("Invalid entry. Please type \"yes\" or \"no\": ")
            continue

confirmation()

''' HISTORY '''

''' 2017/07/13: This version processes all inputs into the pipeine instead of just 3 inputs, so StepMakeRGB can get the desired  3 best inputs for a jpgimage -- Atreyo Pal
    2014/07/29: Original script created by Neil Stilin. Only used in local folder.
    2014/08/07: Added code to allow for use with only one or two FITS files -NS
    2015/01/05: Edited code to automatically find and process all image files from a specified multi-level directory  -- Neil S.
    2015/01/07: This version should be used to reduce all images taken in a given year; added a confirmation function to remove possiblity of unintentional use.  -- Neil
'''
