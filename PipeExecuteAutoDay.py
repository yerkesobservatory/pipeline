#!/usr/local/bin/python3

# Below is the "default" python path, the one above is necessary on stars.
#!/usr/bin/env python

''' This is a version of a script that executes the Data Reduction Pipeline
    from a folder one level above the images. The final goal is to have a script 
    connected to Yerkes' File Manager server that will automatically execute 
    whenever new images are uploaded, or on a set time interval (whichever is 
    more practical).

    This version of the pipeline is optimized for use with i/r/g-band images.
    It will still work if they are not there, but the order of images run through
    the pipeline may be incorrect.  If a different combination of filters are
    being used, the code can be edited appropriately (this version contains
    'dummy' code that is currently commented out that can be easily
    edited/uncommented to work with any set of desired filters).
'''

import os
import sys
import logging
import traceback
import datetime
import re

# Set system variables
logfile = '/data/scripts/DataReduction/PipeLineLog.txt'
#logfile = '/Users/atreyopal/Desktop/pipeline/PipeLineLog.txt'

# Set logging format
logging.basicConfig(filename = logfile, level = logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
log = logging.getLogger('pipe.ExecuteAutoDay')
log.info('Starting up')

# Change directory & import the pipeline settings
#sys.path.append('/Users/atreyopal/Desktop/pipeline/source/')
sys.path.append('/data/scripts/DataReduction/source/')
from darepype.drp.pipeline import PipeLine

today = datetime.date.today()
year = str(today.year)
month = str(today.month)
day = str(today.day)
if len(month) == 1:
    month = '0'+month
if len(day) == 1:
    day = '0'+day

date = year + '-' + month + '-' + day

#datefilepath = '/Users/atreyopal/Desktop/pipeline/Examples'
datefilepath = '/data/images/StoneEdge/0.5meter/'+year+'/'+date

def execute():
    # Call the pipeline configuration
    #pipe = PipeLine(config = '/Users/atreyopal/Desktop/pipeline/pipeconf_stonedge_remote.txt')
    pipe = PipeLine(config = ['/data/scripts/DataReduction/config/pipeconf_SEO.txt',
                              '/data/scripts/DataReduction/config/dconf_stars.txt'])
    # This version only needs to be executed from a terminal. A specific image folder
    # (like the ones on the stars base) is specified for the pipeline.  The pipeline
    # will look in the folder and find any of the sub-folders that contain the FITS images.
    # It will then automatically take the files it finds and run them through the pipeline.
    print(sys.argv)
    if len(sys.argv) > 1:
        # Load the specified directory -- entered as the second argument in the terminal command
        topdirectory = str(sys.argv[1])
    else:
        topdirectory = datefilepath
    # Load a list of everything in the specified directory.
    rawlist = os.listdir(topdirectory)
    objectlist = []
    # make objectlist as a cleaned version of rawlist
    # -- guaranteed to only contain the actual object folders from topdirectory
    for Object in rawlist:
        if not '.' in Object:            # This line makes sure to exlude any stray files
            if not 'itzamna' in Object:  # Exclude unsorted files
                objectlist.append(Object)
    log.info('Object list = %s' %repr(objectlist))
    # Run this loop for each object folder ('entry') found in objectlist
    # THIS IS THE MAIN LOOP OVER ALL OBSERVED OBJECTS
    for entry in objectlist:
        # The empty list is created here so that it keeps getting rewritten for the pipeline
        imagelist = []
        # Add the full path the the path for the observation
        fullentry = os.path.join(topdirectory,entry)
        # Run this loop for each file ('image') found in the object folder ('fullentry')
        for image in os.listdir(fullentry):
            # Makes sure the images collected are FITS images, i.e. avoid reducing
            #     existing data reduction product
            # This searches image filenames which end with .fits preceeded by
            #     '_0', 'RAW', 'seo' or 'SRT' with optional '.gz' at the end
            if not re.search(r'(_0|RAW|seo|SRT)\.fits(?:\.gz)?\Z',image):
                continue
            # Ignore dark, flat or bias images
            if 'dark' in image or 'flat' in image or 'bias' in image or 'pinpoint' in image:
                continue
            # Adds the correct images to imagelist
            imagelist.append(os.path.join(fullentry,image))
        log.info('Object = %s Image list = %s' % (entry, repr(imagelist)))
        if len(imagelist) == 0 :
            log.warning('Image List is Empty, skipping object = %s' % entry)
            continue
        # Now the program will run the files placed in imagelist through the pipeline.
        # It will do this for every entry in objectlist (i.e. for every object)
        pipe.reset()
        # Run the pipeline (return with error message)
        #result = pipe(imagelist)
        try:
            #pass
            result = pipe(imagelist, pipemode = 'seo_server',force=True)
        except Exception as e:
            log.warning("Pipeline for object = %s returned Error" % entry)
            log.warning('Found Error = %s' % repr(e))
            print('Found Error = %s' % repr(e))
            try:
                trb = traceback.format_exc().split('\n').reverse()
                for tr in trb:
                    log.warning(tr)
            except Exception as f:
                log.warning('Unable to print traceback')
                print(traceback.format_exc(),trb)


# Run the setup code in an error with reporting traceback
#execute()
try:
    #pass
    execute()
except Exception as e:
    log.error('Found Error = %s' % repr(e))
    try:
        trb = traceback.format_exc().split('\n').reverse()
        for tr in trb:
            log.error(tr)
    except Exception as f:
        log.warning('Unable to print traceback')
        print(traceback.format_exc(),trb)
    raise e

''' 
HISTORY:
2017/06/23: This version processes all inputs into the pipeine instead of just 3
            inputs, so StepMakeRGB can get the desired  3 best inputs for a jpg
            image -- Atreyo Pal
2015/02/23: This version can be executed without system arguments to process
            files from the current day -- Ben Mahon
2015/01/07: This version of the automatic pipeline should be used to reduce
            all images taken on a specific day.  -- Neil S.
2015/01/05: Edited code to automatically find and process all image files from
            a specified multi-level directory  -- Neil S.
2014/08/07: Added code to allow for use with only one or two FITS files -NS
2014/07/29: Original script created by Neil Stilin. Only used in local folder.
'''

'''
        OLD CODE: for the case we want different red, green and blue images
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
