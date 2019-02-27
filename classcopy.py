#!/usr/local/bin/python

# Below is the "default" python path, the one above is necessary on stars
#!/usr/bin/env python

""" Class Copy
    ==========

    Copy the files in the astroclass folder to the correct date folder in the
    main data depository.

    i.e. the folder
        /astroclass/Marc/M42_..._2018-01-01_...Marc
    will be copied to 
        /StoneEdge/0.5meter/2018/2018-01-01

    Input:
    - A User folder ex: /astroclass/Marc
    
    Steps:
    - Get all subfolders
    - For all files, set the correct observer name
    - Check if any fits files are at the wrong level
    - Copy all folders and fits files to the correct date folder
"""

### Settings
# Logfile
logfile = '/data/scripts/DataReduction/PipeLineLog.txt'
# pythonpath
pypath = '/data/scripts/DataReduction/source'
# Pipeline configuration file
pipeconf = '/data/scripts/DataReduction/pipeconf_stonedge_auto.txt'
# Folder for image database (for images to end up in)
#     program will add yyyy/yyyy-mm-dd
databasefolder = '/data/images/StoneEdge/0.5meter'

### Preparation
# Imports
import sys
import os
sys.path.append(pypath)
import re
import shutil
import logging
from drp.datafits import DataFits

# Set up logging
logging.basicConfig(filename = logfile, level = logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
log = logging.getLogger('ClassCopy')
log.addHandler(logging.StreamHandler())

# Get input folders
if len(sys.argv) > 1:
    input_folder = sys.argv[1] 
else:
    print("""Usage: python classcopy.py input_folder
          where input_folder is the user folder
          ex: /images/astroclass/Berthoud""")
    exit()
log.info('Starting up with folder %s' % input_folder)
# Get observer name (last part of input_folder)
obsname = os.path.split(input_folder)[1]
log.debug('Observer Name = %s' % obsname)

### Get all observation folders 
#   Also make list of misplaced files
obslist = []
missfiles = []
rawlist = os.listdir(input_folder)
for fname in rawlist:
    if os.path.isdir(os.path.join(input_folder,fname)):
        obslist.append(fname)
    else:
        missfiles.append(fname)

### Loop through observation folders
folderfail = []
for obsfolder in obslist:
    log.debug('Start observation folder %s' % obsfolder)
    obspath = os.path.join(input_folder, obsfolder)
    # get data file list
    flist = [f for f in os.listdir(obspath) if '.fits' in f]
    # Go through Files:
    for f in flist:
        # Set observer
        df = DataFits(config = pipeconf)
        df.load(os.path.join(obspath,f))
        fobserver = df.getheadval('OBSERVER')
        if fobserver in ['Remy Prechelt','Matt Nowinkski']:
            log.debug('Set Observer for file %s' % f)
            df.setheadval('OBSERVER', obsname)
            df.save()
        # Check if it's in missfiles
        if f in missfiles:
            missfiles.remove(f)
            os.remove(os.path.join(input_folder,f))
            log.debug('File %s found, removed from input folder' % f)
    # get date in yyyy-mm-dd format
    m = re.search(r'\d+-\d+-\d+',obsfolder)
    dat = m.group()    
    log.debug('Date = %s' % dat)
    # copy folder and files to correct date
    tarfolder = os.path.join(databasefolder,dat[:4],dat,obsfolder)
    log.debug('Copying %s -> %s' % (obspath, tarfolder))
    try:
        shutil.copytree(obspath,tarfolder)
    except:
        log.warn('Unable to copy %s to %s' % (obspath, tarfolder) )
        folderfail.append(obspath)

# Output Missing Files
for f in missfiles:
    log.warn('File %s is missplaced and not found in subfolders' % f)
for f in folderfail:
    log.warn('Folder %s can not be copied' % f)
