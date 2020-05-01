#!/usr/local/bin/python3

# Below is the "default" python path, the one above is necessary on stars
#!/usr/bin/env python

""" QUEUE COPY
    ==========

    This program copies Stoneedge data from the queue folder on stars:
    - Science raw data gets copied into the data folder and if not available
      .RAW.fits is attached to the filename.
    - Flats / Bias / Dark data is copied to the regular data folder if there
      is no duplicate already there

    Usage:
      queuecopy.py source_folders
 
      source_folders: folders where the data to copy is (todays folder if omitted)
                      multiple entries (or globs) are possible

"""

### Settings
# inpath path combined with current date yyyy-mm-dd to grep for folders
#   this is only used for the today option
inpath = '/data/public/queue/*/*%s*'
# bias dark flat folder: folder below which yyyy-mm-dd/flat folders are
bdfpath = '/data/images/StoneEdge/0.5meter/2018' 
# output path: folder below which User/Observation_YYMMDD/rawfile.RAW.fits are copied
outpath = '/data/images/queue'
# piperunpath: folder for the piperun files
piperunpath = '/data/images/queue/A_Test/piperuns'
# pythonpath
pypath = '/data/scripts/DataReduction/source'

### Preparation
# Imports
import os
import sys
sys.path.append(pypath)
import logging
import glob
import time
import string
import re
import shutil
from darepype.drp.datafits import DataFits

# Set up logging
logging.basicConfig(level = logging.DEBUG)
log=logging.getLogger('QueueCopy')
log.setLevel(logging.DEBUG)
log.info('Logging Set Up')

# Read Source Folders
if len(sys.argv) > 1:
    if sys.argv[1].lower() == 'today':
        source_folders = inpath % time.strftime("%Y-%m-%d", time.localtime())
        source_folders = glob.glob(source_folders)
    else:
        source_folders = sys.argv[1:]
else:
    source_folders = [os.getcwd()]
log.debug('Source Folder = %s' % repr(source_folders))

### Loop over source folders
for source_folder in source_folders:

    ### Copy Bias / Dark / Flats - Has been deactivated b/c filenameformat changed.
    for ftype in []: #['bias', 'dark', 'flat']:
        # Get all fitting files: Files must have something like */flat/*flat*.fits
        globstr = os.path.join(source_folder,'*/%s/*%s*.fits' % (ftype, ftype) )
        flist = glob.glob(globstr)
        # Glob target folder for file type
        typefiles = [ f for f in glob.glob(os.path.join(bdfpath,'*/*/*.fits')) if ftype in f]
        # Loop over all files
        for fpathname in flist:
            fpath, fname = os.path.split(fpathname)
            # Get file date - type 1
            ftime = None
            match = re.search(r'^\d{4}-\d\d-\d\d',fname)
            if match:
                # YYYY-MM-DD type string found -> get date
                ftime = time.strptime(match.group(), '%Y-%m-%d')
            # Get file date - type 2
            if not ftime:
                match = re.search(r'\d{4}[A-Za-z]+\d+',fname)
                if match:
                    # YYYYMmmDD type string found -> get date
                    ftime = time.strptime(match.group(), '%Y%b%d')
            # Get file date - get yesterday's date
            #print('%s - %d-%d-%d' % (fname,ftime.tm_year,ftime.tm_mon,ftime.tm_mday) )
            # Check if file is already there
            found = False
            for f in typefiles:
                if fname in f:
                    found = True
                    log.debug('File %s is already in %s' % (fname, os.path.split(f)[0]) )
                    break
            # if not found: copy it there
            if not found:
                # Get target folder
                tpath = os.path.join(bdfpath, time.strftime('%Y-%m-%d',ftime) )
                if not os.path.exists(tpath): os.makedirs(tpath)
                tpath = os.path.join(tpath, ftype)
                if not os.path.exists(tpath): os.makedirs(tpath)
                # Copy file
                log.debug('File %s copied to %s' % (fname, tpath) )
                shutil.copy(fpathname, os.path.join(tpath, fname) )

    ### Copy Raw Data
    # get last part of source_folder i.e. 2018-02-08_galaxieslab1group2_NGC_2129_7K9
    srest, spart = os.path.split(source_folder)
    ssplit = spart.split('_')
    # Get date, username, object name from folderpart
    #sdate = ftime = time.strptime(ssplit[0], '%Y-%m-%d')
    #suser = ssplit[1][0].upper()+ssplit[1][1:] # i.e. Galaxieslab1group2
    suser = os.path.split(srest)[1] # Get from upper level folder -- i.e. rich
    suser = suser[0].upper() + suser[1:] # Uppercase first character -- i.e. Rich
    sobject = ''.join(ssplit[2:-1]) # i.e. NGC_2129
    # Get list of files in raw folder
    sfiles = glob.glob(os.path.join(source_folder,'raw/science','*.fits') )
    if len(sfiles) < 1:
        log.error("No files found for folder %s" % source_folder)
        continue
    # Get exposure time from first file
    df = DataFits()
    df.loadhead(sfiles[0])
    expt = int(df.getheadval('EXPTIME'))
    # Make raw folder for files
    rpath = os.path.join(outpath, suser)
    if not os.path.exists(rpath): os.makedirs(rpath)
    #opath = "%s_%ds_%s" % (sobject, expt, time.strftime('%y%m%d',sdate) ) # Taken out as we go back to queue names
    opath = spart # new version, just take path from above
    rpath = os.path.join(rpath, opath )
    if not os.path.exists(rpath): os.makedirs(rpath)
    os.system('chmod 775 %s' % rpath)
    # Copy files
    for f in sfiles:
        #df.loadhead(f)
        #expt = int(df.getheadval('EXPTIME'))
        #rname = os.path.split(f)[1].replace('.fits', '_%ds.RAW.fits' % expt )
        # Change output name to _RAW.fits
        rname = os.path.split(f)[1].replace('.fits', '_RAW.fits' )
        # Copy the file
        rname = os.path.join(rpath, rname)
        log.debug('Copy %s to %s' % (os.path.split(f)[1], rname) )
        shutil.copy(f, rname )
        os.system('chmod 664 %s' % rname )
        # Change the Observer name in the output file
        df.load(rname)
        if 'rechelt' in df.getheadval('OBSERVER') and not 'rechelt' in suser:
            df.setheadval('OBSERVER',suser)
            df.save(rname)

    ### Make PipeRun file
    #sdate = time.strftime('%y%m%d', sdate) # change sdate to YYMMDD
    # Make piperun filepathname
    #piperun = 'piperun_%s_%s_%d_%s.txt' % (suser, sobject, expt, sdate )
    piperun = spart + '.txt'
    piperun = os.path.join(piperunpath, piperun)
    # Make file header
    #runame = "%s %s %s" % (suser,sobject,sdate)
    runame = piperun
    text = """# === Piperun file for %s ===

# !!! Auto-generated Pipeconf file - may be overwritten !!!
pythonpath = /data/scripts/DataReduction/source
pipeconf = /data/scripts/DataReduction/pipeconf_stonedge_auto.txt
pipemode = stoneedge
loglevel = DEBUG
logfile = /data/scripts/DataReduction/PipeLineLog.txt
""" % (runame)
    # Add custom log file
    #logfile = '%s_%s_%d_%s_pipelog.txt' % (suser, sobject, expt, sdate )
    logfile = spart + '_pipelog.txt'
    logfile = os.path.join(rpath, logfile)
    text += logfile + '\n'
    # Add input and output folders
    text += 'inputfiles = %s\n' % os.path.join(rpath, '*RAW.fits')
    text += 'outputfolder = %s\n' % rpath
    # Save file
    log.info('Writing piperun at %s' % piperun)
    outf = open(piperun,'wt')
    outf.write(text)
    outf.close()
