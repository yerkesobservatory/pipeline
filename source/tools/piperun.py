#! /usr/bin/env python

""" === Pipeline Run ===

    This tool facilitates setting up data reduction using the pipeline.
    The tool requires one (or several) pipeline run description file which
    specify input files, configuration and output location.
    
    Usage:
      piperun.py pipe_run_description_file.txt
    OR
      piperun.py piperunfile1.txt piperunfile2.txt piperunfile3.txt  
    
    Pipe Run Description File:
    This file contains two sections: parameters and input file info.
    The following parameters can be used in the file
        pipeconf = folderpathname for pipeline config file !! REQUIRED !!
        pythonpath = list of paths to append to pythonpath
        pipemode = pipemode to use for reduction
        loglevel = level for logger (DEBUG / INFO (default) / WARN / ERROR)
        logfile = logfile(s) filepathname
        singlefile = T/F if True each file is run through a separate pipeline
                     instead of running all the files through at once (F is default)
        outputfolder = folder to write any results to (default is '.')
    The file can be specified in the following way
        inputfiles = 
        file1.fits
        file2.fits
        file10*.fits # same with wildcards
    Multipe pipeline configuration files and logfiles can be specified
    in the same manner as multiple input files are specified (but no
    "glob" style wildcards are allowed).
    The # character is used to specify comments, all text behind such
    a character is ignored.
    
    Arguments:
    * The only argument is the filename of the pipe run description file(s).
    * If no argument is given piperun looks for a pipe run description
      file called piperun.txt in the current folder.
    * If several pipeline description files are given, piperun will
      start each piperun with itself as a subprocess.
    * This tool has no other arguments, just call pipeline.py directly
      (use 'python pipeline.py -h' for details) to call the pipeline.
"""

### Basic Imports
import os 
import sys
import glob
import logging
import subprocess
import argparse

### Read run description
# Get file from arguments
parser = argparse.ArgumentParser(description = "Pipe Run")
parser.add_argument('prdfile', default='piperun.txt', type=str, nargs = '*',
                    help = 'pipe run description file name (default = piperun.txt)')
args = parser.parse_args()

### Run arguments
if len(args.prdfile) > 1:
    # Set up logger
    loglevel = 'INFO'
    logging.basicConfig(level=loglevel)
    log = logging.getLogger('piperun')
    # Get filepathname of piperun
    progfile = os.path.realpath(__file__)
    # Loop through piperuns
    for prunfile in args.prdfile:
        cmd = 'python %s %s' % (progfile, prunfile)
        log.info('Running Command: %s' % cmd)
        process = subprocess.Popen(cmd, shell=True )
        process.wait()
    # Exit
    exit()

### Process run description (get pythonpath)
# Read run description
rundesc = [l.strip() for l in file(args.prdfile[0])]
# Clear empty lines, strip comments, combine multi-line entries
rundict = {}
l1 = '' # full entry with several lines
rundesc.append('=') # to make sure last entry is added to rundict
for l in rundesc:
    # Search comment location and cut if found
    cloc = l.find('#')
    if cloc > -1: l = l[:cloc]
    # If line is new parameter, append old entry to 
    if '=' in l:
        eloc = l1.find('=')
        if eloc > -1:
            key = l1[:eloc].strip()
            val = l1[eloc+1:].strip()
            rundict[key]=val
        l1 = l
    # append line to previous entry
    else: l1 +='\n'+l

print(rundict)
# Look for pytonpath
if 'pythonpath' in rundict:
    sys.path.append(rundict['pythonpath'])

### Pipeline imports
from drp.pipeline import PipeLine
from drp.dataparent import DataParent

### Go through config options: Logfile Loglevel Pipemode
# Get outputfolder if available
if 'outputfolder' in rundict:
    outputfolder = rundict['outputfolder']
else:
    outputfolder = '.'
# Make outputfolder is required (needed before logging setup)
if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)
# Set Logger: Get loglevel and set up
if 'loglevel' in rundict:
    loglevel = rundict['loglevel']
else:
    loglevel = 'INFO'
logging.basicConfig(level=loglevel)
# Set up logfiles if needed
if 'logfile' in rundict:
    logfiles = rundict['logfile']
    logfiles = logfiles.split('\n')
    logformat = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    for logf in logfiles:
        fhand = logging.FileHandler(logf)
        fhand.setFormatter(logging.Formatter(logformat))
        logging.getLogger().addHandler(fhand)
        print(logf)
# Setup log for piperun
log = logging.getLogger('piperun')
log.info('Setting Up: %s' % os.path.split(args.prdfile[0])[1])
# Get pipeconfig file
if 'pipeconf' in rundict:
    pipeconf = rundict['pipeconf']
else:
    log.error('Missing pipeconf in %s' % args.prdfile)
    raise ValueError('Missing pipeconf in %s' % args.prdfile)
# Load pipeconf, merge if there are multiple pipeconfs
pipeconf = pipeconf.split('\n')
log.debug('Loading Pipeconf=' + repr(pipeconf))
pipeconf = DataParent(config=pipeconf).config

# Get pipemode
if 'pipemode' in rundict:
    pipemode = rundict['pipemode']
else:
    pipemode = None
# Get singlefile
singlefile = False
if 'singlefile' in rundict:
    if rundict['singlefile'][0] in ['T','t']:
        singlefile = True

### Get input file list
# Option 1: Filenames (treat each with glob)
filelist = [] # list of input files to reduce
if 'inputfiles' in rundict:
    inputfiles = rundict['inputfiles'].split('\n')
    for infname in inputfiles:
        # Glob it
        fglob = glob.glob(infname.strip())
        fglob.sort()
        # If none found - Warning
        if len(fglob) == 0:
            log.warn("Glob: No files found for %s" % infname)
        # Add to filelist
        filelist += fglob
        
# Option 2: Folders (below inputfolder) and fileid

# No Files: Issue error and retunt
if len(filelist) < 1:
    log.error('No input files')
    raise ValueError('PipeRun: No input files')

### Copy Files and decompress if needed
filein = []
for fname in filelist:
    log.info("Copy file = " + os.path.split(fname)[1])
    if os.path.exists(fname.strip()):
        # Get new name (remove .bz2 if present)
        nameonly = os.path.split(fname)[1]
        if '.bz2' in fname[-4:]:
            nameonly = nameonly[:-4]
        newname = os.path.join(outputfolder,nameonly)
        # Copy the file
        if not os.path.exists(newname):
            os.system('cp ' + fname + ' ' + outputfolder)
            # Unzip if needed
            if '.bz2' in fname[-4:]:
                log.debug('Unzipping File ' + nameonly)
                os.system('bunzip2 '+ newname+'.bz2')
        # Add log / to filelist
        log.debug('Adding File ' + nameonly)
        filein.append(newname)

### Setup and Run pipeline
if singlefile:
    for f in filein:
        pipe = PipeLine(config=pipeconf)
        pipe(f, pipemode = pipemode, force=True)
else:
    pipe = PipeLine(config=pipeconf)
    pipe(filein, pipemode = pipemode, force=True)

### Final message
log.info('Finished piperun: %s' % os.path.split(args.prdfile[0])[1])
log.info("Reduced files:")
for fname in filein:
    log.info("   " + os.path.split(fname)[1])
