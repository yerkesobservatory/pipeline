#!/usr/bin/env python
""" PIPE STEP LOAD INPUT FILES - Version 1.0.0

    This module provides a pipe step with a simple mechanism to blah blah blah
    
    @author: Matt Merz
    
"""

import os # os library
import glob # glob library
import string # for string.join
import time # time library
import logging
from datetime import datetime
from drp.dataparent import DataParent # Pipeline Data object
from drp.stepniparent import StepNIParent # To check if we have datain or [datain, datain, ...]

class StepLoadInput(StepNIParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """

    def setup(self):
        """
        """
        ### Set names
        # Name of pipeline reduction step
        self.name='loadinput'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname='load'
        # Set Logger for this pipe step
        self.log=logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Set default parameters
        self.paramlist.append(['filelocation', '/*/*/data/images/StoneEdge/0.5meter/2018/%Y-%m-%d/bias/bias*.fits',
            'Filename for input file(s). Can contain * and ? ' +
            'wildcards to match multiple files to be selected using fitkeys '])
        self.paramlist.append(['includeheadvals','OBSERVAT=StoneEdge',
            "List of header keywords to determine if loaded (unused if '', | separated)"])
        self.paramlist.append(['excludeheadvals','',
            "List of header keywords to determine if loaded (unused if '', | separated)"])
        self.paramlist.append(['fileinclude','',
            "List of strings which filename must contain to be loaded (unused if '', | separated)"])
        self.paramlist.append(['fileexclude','MBIAS',
            "List of strings which filename must not contain to be loaded (unused if '' | separated)"])

    def run(self, inpar = '', data = None, multi = False):
        # need os.expandvars like line 110 of steploadaux
        infile = self.getarg('filelocation')
        inglob = (datetime.strftime(datetime.now(), infile))
        indata = glob.glob(inglob)
        self.log.debug('Files found: %s' % indata)
        infilenameinclude = self.getarg('fileinclude').split('|')
        ininclude=[]
        innamelist = []
        if infilenameinclude[0] !='':
            for i in indata:
                split = i.split('/')[len(i.split('/'))-1]
                innamelist.append(split)
            count = 0
            for i in indata:
                for f in infilenameinclude:
                    if f in innamelist[count]:
                        if f not in ininclude:
                            ininclude.append(i)
                count +=1
        else:
            ininclude=indata
        infilenameexclude = self.getarg('fileexclude').split('|')
        inexclude=[]
        exnamelist = []
        for i in indata:
            split = i.split('/')[len(i.split('/'))-1]
            exnamelist.append(split)
        count = 0
        for i in indata:
            for f in infilenameexclude:
                if f in exnamelist[count]:
                    if f not in inexclude:
                        inexclude.append(i)
            count +=1
        self.log.debug('File(s) excluded by filename: %s' % inexclude)
        indatafinal = list(set(ininclude)-set(inexclude))
        headlist = []
        for innam in indatafinal:
            headlist.append(DataParent(config = self.config).loadhead(innam))
        includelist = self.getarg('includeheadvals').split('|')
        keysinclude = []
        inheadinclude = []
        if includelist[0] != '':
            for f in headlist:
                for i in includelist:
                    keysinclude = i.split('=')
                    if str(f.getheadval(keysinclude[0]))==str(keysinclude[1]):
                        if f not in inheadinclude:
                            inheadinclude.append(f)
        else:
            inheadinclude=headlist
        inheadexclude = []
        excludelist = self.getarg('excludeheadvals').split('|')
        if excludelist[0] != '':
            for f in headlist:
                for i in excludelist:
                    keysexclude = i.split('=')
                    if str(f.getheadval(keysexclude[0]))==str(keysexclude[1]):
                        if f not in inheadexclude:
                            inheadexclude.append(f)
        headlistfinal = list(set(inheadinclude)-(set(inheadexclude)))
        finalfiles = []
        for files in headlistfinal:
            finalfiles.append(files.filename)
        self.dataout=[]
        for f in finalfiles:
            self.dataout.append(DataParent(config = self.config).load(f))

    def runend(self,data):
        """ Method to call at the end of pipe the pipe step call
           - Sends final log messages
        """
        # update header (status and history)
        #for d in data:
        #    self.updateheader(d)
        # clear input arguments
        self.arglist = {}
        self.log.info('Finished Reduction: Pipe Step %s' % self.name)

    
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of 
    """
    StepLoadInput().execute()

""" === History ===

"""
