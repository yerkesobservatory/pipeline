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
            'wildcards to match multiple files to be selected using fitkeys ' +
            '(default = %sfolder/*.fits)' % 'bias'])
        self.paramlist.append(['includeheadvals','OBSERVAT',
            'List of header keywords to determine if loaded (unused if '', | separated)'])
        self.paramlist.append(['excludeheadvals','',
            'List of header keywords to determine if loaded (unused if '', | separated)'])
        self.paramlist.append(['fileinclude','',
            'List of strings which filename must contain to be loaded (unused if '', | separated)'])
        self.paramlist.append(['fileexclude','MBIAS',
            'List of strings which filename must not contain to be loaded (unused if '' | separated)'])

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
        inheadlist = []
        for innam in indatafinal:
            inheadlist.append(DataParent(config = self.config).loadhead(innam))
        inkeys = self.getarg('includekeys').split('|')
        inkeysmatch = self.getarg('includevalues').split('|')
        inheadinclude=[]
        if inkeys[0] != '' and inkeysmatch[0] != '':
            for i in inheadlist:
                count = 0
                for f in inkeys:
                    if str(i.getheadval(inkeys[count]))==str(inkeysmatch[count]):
                        if i not in inheadinclude:
                            inheadinclude.append(i)
                    count+=1
        else:
            inheadinclude=inheadlist
        inheadexclude = []
        exkeys = self.getarg('excludekeys').split('|')
        exkeysmatch = self.getarg('excludevalues').split('|')
        if exkeys[0] != '' and exkeysmatch[0] != '':
            for i in inheadinclude:
                count = 0
                for f in exkeys:
                    if str(i.getheadval(exkeys[count]))==str(exkeysmatch[count]):
                        if i not in inheadexclude:
                            inheadexclude.append(i)
                    count+=1
        inheadlistfinal = list(set(inheadinclude)-(set(inheadexclude)))
        finalfiles = []
        for files in inheadlistfinal:
            finalfiles.append(files.filename)
        self.dataout=[]
        for f in finalfiles:
            self.dataout.append(DataParent(config = self.config).loadhead(f))
    
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
