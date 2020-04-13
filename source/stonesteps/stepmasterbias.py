#!/usr/bin/env python
""" PIPE STEP MASTER BIAS - Version 1.0.0

    This module creates master bias frames from inputs.
    Should be run with like data inputs (same binning, etc).

    @author: Matt Merz
"""
import os # os library
import sys # sys library
import numpy # numpy library
import logging # logging object library
import astropy
import ccdproc # package for reducing optical CCD telescope data 
import matplotlib.pyplot as plt
from astropy.io import fits #package to recognize FITS files
from darepype.drp import StepMIParent
from darepype.drp import DataFits

class StepMasterBias(StepMIParent):
    """ Stone Edge Pipeline Step Master Bias Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version
    
    def setup(self):
        """ ### Names and Parameters need to be Set Here ###
            Sets the internal names for the function and for saved files.
            Defines the input parameters for the current pipe step.
            Setup() is called at the end of __init__
            The parameters are stored in a list containing the following
            information:
            - name: The name for the parameter. This name is used when
                    calling the pipe step from command line or python shell.
                    It is also used to identify the parameter in the pipeline
                    configuration file.
            - default: A default value for the parameter. If nothing, set
                       '' for strings, 0 for integers and 0.0 for floats
            - help: A short description of the parameter.
        """
        ### Set Names
        # Name of the pipeline reduction step
        self.name='masterbias'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'mbias'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters !!!! WHAT PARAMETERS ARE NEEDED ????? !!!!!
        self.paramlist.append(['combinemethod','median',
                               'Specifies how the files should be combined - options are median, average, sum'])
        self.paramlist.append(['outputfolder','',
                               'Output directory location - default is the folder of the input files'])

    def run(self):
        """ Runs the combining algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        filelist=[]
        for fin in self.datain:
            self.log.debug("Input filename = %s" % fin.filename)
            filelist.append(fin.filename)
        # Make a dummy dataout
        self.dataout = DataFits(config = self.config)
        if len(self.datain) == 0:
            self.log.error('Bias calibration frame not found.')
            raise RuntimeError('No bias file(s) loaded')
        # self.log.debug('Creating master bias frame...')
        # if there is just one, use it as biasfile or else combine all to make a master bias
        if (len(filelist) == 1):
            self.bias = ccdproc.CCDData.read(filelist[0], unit='adu', relax=True)
        else:
            self.bias = ccdproc.combine(filelist, method=self.getarg('combinemethod'), unit='adu', add_keyword=True)
        # set output header, put image into output
        self.dataout.header=self.datain[0].header
        self.dataout.imageset(self.bias)
        # rename output filename
        outputfolder = self.getarg('outputfolder')
        if outputfolder != '':
            outputfolder = os.path.expandvars(outputfolder)
            self.dataout.filename = os.path.join(outputfolder, os.path.split(filelist[0])[1])
        else:
            self.dataout.filename = filelist[0]
        # Add history
        self.dataout.setheadval('HISTORY','MasterBias: %d files used' % len(filelist))

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
    StepMasterBias().execute()
        
        
""" === History ===
    2018-07-23 New step created based on StepRGB - Matt Merz
    2018-08-02 Updates to documentation, step functionality - Matt Merz
"""
