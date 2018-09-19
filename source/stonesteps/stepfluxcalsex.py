#!/usr/bin/env python
""" 
    Pipestep FluxCalSex

    This module defines the pipeline step to flux calibrate data files.
    The pipe step runs sextractor on the data and compares itentified
    sources with values from the StSci guide star catalog.
    
    Requirements: This step requires the source extractor program see
        https://www.astromatic.net/software/sextractor
      for details.
    
    Author: Amanda Pagul / Marc Berthoud
    
    2DO:
    - Init
    - Setup - determine parameters
    - Run:
      - sextractor
      - get database
      - fit
      - add keywords
      - put results in table (optional)
"""
import os # os library
import sys # sys library
import numpy # numpy library
import logging # logging object library
import shutil # library to provide operations on collections of files
from drp.pipedata import PipeData # pipeline data object
from drp.stepparent import StepParent # pipestep stepparent object

class StepFluxCalSex(StepParent):
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
    """
    
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepBiasDarkFlat,self).__init__()
        # bias values
        self.biasloaded = False # indicates if bias has been loaded
        self.bias = None # CCD data object containing arrays with bias values
        self.biasdata = PipeData() # Pipedata object containing the bias file
        # bias file info and header keywords to fit
        self.biasname = '' # name of selected bias file
        self.biasfitkeys = [] # FITS keywords that are present in bias      
        self.biaskeyvalues = [] # values of FITS keywords (from data file)  
        # dark values
        self.darkloaded = False # indicates if dark has been loaded
        self.dark = None # CCD data object containing arrays with dark values
        self.darkdata = PipeData() # Pipedata object containing the dark file
        # dark file info and header keywords to fit
        self.darkname = '' # name of selected dark file
        self.darkfitkeys = [] # FITS keywords that have to fit for dark     
        self.darkkeyvalues = [] # values of FITS keywords (from data file)  
        # flat values
        self.flatloaded = False # indicates if flat has been loaded
        self.flat = None # CCD data object containing arrays with flat values
        self.flatdata = PipeData() # Pipedata object containing the flat file
        # flat file info and header keywords to fit
        self.flatname = '' # name of selected flat file
        self.flatfitkeys = [] # FITS keywords that have to fit for flat
        self.flatkeyvalues = [] # values of flat keywords (from data file)
        # set configuration
        self.log.debug('Init: done')
        
    def setup(self):
        """ ### Names and Parameters need to be Set Here ###
            Sets the internal names for the function and for saved files.
            Defines the input parameters for the current pipe step.
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
        self.name='biasdarkflat'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'bdf'
        # Set Logger for this pipe step
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['reload', False,
            'Set to True to look for new bias files for every input'])
        # Get parameters for StepLoadAux, replace auxfile with biasfile
        self.loadauxsetup('bias')
        # Get parameters for StepLoadAux, replace auxfile with darkfile
        self.loadauxsetup('dark')
        # Get parameters for StepLoadAux, replace auxfile with flatfile
        self.loadauxsetup('flat')
        # confirm end of setup
        self.log.debug('Setup: done')
    
    '''# Looking for similar exptime
    def closestExp(self):
        input_exptime = self.datain.getheadval('EXPTIME')
        dark_exptime = self.loadauxname('dark', multi = True).getheadval('EXPTIME')
        nearexp = {abs(dark_ave_exptime - exp): exp for exp in dark_exptime} 
        return nearexp[min(nearexp.keys())]
    '''
   
    def run(self):
        """ Runs the calibrating algorithm. The calibrated data is
            returned in self.dataout
        """
        ### Preparation
        self.dataout = self.datain

if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
        python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
        --config=ConfigFilePathName.txt : name of the configuration file
        --test : runs the functionality test i.e. pipestep.test()
        --loglevel=LEVEL : configures the logging output for a particular level
    """
    StepFluxCalSex().execute()
    
'''HISTORY:
2018-09-019 - Started based on Amanda's code. - Marc Berthoud
'''