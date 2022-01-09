#!/usr/bin/env python
""" 
    Pipestep Flat HDR (High Dynamic Range)

    This module defines the pipeline step that corrects a pair raw image low
    and high gain images files for detector dark and flat effects.
    
    The step requires as input two files, a low gain file and a high gain file.
    It produces one output file.
    
    It uses StepLoadAux functions to call the following files:
        - masterpfit: ADD HERE WHAT THIS FILE IS
        - masterdark: ADD HERE WHAT THIS FILE IS
        - masterflat: ADD HERE WHAT THIS FILE IS
    
    Authors: Carmen Choza, Al Harper, Marc Berthoud
"""

from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class

class StepFlatHdr(StepLoadAux, StepMIParent):
    """ Pipeline Step Object to calibrate Flatfield High Dynamic Range files
    """
    
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # Call superclass constructor (calls setup)
        super(StepFlatHdr,self).__init__()

        # Pfit values
        self.pfitloaded = False # indicates if bias has been loaded
        self.pfit = None # Numpy array object containing bias values
        self.pfitname = '' # name of selected bias file
        
        # Dark values
        self.darkloaded = False # indicates if dark has been loaded
        self.dark = None # Numpy array object containing dark values
        self.darkname = '' # name of selected dark file
        
        # Flat values
        self.flatloaded = False # indicates if flat has been loaded
        self.flat = None # Numpy array object containing flat values
        self.flatname = '' # name of selected flat file
        
        # Finish up.
        self.log.debug('Init: done')
 
    def setup(self):
        """ ### Internal Names and Parameters need to be set Here ###
            
            Sets the internal names for the function and for saved files.
            - name: The name for the parameter. This name is used when
                    calling the pipe step from the command line or python shell.
                    It is also used to identify the step in the configuration file.
                    The name should be lower-case.
            - procname: Used to identify the step in saved output files.
                    The procname should be upper case.
                        
            Defines the input parameters for the current pipe step.
            The parameters are stored in a list containing the following
            information:
            - parameter name
            - default: A default value for the parameter. If nothing, set
                       to '', for strings, to 0 for integers and to 0.0 for floats.
            - help: A short description of the parameter.
        """
        ### SET NAMES

        # Set internal name of the pipeline reduction step.
        self.name='flathdr'
        # Set procname.
        self.procname = 'FLT'
        
        ## SET UP PARAMETER LIST AND DEFAULT VALUES
        
        # Clear Parameter list.
        self.paramlist = []
        # Append parameters.
        self.paramlist.append(['reload', False,
            'Set to True to look for new pfit/flat/dark files for every input'])
        self.paramlist.append(['intermediate', False,
            'Set to T to include the result of the step'])
        # Set root names for loading parameters with StepLoadAux.
        self.loadauxsetup('pfit')
        self.loadauxsetup('dark')
        self.loadauxsetup('flat')
        
        ## SET LOGGER AND FINISH UP
        
        # Set Logger for this pipe step.
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)  
        # Confirm end of setup.
        self.log.debug('Setup: done')
        
    def run(self):
        """ Runs the correction algorithm. The corrected data is
            returned in self.dataout
        """
        ### Load pfit, dark and flat files
        # Set loaded flags to false if reload flag is set
        if self.getarg('reload'):
            self.pfitloaded = False
            self.darkloaded = False
            self.flatloaded = False
        # Load pfit file
        if not self.pfitloaded:
            self.pfitname = self.loadauxfile('pfit', multi = False)
            self.pfit = DataFits(config = self.config)
            self.pfit.load(self.pfitname)
            self.pfitloaded = True
        # Load mdark file
        if not self.darkloaded:
            self.darkname = self.loadauxfile('dark', multi = False)
            self.dark = DataFits(config = self.config)
            self.dark.load(self.darkname)
            self.darkloaded = True
        # Load mflat file
        if not self.flatloaded:
            self.flatname = self.loadauxfile('flat', multi = False)
            self.flat = DataFits(config = self.config)
            self.flat.load(self.flatname)
            self.flatloaded = True
        
        ### THE IMAGES ARE in DataFits objects
        #   IT IS UNCLEAR Which is high and which is low gain image.
        firstimage = self.datain[0]
        secondimage = self.datain[1]
        
        # Make dataout
        self.dataout = self.datain[0].copy() # could also be new DataFits() or copy of datain[1]]
        
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.pfitloaded = False
        self.pfit = None
        self.darkloaded = False
        self.dark = None
        self.flatloaded = False
        self.flat = None
        self.log.debug('Reset: done')

if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
        python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
        --config=ConfigFilePathName.txt : name of the configuration file
        --test : runs the functionality test i.e. pipestep.test()
        --loglevel=LEVEL : configures the logging output for a particular level
    """
    StepFlatHdr().execute()

'''HISTORY:
2022-1-5 - Set up file, most code copied from StepBiasDarkFlat - Marc Berthoud
'''