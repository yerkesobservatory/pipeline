#!/usr/bin/env python
""" PipeStep OverCut

    This pipe step cuts the overscan lines from gsense detectors.
    It also pusts statistical information about the oversampled
    lines into the FITS header.
    
    Statistical information: OVSCMEAN - mean of the overscan data
    
    @author: Marc Berthoud
    
    8  2048
    16 4096

"""

import logging # logging object library
import numpy as np # numpy library
from darepype.drp.stepparent import StepParent

class StepOverCut(StepParent):
    """ DarePype Step ReSample Object
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
        # Name of the pipeline reduction step (same case as config file)
        self.name='overcut'
        # Identifier for saved file names.
        self.procname = 'OVC'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['ovrscncolsfrac', 257, 'Fraction of the detector colums which are overscaned'])
        #self.paramlist.append(['subframe', ':,:', 'subframe to keep in numpy ( i:j, k:l ) notation. Default is :,: i.e. full image'])
        #self.paramlist.append(['keysmulti','', '| separated list of fits keywords to multiply with downsample'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        self.log.debug('Running step %s' % self.name)
        # Set dataout to datain
        self.dataout = self.datain
        # Get statistics from the overscan
        oscfrac = self.getarg('ovrscncolsfrac')
        overscancoln = int(self.datain.image.shape[1]/oscfrac)
        oscnmed = np.median(self.datain.image[:,-overscancoln:])
        self.log.debug('Overscan has %d columns med = %f' % (overscancoln, oscnmed))
        self.dataout.setheadval('OVSCNMED',oscnmed,'Detector overscan columns median value')
        # Remove the overscan
        self.dataout.image = self.datain.image[:,:-overscancoln]
        self.dataout.setheadval('NAXIS1',self.dataout.image.shape[1])
        
            
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of parameters
    """
    StepOverCut().execute()

""" === History ===
    2021-4: Wrote first version - Marc Berthoud
"""
