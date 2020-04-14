#!/usr/bin/env python
""" PIPE HOTPIX - Version 1.0.0

    Template for StepHotPix in pipeline.

    This module defines the HAWC pipeline step parent object. Pipe steps are
    the modules responsible for all HAWC data reduction. They are called by
    the pipeline and work with pipedata objects. All pipe step objects are
    descendants from this one. Pipe steps are callable objects that return
    the reduced data product (as pipedata object).
    
    @author: berthoud
"""

import os # os library
import numpy # numpy library
import scipy.ndimage #scipy sublibrary
import logging # logging object library
from scipy.ndimage import median_filter #Used to filter hot pixels
from darepype.drp import StepParent # pipe step parent object

class StepHotpix(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version

    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
	# call superclass constructor (calls setup)
        super(StepHotpix,self).__init__()
	# list of data
        self.datalist = [] # used in run() for every new input data file
	# set configuration
        self.log.debug('Init: done')
    
    def setup(self):
        """ ### Names and Prameters need to be Set Here ###
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
        self.name='hotpix'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'hpx'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['hotpixfile', 'search',
            'Filename for clean file or "search" for searching ' +
            'file in cleanfolder (default = search)'])

    def run(self):
        """ Runs the hot pix removal algorithm. The self.datain is run
            through the code, the result is in self.dataout.
            Tolerance is the number of standard deviations used to cutoff
            the hot pixels.
        """
        # Copy input to output data
        self.dataout = self.datain.copy()
        img = self.datain.image
        ''' Cleaning Algorithm '''
        #Apply a filter that creates a threshold for hotpixels
        blurred = median_filter(img, size=2)
        difference = img - blurred
        threshold = 10*numpy.std(difference)
        #Find the hotpixels
        hot_pixels = numpy.nonzero((numpy.abs(difference[1:-1,1:-1])>threshold))
        hot_pixels = numpy.array(hot_pixels) +1 #ignored the edges
        #This is the image with the hot pixels removed
        for y,x in zip(hot_pixels[0],hot_pixels[1]):
            img[y,x]=blurred[y,x]
        ''' Cleaning Algorithm (end) '''
        self.dataout.image = img
        # Set complete flag
        self.dataout.setheadval('COMPLETE',1,
                                'Data Reduction Pipe: Complete Data Flag')

    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.log.debug('Reset: done')
        
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step hotpix')

        # log message
        self.log.info('Testing pipe step hotpix - Done')
    
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
    StepHotpix().execute()

""" === History ===
    2014-06-30 New file created by Neil Stilin from template file by Nicolas Chapman
"""
