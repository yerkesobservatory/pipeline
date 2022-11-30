#!/usr/bin/env python
""" 
        Pipestep HDR (High Dynamic Range)

    This module corrects a pair of raw image low and high gain image files for 
    detector dark, bias, and flat effects. The step crops the overscan from
    input images, reduces high- and low-gain images individually, then
    performs a Gaussian interpolation to remove nans from the data before
    combining into one HDR image and downsampling. Optionally, the step can also
    apply the astroscrappy function to remove cosmic ray defects from the
    downsampled image.
    
    The step requires as input two files, a low gain file and a high gain file.
    It produces one output file.
    
    It uses StepLoadAux functions to access the following auxiliary files:
        - masterpfit (x2): A 3D image containing one HDU with two planes, constructed
            by performing a polynomial fit (currently linear) to a set of 
            short-exposure darks taken at different exposure times (currently between 1s
            and 64s). The first image plane has pixel values which are the slopes of
            linear fits to data from each pixel; the second image contains
            the intercepts and takes the place of a zero-exposure image (bias image)
            for the CMOS camera. Here, the intercept image is used for bias subtraction.
            The step uses two of these, one each for high- and low-gain.
        - masterdark (x2): A matched pair of high- and low-gain master darks of 
            the same exposure length as the raw sky image. MDARK.fits or MMDARK.fits
            files. Each contains one HDU.
        - mastertfit (x2): A matched pair of high- and low-gain master darks of 
            the same exposure length as the raw sky image. TFIT files.
            Each contains one HDU containing a 3D image with two planes constructed
            by performing a polynomial fit (currently linear) to a set of long-exposure
            darks taken at different heatsink temperatures. The first image plane has
            pixel values which are the slopes of linear fits to the data from each pixel.
            The second plane contains the intercepts (this is equivalent to a dark
            image taken at a heatsink temperature of 0 C).    
        - masterflat: A file consisting of three HDUs. The first contains a 3D
            flat-field image with a plane each for the high and low gain data
            (used for flat correction); the second contains a gain ratio image
            (used for creating an HDR output image); the third contains a table of
            statistical information about the flats (not used here).
            
    The high-gain mode of the CMOS detector currently causes low-gain images to be
    loaded into the second HDU, and the header data to be loaded into an XTENSION
    header in the second image position. The step locates the low-gain file and
    swaps the data into the correct positions.
    
    The current version of this code assumes that there is a unique master file in the
    master dark directory with exactly the same exposure time as the file to be reduced.
    That is, it assumes that at each exposure time, there is either an MDARK.fits or
    TFIT.fits file, but not both. It then automatically detects whether the file with
    the matching exposure time is an MDARK.fits or TFIT.fits.
    
    Authors: Carmen Choza, Al Harper, Marc Berthoud
"""

import os # os library
import sys # sys library
import logging # logging object library
import numpy as np
from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans # For masking/replacing
from skimage.measure import block_reduce # for optional downsampling
import astroscrappy # for optional cosmic-ray removal


class StepBdf(StepLoadAux, StepParent):
    """ Pipeline Step Object to calibrate CCD camera files
    """
    
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # Call superclass constructor (calls setup)
        super(StepBdf,self).__init__()
        
        # Bias values
        self.biasloaded = False # indicates if bias has been loaded
        self.bias = None # Numpy array object containing bias values
        self.biasname = '' # name of selected bias file
               
        # Dark values
        self.darkloaded = False # indicates if high-gain dark has been loaded
        self.dark = None # Numpy array object containing high-gain dark values
        self.darkname = '' # name of selected high-gain dark file
        
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
        self.name='bdf'
        # Set procname.
        self.procname = 'BDF'
        
        ## SET UP PARAMETER LIST AND DEFAULT VALUES
        
        # Clear Parameter list.
        self.paramlist = []
        # Append parameters.
        self.paramlist.append(['reload', False,
            'Set to True to look for new pfit/flat/dark files for every input'])
        self.paramlist.append(['remove_cosmics', True, 
        	'Set to False to omit cosmic ray removal'])
        self.paramlist.append(['psffwhm', 5.0, 
        	'Set to psf full width half maximum'])
        self.paramlist.append(['fsmode', 'convolve', 
        	'Set to fsmode for astroscrappy'])
        self.paramlist.append(['cleantype', 'medmask', 
        	'Set to cleantype for astroscrappy'])
            
            
        # Set root names for loading parameters with StepLoadAux.
        self.loadauxsetup('bias')
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
        ### Load bias, dark, and flat files and create numpy arrays with their data
        
        # Set loaded flags to false if reload flag is set
        if self.getarg('reload'):
            self.biasloaded = False
            self.darkloaded = False
            self.flatloaded = False
            
        # Load bias file
        self.log.debug('LoadBias: Start')
        if not self.biasloaded:
            self.biasname = self.loadauxname('bias', multi = False)
            self.bias = DataFits(config = self.config)
            #self.log.info('File loaded: %s' % biasname)
            #if(self.biasname == None ):
                #self.log.error('Bias calibration image not found.')
            	#raise RuntimeError('No bias file loaded')
            self.bias.load(self.biasname)
            self.biasloaded = True
            
        # Load dark file
        if not self.darkloaded:
            self.darkname = self.loadauxname('dark', multi = False)
            self.dark = DataFits(config = self.config)
        	#self.log.info('File loaded: %s' % darkname)
        	#if(self.darkname == None):
                #self.log.error('Dark calibration image not found.')
                #raise RuntimeError('No dark file loaded')
            self.dark.load(self.darkname)
            self.darkloaded = True
            self.log.debug('LoadDark: done')

        # Load flat file
        if not self.flatloaded:
            self.flatname = self.loadauxname('flat', multi = False)
            self.flat = DataFits(config = self.config)
            #if(self.flatname == None):
            	#self.log.error('Flat calibration frame not found.')
            	#raise RuntimeError('No flat file loaded')
            self.flat.load(self.flatname)
            self.flatloaded = True
        	#self.log.debug('LoadFlat: done')
        	
        
        # Create numpy arrays from calibration image data
        bias = self.bias.image
        dark = self.dark.image
        flat = self.flat.image
        
        
        # Set up to handle either TFIT or MDARK files
        dname = os.path.split(self.darkname)
        if '_TFIT.fit' in dname[1]:    # Make temperature-matched darks from TFIT data        
            temp = self.datain.getheadval('secondar')        
            dark = dark[1] + dark[0]*temp
            
            
        # Set up to handle either BTFIT or MBIAS files
        bname = os.path.split(self.biasname)
        if '_BTFIT.fit' in bname[1]:    # Make temperature-matched biases from BTFIT data        
            temp = self.datain.getheadval('secondar')        
            bias = bias[1] + bias[0]*temp  
                                     

        # Get the image to be corrected, convert to float, and find exposure time.        
        image = self.datain.image * 1.0
        image_exp = self.datain.getheadval('EXPTIME')
        
        # Find exposure time of dark image
        dark_exp = self.dark.getheadval('EXPTIME')
                

        # Process data
        dark_current = (dark - bias) / dark_exp
        databdf = ((image - bias) - (dark_current * image_exp)) / flat
        
        
        # Replace NANs with interpolated values
        kernel = Gaussian2DKernel(x_stddev=2) # Make kernel for Gaussian interpolation
        nanmask = np.isnan(databdf)
        self.log.debug(np.sum(nanmask)) 
        databdf = interpolate_replace_nans(databdf, kernel)
        self.log.debug(np.sum(np.isnan(databdf)))
        
        
        # Option to use astroscrappy to remove cosmic ray artifacts
        psffwhm = self.getarg('psffwhm')
        fsmode = self.getarg('fsmode')
        cleantype = self.getarg('cleantype')        
        if self.getarg('remove_cosmics'):
            crarray, outdata = astroscrappy.detect_cosmics(databdf, psffwhm=psffwhm, fsmode=fsmode, cleantype=cleantype)
        else:
            outdata = databdf

        # Make dataout
        self.dataout = self.datain.copy()
        self.dataout.image = outdata
        
        # Add bias, dark, and flat files to History
        self.dataout.setheadval('HISTORY', 'BIAS: %s' % self.biasname)
        self.dataout.setheadval('HISTORY', 'DARK: %s' % self.darkname)
        self.dataout.setheadval('HISTORY', 'FLAT: %s' % self.flatname)
            
        
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.biasloaded = False
        self.bias = None
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
    StepHdr().execute()

'''HISTORY:
2022-1-5 - Set up file, most code copied from StepBiasDarkFlat - Marc Berthoud
'''
