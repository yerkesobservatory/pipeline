#!/usr/bin/env python
""" 
    Pipestep HDR (High Dynamic Range)

    This module corrects a pair raw image low and high gain images files for 
    detector dark, bias, and flat effects. The step crops the overscan from
    input images, reduces high- and low-gain images individually, then
    performs a Gaussian interpolation to remove nans from the data before
    combining into one HDR image and downsampling.
    
    The step requires as input two files, a low gain file and a high gain file.
    It produces one output file.
    
    It uses StepLoadAux functions to call the following files:
        - masterpfit (x2): A 3D image containing one HDU with two planes, constructed
            by performing a polynomial fit (currently linear) to a set of 
            short-exposure darks. The first image has pixel values which are the 
            slopes of linear fits to data from each pixel; the second image contains
            the intercepts and takes the place of a zero-exposure for the CMOS camera. 
            Here, the intercept image is used for bias subtraction. The step uses two
            of these, one for each high- and low-gain.
        - masterdark (x2): A matched pair of high- and low-gain master darks of 
            the same exposure length as the raw sky image. MDARK or MMDARK files.
            Each contains one HDU.
        - masterflat: A file consisting of three HDUs. The first contains a 3D
            flat-field image with a plane each for the high and low gain data
            (used for flat correction); the second contains a gain ratio image
            (used for creating an HDR output); the third contains a table of
            statistical information about the flats (not used here).
            
    The high-gain mode of the CMOS detector currently causes low-gain images to be
    loaded into the second HDU, and the header data to be loaded into an XTENSION
    header in the second image position. The step locates the low-gain file and
    swaps the data into the correct positions.
    
    Authors: Carmen Choza, Al Harper, Marc Berthoud
"""

from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans # For masking/replacing
import scipy.ndimage as nd
import numpy as np
import logging

class StepHdr(StepLoadAux, StepMIParent):
    """ Pipeline Step Object to calibrate Flatfield High Dynamic Range files
    """
    
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # Call superclass constructor (calls setup)
        super(StepHdr,self).__init__()

        # Pfit values
        self.hpfitloaded = False # indicates if bias has been loaded
        self.hpfit = None # Numpy array object containing bias values
        self.hpfitname = '' # name of selected bias file
        
        self.lpfitloaded = False # indicates if bias has been loaded
        self.lpfit = None # Numpy array object containing bias values
        self.lpfitname = '' # name of selected bias file
        
        # Dark values
        self.hdarkloaded = False # indicates if high-gain dark has been loaded
        self.hdark = None # Numpy array object containing high-gain dark values
        self.hdarkname = '' # name of selected high-gain dark file
        
        self.ldarkloaded = False # indicates if low-gain dark has been loaded
        self.ldark = None # Numpy array object containing low-gain dark values
        self.ldarkname = '' # name of selected low-gain dark file
        
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
        self.name='hdr'
        # Set procname.
        self.procname = 'HDR'
        
        ## SET UP PARAMETER LIST AND DEFAULT VALUES
        
        # Clear Parameter list.
        self.paramlist = []
        # Append parameters.
        self.paramlist.append(['reload', False,
            'Set to True to look for new pfit/flat/dark files for every input'])
        self.paramlist.append(['intermediate', False,
            'Set to T to include the result of the step'])
        self.paramlist.append(['splice_thresh', 3000.0,
            'Change to alter the cutoff threshold for combining high- and low-gain images'])
        self.paramlist.append(['divkeys', [], 
            'list of header keywords to divide by sample factor'])
        self.paramlist.append(['multkeys', [], 
            'list of header keywords to multiply by sample factor'])
       
        # Set root names for loading parameters with StepLoadAux.
        self.loadauxsetup('lpfit')
        self.loadauxsetup('hpfit')
        self.loadauxsetup('ldark') 
        self.loadauxsetup('hdark')
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
            self.lpfitloaded = False
            self.hpfitloaded = False
            self.hdarkloaded = False
            self.ldarkloaded = False
            self.flatloaded = False
        # Load pfit file
        if not self.hpfitloaded:
            self.hpfitname = self.loadauxname('hpfit', multi = False)
            self.hpfit = DataFits(config = self.config)
            self.hpfit.load(self.hpfitname)
            self.hpfitloaded = True
        if not self.lpfitloaded:
            self.lpfitname = self.loadauxname('lpfit', multi = False)
            self.lpfit = DataFits(config = self.config)
            self.lpfit.load(self.lpfitname)
            self.lpfitloaded = True
        # Load mdark file
        if not self.hdarkloaded:
            self.hdarkname = self.loadauxname('hdark', multi = False)
            self.hdark = DataFits(config = self.config)
            self.hdark.load(self.hdarkname)
            self.hdarkloaded = True
        if not self.ldarkloaded:
            self.ldarkname = self.loadauxname('ldark', multi = False)
            self.ldark = DataFits(config = self.config)
            self.ldark.load(self.ldarkname)
            self.ldarkloaded = True
        # Load mflat file
        if not self.flatloaded:
            self.flatname = self.loadauxname('flat', multi = False)
            self.flat = DataFits(config = self.config)
            self.flat.load(self.flatname)
            self.flatloaded = True
        
        # Get all relevant image data
        hbias = self.hpfit.image[1]   # high-gain bias from polyfit
        lbias = self.lpfit.image[1]   # low-gain bias from polyfit
        
        hdark = self.hdark.image      # high-gain dark
        ldark = self.ldark.image      # low-gain dark
        
        gain = self.flat.imageget('gain ratio')    # high-gain flat divided by low-gain flat, ratio between modes
        hflat = self.flat.image[1]                 # high-gain flat
        lflat = self.flat.image[0]                 # low-gain flat
        
        
        ### The images are now in DataFits objects
        
        # Get the filename to determine gain
        filename1 = self.datain[0].filenamebegin
        filename2 = self.datain[1].filenamebegin
        
        self.log.debug(filename1)
        self.log.debug(filename2)

    #   if not ('bin1H' in filename1 or 'bin1H' in filename2)
        if 'bin1H' not in filename1 and 'bin1H' not in filename2:
            self.log.error('High-gain file not found.')
            raise RuntimeError('No high-gain file loaded')
        elif 'bin1L' not in filename1 and 'bin1L' not in filename2:
            self.log.error('Low-gain image not found.')
            raise RuntimeError('No low-gain file loaded')
        elif 'bin1L' in filename1:
            hdata_df = self.datain[1]       # Set high-gain data
            ldata_df = DataFits(config = self.config)
            #ldata_df = df.imageset(self.datain[0].imageget(self.datain[0].imgnames[1]))
            dL = self.datain[0]                                  
            ldata_df.imageset(dL.imageget(dL.imgnames[1]))            # Swap data from second HDU to first
        elif 'bin1H' in filename1:
            hdata_df = self.datain[0]       # Set high-gain data
            ldata_df = DataFits(config = self.config)
            # ldata_df = ldata_df.imageset(self.datain[1].imageget(self.datain[1].imgnames[1]))
            dL = self.datain[1]
            ldata_df = DataFits(config = self.config)
            ldata_df.imageset(dL.imageget(dL.imgnames[1]))
        
        
                
        # dataL_df now contains the low-gain file, dataH_df now contains the high-gain file:
        hdata = hdata_df.image[:,:4096] * 1.0       # Crop overscan and convert to float
        ldata = ldata_df.image[:,:4096] * 1.0       
        
        self.log.debug('Shape of Hdata: %s' % repr(np.shape(hdata)))
        self.log.debug('Shape of Ldata: %s' % repr(np.shape(ldata)))
        
        hdata_df.setheadval('CRPIX1', len(hdata)/2)  # Update due to overscan crop
        
        kernel = Gaussian2DKernel(x_stddev=2) # Make kernel for Gaussian interpolation
        
        # Process high-gain data
        hdatabdf = ((hdata - hbias) - (hdark - hbias))/hflat
        
        nanmask = np.isnan(hdatabdf)
        self.log.debug(np.sum(nanmask))
        
        hdatabdf = interpolate_replace_nans(hdatabdf, kernel)
        
        self.log.debug(np.sum(np.isnan(hdatabdf)))
        
        '''Process low-gain data'''
        
        ldatabdf = (((ldata - lbias) - (ldark - lbias))/lflat) * gain
        
        nanmask = np.isnan(ldatabdf)
        self.log.debug(np.sum(nanmask))
        
        ldatabdf = interpolate_replace_nans(ldatabdf, kernel)
        
        self.log.debug(np.sum(np.isnan(ldatabdf)))
        
        # Combine high- and low-gain data into HDR image
        splice_thresh = self.getarg('splice_thresh') # Get crossover threshold
        lupper = np.where(ldatabdf > splice_thresh) # Choose all pixels in low-gain data above a certain threshold parameter
        ldata = ldatabdf.copy()
        HDRdata = hdatabdf.copy()
        HDRdata[lupper] = ldata[lupper]      # Replace upper range of high-gain image with low-gain * gain values
        
        # Downsample image by factor of two
        outdata = nd.zoom(HDRdata,0.5)
        
        # Make dataout
        self.dataout = hdata_df.copy() # could also be new DataFits() or copy of datain[1]]
        self.dataout.image = outdata
            
        # Edit output header keys according to downsample 
        for key in self.getarg('divkeys'):                                # Divide keywords by sample factor
            if not key in self.dataout.header: continue
            self.dataout.setheadval(key,self.dataout.getheadval(key)/2)
        for key in self.getarg('multkeys'):                               # Multiply keywords by sample factor
            if not key in self.dataout.header: continue
                if key == 'GAIN':
                    self.dataout.setheadval(key,self.dataout.getheadval(key)*4)
                else:
                    self.dataout.setheadval(key,self.dataout.getheadval(key)*2)
        
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
    StepHdr().execute()

'''HISTORY:
2022-1-5 - Set up file, most code copied from StepBiasDarkFlat - Marc Berthoud
'''