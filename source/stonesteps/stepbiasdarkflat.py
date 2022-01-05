#!/usr/bin/env python
""" 
    Pipestep BiasDarkFlat

    This module defines the pipeline step that corrects raw image files by subtracting
    a bias image, subtracting a dark-charge image, and dividing by a flat image.
    
    It uses StepLoadAux functions to call the following files:
        - masterbias: a median or mean average of zero-exposure dark images
        - masterdark : the difference between a median or mean average of  dark images
            with exposure similar to the exposure of the raw image and a masterbias image
            (the result is a dark-charge image that can then be scaled to match the
            exact exposure of the raw image)
        - masterflat: a median of flat images with exposures that result in data values
            in the middle to high linear gain range of the detector pixels. It is
            assumed in this pipe step that the masterflat image has been normalized
            to its own median or mean.
    
    Authors: Atreyo Pal, Lorenzo Orders
"""
import os # os library
import sys # sys library
import numpy as np # numpy library
import logging # logging object library
from astropy.io import fits #package to recognize FITS files
from darepype.drp import DataFits # pipeline data object class
from darepype.drp import StepParent # pipestep stepparent object class
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class

class StepBiasDarkFlat(StepLoadAux, StepParent):
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
    """
    
    stepver = '1.0' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # Call superclass constructor (calls setup)
        super(StepBiasDarkFlat,self).__init__()

        # Bias values
        self.biasloaded = False # indicates if bias has been loaded
        self.bias = None # Numpy array object containing bias values
        self.biasname = '' # name of selected bias file
        
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
        self.name='biasdarkflat'
        # Set procname.
        self.procname = 'BDF'
        
        ## SET UP PARAMETER LIST AND DEFAULT VALUES
        
        # Clear Parameter list.
        self.paramlist = []
        # Append parameters.
        self.paramlist.append(['reload', False,
            'Set to True to look for new bias files for every input'])
        self.paramlist.append(['intermediate', False,
            'Set to T to include the result of bias, dark, and flat'
            'subtraction'])
        # Set root names for loading parameters with StepLoadAux.
        self.loadauxsetup('bias')
        self.loadauxsetup('dark')
        self.loadauxsetup('flat')
        
        ## SET LOGGER AND FINISH UP
        
        # Set Logger for this pipe step.
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)  
        # Confirm end of setup.
        self.log.debug('Setup: done')
    
    def loadbias(self):
        """ Loads the bias image corresponding to the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        self.log.debug('LoadBias: Start')

        # Search for the closest master bias image
        name = self.loadauxname('bias', multi = False)
        self.log.info('File loaded: %s' % name)
        if(name == None):
            self.log.error('Bias calibration image not found.')
            raise RuntimeError('No bias file loaded')
        
        # Load bias object
        bias = DataFits(name)
        bias.load()
        self.bias = bias.imageget()

        # Finish up
        self.biasloaded = True
        self.biasname = name
        self.log.debug('LoadBias: done')

    def loaddark(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            It is assumed that the dark images are already bias corrected
        """
        self.log.debug('LoadDark: Start')
        
        # Find name of closest master dark image.
        name = self.loadauxname('dark', multi = False)
        self.log.info('File loaded: %s' % name)
        if(name == None):
            self.log.error('Dark calibration image not found.')
            raise RuntimeError('No dark file loaded')
        
        # Load dark data.
        dark = DataFits(name)
        dark.load()
        self.dark = dark.imageget()

        # Get the dark exposure time from the header of the dark object
        self.dark_exp_length = dark.getheadval('EXPTIME')
        
        # Finish up 
        self.darkloaded = True 
        self.darkname = name
        self.log.debug('LoadDark: done')
            
    def loadflat(self):
        """ Loads the flat image information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            It is assumed the flat file is already bias and dark corrected
        """
        self.log.debug('LoadFlat: Start')

        # Find closest master flat image.
        name = self.loadauxname('flat', multi = False)
        if(name == None):
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file loaded')
        
        # Load flat data.
        flat = DataFits(name)
        flat.load()
        self.flat = flat.imageget()

        # Finish up
        self.flatloaded = True  
        self.flatname = name 
        self.log.debug('LoadFlat: done')
    
    def subtract_bias(self, image, bias):
        """
        Subtract master bias from image.
        Arguments:
        - image: numpy.ndarray, a raw image from which bias will be subtracted
        - bias : numpy.ndarray, a master bias image to subtract
        Returns: numpy.ndarray, bias-subtracted image
        """

        # Subtract bias.
        self.log.debug('Subtracting bias...')
        result = image - bias
        self.log.debug('Subtracted bias.')
        return result

    def subtract_dark(self, image, dark, scale=False, img_exposure=None, 
                      dark_exposure=None):
        """
        Subtract dark charge from an image.
        Arguments:
        - image:         numpy.ndarray, raw image from which dark will be subtracted
        - dark:          numpy.ndarray`, dark image with bias already subtracted
        - scale:         bool, if True, scale dark to match exposure of raw image
        - img_exposure:  double, optional, exposure time of raw image
        - dark_exposure: double, optional, exposure time for dark image        
        Returns: numpy.ndarray, dark-subtracted image.
        """
        
        self.log.debug('Subtracting dark...')
        # If scale = True, scale the dark to match exposure time of raw image.
        # Assumes dark current is linear.
        if scale:
            dark_scaled = dark * img_exposure / dark_exposure
            result = image - dark_scaled
        else:
            result = image - dark
        
        # Finish up.
        self.log.debug('Subtracted dark.')
        return result
    
    def flat_correct(self, image, flat):
        """
        Correct the image with a flat-field image.
        Arguments:
        - image: numpy.ndarray, data to be corrected.
        - flat:  numpy.ndarray, flat-field image normalized to its mean or median.
        Returns: numpy.ndarray, flat-corrected image.
        """
        self.log.debug('Correcting flat...')
        
        # Divide by the flat
        flat_corrected = image.copy()
        
        flat_median = np.median(flat)

        # Normalize the flat.
        flat_normed = flat / flat_median

        # divide through the flat
        flat_corrected = image / flat_normed
        self.log.debug('Corrected flat.')
        return flat_corrected

    def run(self):
        """ Runs the correction algorithm. The corrected data is
            returned in self.dataout
        """
        # Load auxiliary files, unless they have already been loaded.
        # Warn the user if any keywords listed in biasfitkeys are different
        # in the current data than in the loaded file. The warning code
        # is only used if multiple files are reduced and reload==False.

        # Load bias file
        if not self.biasloaded or self.getarg('reload'):
            self.loadbias()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.biasfitkeys)):
                if self.biaskeyvalues[keyind] != self.datain.getheadval(self.biasfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.biasfitkeys[keyind])
        # Load dark file
        if not self.darkloaded or self.getarg('reload'):
            self.loaddark()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.darkfitkeys)):
                if self.darkkeyvalues[keyind] != self.datain.getheadval(self.darkfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.darkfitkeys[keyind])
        # Load flat file
        if not self.flatloaded or self.getarg('reload'):
            self.loadflat()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.flatfitkeys)):
                if self.flatkeyvalues[keyind] != self.datain.getheadval(self.flatfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.flatfitkeys[keyind])
        
        ## PREPARE TO CORRECT IMAGE 
        
        # in the config file, set the 'intermediate' variable to either T or F
        # to enable saving of intermediate steps as additional HDUs
        save_intermediate_steps = self.getarg('intermediate')
        
        # Get the image to be corrected, convert to float, and find exposure time.        
        image = self.datain.image * 1.0
        image_exp = self.datain.getheadval('EXPTIME')
        
        # Create self.dataout by copying self.datain. This loads the output object
        # with the config and header of the input file and a placeholder image
        # for the eventual output image.
        self.dataout = self.datain.copy()

        ## CORRECT THE IMAGE

        # Subtract bias from image
        bias_subtracted = self.subtract_bias(image, self.bias)
        
        # Subtract dark from image
        dark_subtracted = self.subtract_dark(bias_subtracted, self.dark, scale=True,
                                   img_exposure=image_exp,
                                   dark_exposure=self.dark_exp_length)

        # Apply flat correction to image
        flat_correct = self.flat_correct(dark_subtracted, self.flat)
        
        # Copy corrected image into self.dataout
        self.dataout.image = flat_correct
        
        # Add bias, dark, and flat files to History
        self.dataout.setheadval('HISTORY', 'BIAS: %s' % self.biasname)
        self.dataout.setheadval('HISTORY', 'DARK: %s' % self.darkname)
        self.dataout.setheadval('HISTORY', 'FLAT: %s' % self.flatname)

        # if saving intermediate steps, save as individual HDUs
        # if not, the final image should be saved as the primary image
        if (save_intermediate_steps):         
            dataname = "DARK_SUBTRACT"
            self.dataout.imageset(dark_subtracted, imagename=dataname)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                    dataname=dataname)
            self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname,
                                    dataname=dataname)

            dataname = "BIAS_SUBTRACT"
            self.dataout.imageset(bias_subtracted, imagename=dataname)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                    dataname=dataname)

        # Set the base filename for the DataFits output object.
        self.dataout.filename = self.datain.filename
            
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
    StepBiasDarkFlat().execute()
    
'''HISTORY:
2020-12-04 - Reintroduce the normalized flat code to maintain continuity with prior versions
             of the BDF step. Reorganized the order of saving the HDUs  - Lorenzo Orders
2020-09-22 - Eliminate unnecesary code and update documentation - Al Harper and Lorenzo Orders
2020-09-19 - Removed code to normalize flat and made other minor edits - Al Harper
2020-07-26 - Saving intermediate BDF steps and removal of CCDProc code - Lorenzo Orders
2018-08-02 - Bias/Dark correction of darks/flats moved to stepmasterbias/dark/flat - Matt Merz
07/28/2017 - Script created by Atreyo Pal
'''
