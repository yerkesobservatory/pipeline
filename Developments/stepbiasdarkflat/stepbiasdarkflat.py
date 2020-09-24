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
import numpy # numpy library
import logging # logging object library
from astropy.io import fits #package to recognize FITS files
from darepype.drp import DataFits # pipeline data object class
from darepype.drp import StepParent # pipestep stepparent object class
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class

class StepBiasDarkFlat(StepLoadAux, StepParent):
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
    """
    
    stepver = '0.2' # pipe step version
    
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
        ### Set Names
        # Name of the pipeline reduction step
        self.name='biasdarkflat'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names. The output file will be file_BDF.fits
        self.procname = 'BDF'
        # Set Logger for this pipe step
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)
        
        # Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['reload', False,
            'Set to True to look for new bias files for every input'])
        self.paramlist.append(['intermediate', False,
            'Set to T to include the result of bias, dark, and flat'
            'subtraction'])
        # Get parameters for StepLoadAux
        self.loadauxsetup('bias')
        # Get parameters for StepLoadAux
        self.loadauxsetup('dark')
        # Get parameters for StepLoadAux
        self.loadauxsetup('flat')
        # confirm end of setup
        self.log.debug('Setup: done')
    
    def loadbias(self):
        """ Loads the bias image corresponding to the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        # Search for the closest master bias image
        self.log.debug('LoadBias: Start')
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
        
        # Load dark DataFits object.
        dark = DataFits(name)

        # load the dark image data
        dark.load()
        self.dark = dark.imageget()

        # Get the dark exposure time from the header of the DataFits object
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
        #create master flat frame
        self.log.debug('LoadFlat: Start')
        name = self.loadauxname('flat', multi = False)
        if(name == None):
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file loaded')
        
        flat = DataFits(name)
        flat.load()
        self.flat = flat.imageget()

        # Finish up
        self.flatloaded = True  
        self.flatname = name 
        self.log.debug('LoadFlat: done')
    
    # The subtract_bias, subtract_dark, and subtract_flat functions are derived
    # from https://github.com/astropy/ccdproc/blob/master/ccdproc/core.py
    def subtract_bias(self, image, bias):
        """
        Subtract master bias from image.
        
        Parameters
        ----------
        - image: `~numpy.ndarray`, raw image from which bias will be subtracted
        - bias : `~astropy.ndarray`, master bias image to subtract
        
        Returns
        -------
        result : `~numpy.ndarray`, bias-subtracted iamge
        """
        self.log.debug('Subtracting bias...')
        result = image - bias
        self.log.debug('Subtracted bias.')
        return result

    def subtract_dark(self, image, dark, scale=False, img_exposure=None, 
                      dark_exposure=None):
        """
        Subtract dark charge from an image.
        
        Parameters
        ----------
        - image : `~numpy.ndarray`
            Raw image from which the dark will be subtracted.
        - dark : `~numpy.ndarray`
            Dark image with bias already subtracted.
        - scale: `bool`, optional.
            If True, scale the dark image by the exposure time to match the raw
            image's exposure. Default is ``False``.
        - img_exposure : `double`, optional.
            Exposure time of the raw image. Only used if scale = True.
            Default is ``None``.
        - dark_exposure : `double`, optional.
            Exposure time for the dark image. Only used if scale = True.
            Default is ``None``.
        
        Returns
        -------
        result : `~numpy.ndarray`
            Dark-subtracted image.
        """
        
        self.log.debug('Subtracting dark...')
        # if dark current is linear, then this first step scales the
        # provided dark to match the exposure time of the raw image
        if scale:
            dark_scaled = dark * img_exposure / dark_exposure
            result = image - dark_scaled
        else:
            result = image - dark
        
        self.log.debug('Subtracted dark.')
        return result
    
    # This code is also from ccdproc
    def flat_correct(self, image, flat):
        """Correct the image with a flat field.
        
        Parameters
        ----------
        - image : `~numpy.ndarray`
            Data to be corrected.
        - flat : `~numpy.ndarray`
            flat-field image
        
        Returns
        -------
        flat_corrected : `~numpy.ndarray`.
            Flat corrected image
        """
        self.log.debug('Correcting flat...')
        
        # Divide by the flat
        flat_corrected = image.copy()
        
        flat_mean = flat.mean()

        # Normalize the flat.
        flat_normed = flat / flat_mean

        # divide through the flat
        flat_corrected = image / flat_normed
        self.log.debug('Corrected flat.')
        return flat_corrected

    def run(self):
        """ Runs the correction algorithm. The corrected data is
            returned in self.dataout
        """
        
        ### Preparation

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
        # in the config file, set the 'intermediate' variable to either T or F
        # to enable saving of intermediate steps as additional HDUs
        save_intermediate_steps = self.getarg('intermediate')

        # Set self.dataout with configuration object of self.datain.
        self.dataout = DataFits(config=self.datain.config)
        
        image = self.datain.image * 1.0
        header = self.datain.header
        image_exp = self.datain.getheadval('EXPTIME')
        
        ## Correct the image

        # Subtract bias from image
        bias_subtracted = self.subtract_bias(image, self.bias)
        
        # Subtract dark from image
        dark_subtracted = self.subtract_dark(bias_subtracted, self.dark, scale=True,
                                   img_exposure=image_exp,
                                   dark_exposure=self.dark_exp_length)

        # Apply flat correction to image
        flat_correct = self.flat_correct(dark_subtracted, self.flat)
        
        # if saving intermediate steps, save as individual HDUs
        # if not, the final image should be saved as the primary image
        if (save_intermediate_steps):
            dataname = "FLAT_CORRECT"
            self.dataout.imageset(flat_correct, imagename=dataname)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                dataname=dataname)
            self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname,
                                dataname=dataname)
            self.dataout.setheadval('HISTORY','FLAT: %s' % self.flatname,
                                dataname=dataname)
            
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
            

        else:
            # Copy corrected image into self.dataout
            self.dataout.image = flat_correct
            self.dataout.header = self.datain.header
            
            # Add bias, dark, and flat files to History
            self.dataout.setheadval('HISTORY', 'BIAS: %s' % self.biasname)
            self.dataout.setheadval('HISTORY', 'DARK: %s' % self.darkname)
            self.dataout.setheadval('HISTORY', 'FLAT: %s' % self.flatname)

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
2020-09-22 - Eliminate unnecesary code and update documentation - Al Harper and Lorenzo Orders
2020-09-19 - Removed code to normalize flat and made other minor edits - Al Harper
2020-07-26 - Saving intermediate BDF steps and removal of CCDProc code - Lorenzo Orders
2018-08-02 - Bias/Dark correction of darks/flats moved to stepmasterbias/dark/flat - Matt Merz
07/28/2017 - Script created by Atreyo Pal
'''