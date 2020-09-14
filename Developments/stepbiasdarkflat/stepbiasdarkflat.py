#!/usr/bin/env python
""" 
    Pipestep BiasDarkFlat
    This module defines the pipeline step to calibrate data files w.r.t. bias, 
    dark and flat frames. This pipestep works with pipedata object lists called 
    by StepLoadAux and removes bad pixels from the final RGB image. This is very 
    important to smoothen the background texture of the image and the final image
    looks far more clearer once calibration is done.
    
    Authors: Atreyo Pal, Lorenzo Orders
"""
import os # os library
import sys # sys library
import numpy # numpy library
import logging # logging object library
import shutil # library to provide operations on collections of files
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits #package to recognize FITS files
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepParent # pipestep stepparent object
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object

class StepBiasDarkFlat(StepLoadAux, StepParent):
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
    """
    
    stepver = '0.2' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepBiasDarkFlat,self).__init__()
        # bias values
        self.biasloaded = False # indicates if bias has been loaded
        self.bias = None # CCD data object containing arrays with bias values
        self.biasdata = DataFits() # Pipedata object containing the bias file
        # bias file info and header keywords to fit
        self.biasname = '' # name of selected bias file
        self.biasfitkeys = [] # FITS keywords that are present in bias      
        self.biaskeyvalues = [] # values of FITS keywords (from data file)  
        # dark values
        self.darkloaded = False # indicates if dark has been loaded
        self.dark = None # CCD data object containing arrays with dark values
        self.darkdata = DataFits() # Pipedata object containing the dark file
        # dark file info and header keywords to fit
        self.darkname = '' # name of selected dark file
        self.darkfitkeys = [] # FITS keywords that have to fit for dark     
        self.darkkeyvalues = [] # values of FITS keywords (from data file)  
        # flat values
        self.flatloaded = False # indicates if flat has been loaded
        self.flat = None # CCD data object containing arrays with flat values
        self.flatdata = DataFits() # Pipedata object containing the flat file
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
        self.paramlist.append(['intermediate', False,
            'Set to True to include the result of bias, dark, and flat'
            'subtraction'])
        # Get parameters for StepLoadAux, replace auxfile with biasfile
        self.loadauxsetup('bias')
        # Get parameters for StepLoadAux, replace auxfile with darkfile
        self.loadauxsetup('dark')
        # Get parameters for StepLoadAux, replace auxfile with flatfile
        self.loadauxsetup('flat')
        # confirm end of setup
        self.log.debug('Setup: done')
    
   
    def run(self):
        """ Runs the calibrating algorithm. The calibrated data is
            returned in self.dataout
        """
        
        ### Preparation
        # Load bias files if necessary
        if not self.biasloaded or self.getarg('reload'):
            self.loadbias()
        # Else: check data for correct instrument configuration
        # - currently not in use(need improvement)
        else:
            for keyind in range(len(self.biasfitkeys)):
                if self.biaskeyvalues[keyind] != self.datain.getheadval(self.biasfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.biasfitkeys[keyind])
        # Load dark files if necessary
        if not self.darkloaded or self.getarg('reload'):
            self.loaddark()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.darkfitkeys)):
                if self.darkkeyvalues[keyind] != self.datain.getheadval(self.darkfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.darkfitkeys[keyind])
        # Load flat files if necessary
        if not self.flatloaded or self.getarg('reload'):
            self.loadflat()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.flatfitkeys)):
                if self.flatkeyvalues[keyind] != self.datain.getheadval(self.flatfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.flatfitkeys[keyind])
        # in the config file, set the 'intermediate' variable to either true or false to enable
        # saving of intermediate steps
        saveIntermediateSteps = self.getarg('intermediate')
        self.dataout = DataFits(config=self.datain.config)
        
        image = self.datain.image
        image_exp = self.datain.header['EXPTIME']
        
        # subtract bias from image
        image = self.subtract_bias(image, self.bias)
        if (saveIntermediateSteps == True):
            dataName = "BIAS_SUBTRACT"
            self.dataout.imageset(image, imagename=dataName)
            self.dataout.setheadval('DATATYPE','IMAGE',
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                    dataname=dataName)
        
        # subtract dark from image
        image = self.subtract_dark(image, self.dark, scale=True,
                                   img_exposure=image_exp,
                                   dark_exposure=self.dark_exp_length)
        if (saveIntermediateSteps == True):
            dataName = "DARK_SUBTRACT"
            self.dataout.imageset(image, imagename=dataName)
            self.dataout.setheadval('DATATYPE','IMAGE',
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname,
                                    dataname=dataName)
        
        # apply flat correction to image
        image = self.flat_correct(image, self.flat)

        # if separating bias,dark,flat steps, save the flat-corrected portion
        # into its own hdu
        if (saveIntermediateSteps == True):
            dataName = "FLAT_CORRECT"
            self.dataout.imageset(image, imagename=dataName)
            self.dataout.setheadval('DATATYPE','IMAGE',
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname,
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname,
                                    dataname=dataName)
            self.dataout.setheadval('HISTORY','FLAT: %s' % self.flatname,
                                    dataname=dataName)
        else:
            # copy calibrated image into self.dataout
            self.dataout.image = image
            self.dataout.header = self.datain.header
            ### Finish - cleanup
            # Update DATATYPE
            self.dataout.setheadval('DATATYPE','IMAGE')
            # Add bias, dark files to History
            self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname)
            self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname)
            self.dataout.setheadval('HISTORY','FLAT: %s' % self.flatname)

        self.dataout.filename = self.datain.filename


    # This function is directly lifted from CCDProc
    #     https://github.com/astropy/ccdproc/blob/master/ccdproc/core.py
    # Instead of directly calling CCDProc, we have included the function here
    # to increase educational value and to decrease reliance on external
    # libraries.
    def subtract_bias(self, image, bias):
        """
        Subtract master bias from image.
        Parameters
        ----------
        image : `~numpy.ndarray`
            Image from which bias will be subtracted.
        bias : `~astropy.ndarray`
            Master image to be subtracted from ``image``.
        {log}
        Returns
        -------
        result : `~numpy.ndarray`
            numpy array object with bias subtracted.
        """
        self.log.debug('Subtracting bias...')
        result = image - bias
        self.log.debug('Subtracted bias.')
        return result

    # this code is also lifted from ccdproc
    # https://github.com/astropy/ccdproc/blob/master/ccdproc/core.py
    # some of the code is removed from the original ccdproc because it is not
    # relevant to how SEO currently processes data. If you are looking at this
    # code in the future, there is more code available to draw from
    def subtract_dark(self, image, dark, scale=False, img_exposure=None, dark_exposure=None):
        """
        Subtract dark current from an image.
        Parameters
        ----------
        image : `~numpy.ndarray`
            Image from which dark will be subtracted.
        dark : `~anumpy.ndarray`
            Dark image.
        scale: bool, optional
            If True, scale the dark frame by the exposure times.
            Default is ``False``.
        img_exposure : double, optional
            Exposure time for the base image. Only used if scale = True
            Default is ``None``.
        dark_exposure : double, optional
            Exposure time for the dark image. Only used if scale = True
            Default is ``None``.
        
        {log}
        Returns
        -------
        result : `~numpy.ndarray`
            Dark-subtracted image.
        """
        
        self.log.debug('Subtracting dark...')
        # if dark current is linear, then this first step scales the
        # provided dark to match the exposure time
        if scale:
            dark_scaled = dark.copy()

            dark_scaled = dark_scaled * img_exposure / dark_exposure
            result = image - dark_scaled
        else:
            result = image - dark
    
        self.log.debug('Subtracted dark.')
        return result
    
    # This code is also from ccdproc. A notable removal is the option to 
    # manually choose maximum and minimum flat values.
    def flat_correct(self, image, flat):
        """Correct the image for flat fielding.
        The flat field image is normalized by its mean before flat correcting.
        Parameters
        ----------
        image : `~numpy.ndarray`
            Data to be transformed.
        flat : `~numpy.ndarray`
            Flatfield to apply to the data.
        {log}
        Returns
        -------
        flat_corrected : `~numpy.ndarray`
            Numpy array object with flat corrected.
        """
        self.log.debug('Correcting flat...')
        # Use the min_value to replace any values in the flat
        flat_corrected = image.copy()
        
        flat_mean = flat.mean()

        # Normalize the flat.
        flat_normed = flat / flat_mean

        # divide through the flat
        flat_corrected = image / flat_normed

        self.log.debug('Corrected flat.')
        return flat_corrected
        
    def loadbias(self):
        """ Loads the bias information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        # master bias frame
        # Search for bias and load it into data object
        namelist = self.loadauxname('bias', multi = False)
        self.log.info('File loaded: %s' % namelist)
        if(len(namelist) == 0):
            self.log.error('Bias calibration frame not found.')
            raise RuntimeError('No bias file loaded')
        self.log.debug('Creating master bias frame...')
        # if there is just one, use it as biasfile or else combine
        # all to make a master bias
        biasData = CCDData.read(namelist, unit='adu', relax=True)
        self.bias = biasData.data
        # Finish up
        self.biasloaded = True
        self.biasname = namelist
        self.log.debug('LoadBias: done')

    def loaddark(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            Also bias corrects dark files if not already done.
        """
        #master dark frame
        dark_is_bias_corrected = False
        dark_bias = None
        namelist = self.loadauxname('dark', multi = False)
        if(len(namelist) == 0):
            self.log.error('Dark calibration frame(s) not found.')
            raise RuntimeError('No dark file loaded')
        
        self.log.debug('Creating master dark frame...')
        #if there is just one, use it as darkfile or else combine all to make a master dark
        darkFile = CCDData.read(namelist, unit='adu', relax=True)
        self.dark = darkFile.data
        self.dark_exp_length = darkFile.header['EXPTIME']
        # Finish up 
        self.darkloaded = True 
        self.darkname = namelist
        self.log.debug('LoadDark: done')
            
    def loadflat(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            Also biascorrects and dark corrects flat files if not already done.
        """
        #create master flat frame
        namelist = self.loadauxname('flat', multi = False)
        if(len(namelist) == 0):
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file loaded')
        count = 0
        datalist = []
        flat_corrected = None

        self.log.debug('Creating master flat frame...')
        #if there is just one, use it as flatfile or else combine all to make a master flat
        flatFile = CCDData.read(namelist, unit='adu', relax=True)
        self.flat = flatFile.data

        # Finish up
        self.flatloaded = True  
        self.flatname = namelist 
        self.log.debug('LoadFlat: done')
            
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
2020-07-26 - Saving intermediate BDF steps and removal of CCDProc code - Lorenzo Orders
2018-08-02 - Bias/Dark correction of darks/flats moved to stepmasterbias/dark/flat - Matt Merz
07/28/2017 - Script created by Atreyo Pal
'''