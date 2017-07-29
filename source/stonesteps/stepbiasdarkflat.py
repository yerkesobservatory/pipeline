#!/usr/bin/env python
""" 
    Pipestep BiasDarkFlat

    This module defines the pipeline step to calibrate data files w.r.t. bias, 
    dark and flat frames. This pipestep works with pipedata object lists called 
    by StepLoadAux and removes bad pixels from the final RGB image.
    
    Author: Atreyo Pal
"""
import os # os library
import sys # sys library
import numpy # numpy library
import logging # logging object library
import shutil # library to provide operations on collections of files
from astropy import units as u
import ccdproc # package for reducing optical CCD telescope data 
from astropy.io import fits #package to recognize FITS files
from drp.pipedata import PipeData # pipeline data object
from drp.stepparent import StepParent # pipestep stepparent object
from drp.steploadaux import StepLoadAux # pipestep steploadaux object

class StepBiasDarkFlat(StepLoadAux, StepParent):
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
        self.bias = [] # CCD data object containing arrays with bias values
        self.biasdata = PipeData() # Pipedata object containing the bias file
        # bias file info and header keywords to fit
        self.biasfile = '' # name of selected bias file
        self.biasfitkeys = [] # FITS keywords that are present in bias      
        self.biaskeyvalues = [] # values of FITS keywords (from data file)  
        # dark values
        self.darkloaded = False # indicates if dark has been loaded
        self.dark = [] # CCD data object containing arrays with dark values
        self.darkdata = PipeData() # Pipedata object containing the dark file
        # dark file info and header keywords to fit
        self.darkfile = '' # name of selected dark file
        self.darkfitkeys = [] # FITS keywords that have to fit for dark     
        self.darkkeyvalues = [] # values of FITS keywords (from data file)  
        # flat values
        self.flatloaded = False # indicates if flat has been loaded
        self.flat = [] # CCD data object containing arrays with flat values
        self.flatdata = PipeData() # Pipedata object containing the flat file
        # flat file info and header keywords to fit
        self.flatfile = '' # name of selected flat file
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
        dark_ave_exptime = self.loadauxname('dark').getheadval('EXPTIME')
        nearexp = {abs(dark_ave_exptime - exp): exp for exp in input_exptime} 
        return nearexp[min(nearexp.keys())]
    '''
   
    def run(self):
        """ Runs the calibrating algorithm. The calibrated data is
            returned in self.dataout
        """
        ### Preparation
        # Load bias files if necessary
        if not self.biasloaded or self.getarg('reload'):
            self.loadbias()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.biasfitkeys)):
                if self.biaskeyvalues[keyind] != self.datain.getheadval(self.biasfitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.biaskeys[keyind])
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
        #convert self.datain to CCD Data object
        image = ccdproc.CCDData(self.datain.image, unit='adu')
        image.header = self.datain.header
        #subtract bias from image    
        image = ccdproc.subtract_bias(image, self.bias, add_keyword=False)
        #subtract dark from image
        image = ccdproc.subtract_dark(image, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second, add_keyword=False)
        #apply flat correction to image
        image = ccdproc.flat_correct(image, self.flat, add_keyword=False)
        # copy calibrated image into self.dataout - make sure self.dataout is a pipedata object
        self.dataout = PipeData(config=self.config)
        self.dataout.image = image.data
        self.dataout.header = image.header
        self.dataout.filename = self.datain.filename
        self.dataout.config = self.datain.config
        ### Finish - cleanup
        # Update DATATYPE
        self.dataout.setheadval('DATATYPE','IMAGE')
        # Add bias, dark files to History
        self.dataout.setheadval('HISTORY','BIAS: %s' % self.dataout.filename)
        self.dataout.setheadval('HISTORY','DARK: %s' % self.dataout.filename)
        self.dataout.setheadval('HISTORY','FLAT: %s' % self.dataout.filename)

    def loadbias(self):
        """ Loads the bias information for the instrument settings
            described in the header of self.datain.

            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        #master bias frame
        #Search for bias and load it into data object
        img = self.loadauxname('bias', multi = True)
        if(len(img) == 0):
            self.log.error('Bias calibration frame not found.')
            raise RuntimeError('No bias file loaded')
        biases = None
        for i in range(0,len(img)):
            if(biases):
                biases += ','+img[i]
            else:
                biases = img[i]
        #if there is just one, make it two of the same for the combine!
        if (len(img) == 1):
            biases += ','+img[0]
        self.log.debug(biases)
        self.log.warn('Creating master bias frame...')
        self.bias = ccdproc.combine(biases, method='median', unit='adu', add_keyword=False)
        # Finish up
        self.biasloaded = True
        self.log.debug('LoadBias: done')
        
    def loaddark(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.

            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        #master dark frame
        dark_is_bias_corrected = False
        dark_bias = None
        img = self.loadauxname('dark', multi = True)
        if(len(img) == 0):
            self.log.error('Dark calibration frame(s) not found.')
            sys.exit(-1)
        darks = None
        for i in range(0,len(img)):
            #is (any) dark file bias corrected?
            header = fits.getheader(img[i])
            if(header.get('BIAS') != None):
                dark_is_bias_corrected = True
                dark_bias = header.get('BIAS')
            elif(header.get('BIASCORR') != None):
                dark_is_bias_corrected = True
                dark_bias = header.get('BIASCORR')
            if(darks):
                darks += ','+img[i]
            else:
                darks = img[i]
        #if there is just one, make it two of the same for the combine!
        if (len(img) == 1):
            darks += ','+img[0]
        self.log.warn('Creating master dark frame...')
        self.dark = ccdproc.combine(darks, method='median', unit='adu', add_keyword=False, **{'verify': 'ignore'})
        #bias correct, if necessary
        if(not dark_is_bias_corrected):
            #Subtracting master bias frame from master dark frame
            self.dark = ccdproc.subtract_bias(self.dark, self.bias, add_keyword=False)
        else:
            self.log.warn('Master dark frame is *already* bias corrected (%s).' % dark_bias) 
        # Finish up 
        self.darkloaded = True 
        self.log.debug('LoadDark: done')
            
    def loadflat(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.

            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        #create master flat frame
        flat_is_bias_corrected = False
        flat_bias = None
        flat_is_dark_corrected = False
        flat_dark = None
        flat_ave_exptime = 0
        img = self.loadauxname('flat', multi = True)
        if(len(img) == 0):
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file loaded')
        flats = None
        count = 0
        flat_corrected = None
        #check a few things in these flat component frames
        for i in range(0,len(img)):
            header = fits.getheader(img[i])
            #is this flat bias corrected?
            if(header.get('BIAS') != None):
                flat_is_bias_corrected = True
                flat_bias = header.get('BIAS')
            elif(header.get('BIASCORR') != None):
                flat_is_bias_corrected = True
                flat_bias = header.get('BIASCORR')
            #is this flat dark corrected?
            if(header.get('DARK') != None):
                flat_is_dark_corrected = True
                flat_dark = header.get('DARK')
            elif(header.get('DARKCORR') != None):
                flat_is_dark_corrected = True
                flat_dark = header.get('DARKCORR')
            flat_corrected = "%s"%(img[i].rsplit('.',1)[0])+".corrected"	
            shutil.copy(img[i], flat_corrected)
            self.flat = ccdproc.CCDData.read(flat_corrected, unit='adu', relax=True)
            #bias correct, if necessary
            if(not flat_is_bias_corrected):
                self.log.warn('Subtracting master bias frame from flat frame...')
                self.flat = ccdproc.subtract_bias(self.flat, self.bias, add_keyword=False)
            else:
                self.log.warn('Flat frame (%s) is *already* bias corrected.'%flat_bias)
            #dark correct, if necessary
            if(not flat_is_dark_corrected):
                #self.log.warn('Subtracting master dark frame from flat frame...')
                self.flat = ccdproc.subtract_dark(self.flat, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second, add_keyword=False)
            else:
                self.log.warn('Flat frame (%s) is *already* dark corrected.'%flat_dark)      
            if(flats):
                flats += ','+flat_corrected
            else:
                flats = flat_corrected
            #calc average exposure time for potential dark correction
            if(header.get('EXPTIME') != None):
                try:
                    exptime = float(header.get('EXPTIME'))
                    flat_ave_exptime += exptime
                except ValueError:
                    self.log.error('Exposure time (EXPTIME) is not a float (%s).'%(header.get('EXPTIME')))
                count += 1
        #calc average exposure time
        if(count > 0):
            flat_ave_exptime = flat_ave_exptime/count
            self.flat.header['EXPTIME'] = flat_ave_exptime
            self.log.warn("Average exposure time for flats is %f"%flat_ave_exptime)
        #if there is just one, make it two of the same!
        if (len(img) == 1):
            flats += ','+flat_corrected       
        self.log.warn('Creating master flat frame...')
        #scale the flat component frames to have the same mean value, 10000.0
        scaling_func = lambda arr: 10000.0/numpy.ma.median(arr)
        #combine them
        self.flat = ccdproc.combine(flats, method='median', scale=scaling_func, unit='adu', add_keyword=False)
        # Finish up
        self.flatloaded = True   
        self.log.debug('LoadFlat: done')
            
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
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
07/28/2017 - Script created by Atreyo Pal
'''
