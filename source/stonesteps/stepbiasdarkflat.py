#!/usr/bin/env python
""" 
    Pipestep BiasDarkFlat

    This module defines the pipeline step to calibrate data files w.r.t. bias, 
    dark and flat frames. This pipestep works with pipedata object lists called 
    by StepLoadAux and removes bad pixels from the final RGB image. This is very 
    important to smoothen the background texture of the image and the final image
    looks far more clearer once calibration is done.
    
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
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepParent # pipestep stepparent object
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
        # Load bias files if necessary
        if not self.biasloaded or self.getarg('reload'):
            self.loadbias()
        # Else: check data for correct instrument configuration - currently not in use(need improvement)
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
        self.dataout = DataFits(config=self.datain.config)
        self.dataout.image = image.data
        self.dataout.header = image.header
        self.dataout.filename = self.datain.filename
        ### Finish - cleanup
        # Update DATATYPE
        self.dataout.setheadval('DATATYPE','IMAGE')
        # Add bias, dark files to History
        self.dataout.setheadval('HISTORY','BIAS: %s' % self.biasname)
        self.dataout.setheadval('HISTORY','DARK: %s' % self.darkname)
        self.dataout.setheadval('HISTORY','FLAT: %s' % self.flatname)

    def loadbias(self):
        """ Loads the bias information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            If multiple matching files are found, they are combined into a single 
            master bias frame by ccdproc.
        """
        #master bias frame
        #Search for bias and load it into data object
        namelist = self.loadauxname('bias', multi = False)
        self.log.info('File loaded: %s' % namelist)
        if(len(namelist) == 0):
            self.log.error('Bias calibration frame not found.')
            raise RuntimeError('No bias file loaded')
        self.log.debug('Creating master bias frame...')
        #if there is just one, use it as biasfile or else combine all to make a master bias
        self.bias = ccdproc.CCDData.read(namelist, unit='adu', relax=True)
        # Finish up
        self.biasloaded = True
        self.biasname = namelist
        self.log.debug('LoadBias: done')
        
    def loaddark(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            If multiple matching files are found, they are combined into a single 
            master dark frame by ccdproc.
            Also bias corrects dark files if not already done.
        """
        #master dark frame
        dark_is_bias_corrected = False
        dark_bias = None
        namelist = self.loadauxname('dark', multi = False)
        if(len(namelist) == 0):
            self.log.error('Dark calibration frame(s) not found.')
            raise RuntimeError('No dark file loaded')
        # This has been commented out as it is now in StepMasterDark
        # darks = None
        # for name in namelist:
        #     #is (any) dark file bias corrected?
        #     header = fits.getheader(name)
        #     if(header.get('BIAS') != None):
        #         dark_is_bias_corrected = True
        #         dark_bias = header.get('BIAS')
        #     elif(header.get('BIASCORR') != None):
        #         dark_is_bias_corrected = True
        #         dark_bias = header.get('BIASCORR')
        #     if(darks):
        #         darks += ','+name
        #     else:
        #         darks = name
        self.log.debug('Creating master dark frame...')
        #if there is just one, use it as darkfile or else combine all to make a master dark
        self.dark = ccdproc.CCDData.read(namelist, unit='adu', relax=True)
        #bias correct, if necessary
        # if(not dark_is_bias_corrected):
        #     #Subtracting master bias frame from master dark frame
        #     self.dark = ccdproc.subtract_bias(self.dark, self.bias, add_keyword=False)
        # else:
        #     self.log.debug('Master dark frame is *already* bias corrected (%s).' % dark_bias) 
        # Finish up 
        self.darkloaded = True 
        self.darkname = namelist
        self.log.debug('LoadDark: done')
            
    def loadflat(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
            If multiple matching files are found, they are combined into a single 
            master flat frame by ccdproc.
            Also biascorrects and dark corrects flat files if not already done.
        """
        #create master flat frame
        flat_is_bias_corrected = False
        flat_bias = None
        flat_is_dark_corrected = False
        flat_dark = None
        flat_ave_exptime = 0
        namelist = self.loadauxname('flat', multi = False)
        if(len(namelist) == 0):
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file loaded')
        count = 0
        datalist = []
        flat_corrected = None
        # This has been commented out as it is now in StepMasterFlat
        #check a few things in these flat component frames
        # for name in namelist:
            # header = fits.getheader(name)
            #is this flat bias corrected?
            # if(header.get('BIAS') != None):
            #     flat_is_bias_corrected = True
            #     flat_bias = header.get('BIAS')
            # elif(header.get('BIASCORR') != None):
            #     flat_is_bias_corrected = True
            #     flat_bias = header.get('BIASCORR')
            # #is this flat dark corrected?
            # if(header.get('DARK') != None):
            #     flat_is_dark_corrected = True
            #     flat_dark = header.get('DARK')
            # elif(header.get('DARKCORR') != None):
            #     flat_is_dark_corrected = True
            #     flat_dark = header.get('DARKCORR')
            # flat_corrected = "%s"%(name.rsplit('.',1)[0])+".corrected"	
            # shutil.copy(name, flat_corrected)
            # self.log.debug('Copying %s to %s' % (name, flat_corrected))
            # self.flat = ccdproc.CCDData.read(flat_corrected, unit='adu', relax=True)
            # #bias correct, if necessary
            # if(not flat_is_bias_corrected):
            #     self.log.debug('Subtracting master bias frame from flat frame...')
            #     self.flat = ccdproc.subtract_bias(self.flat, self.bias, add_keyword=False)
            # else:
            #     self.log.debug('Flat frame (%s) is *already* bias corrected.'%flat_bias)
            # #dark correct, if necessary
            # if(not flat_is_dark_corrected):
            #     self.log.debug('Subtracting master dark frame from flat frame...')
            #     self.flat = ccdproc.subtract_dark(self.flat, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second, add_keyword=False)
            # else:
            #     self.log.debug('Flat frame (%s) is *already* dark corrected.'%flat_dark)      
            # #create CCD Data object list with corrected flat files
            # datalist.append(self.flat)
            # #calc average exposure time for potential dark correction
        #     if(header.get('EXPTIME') != None):
        #         try:
        #             exptime = float(header.get('EXPTIME'))
        #             flat_ave_exptime += exptime
        #         except ValueError:
        #             self.log.error('Exposure time (EXPTIME) is not a float (%s).'%(header.get('EXPTIME')))
        #         count += 1
        # #calc average exposure time
        # if(count > 0):
        #     flat_ave_exptime = flat_ave_exptime/count
        #     self.flat.header['EXPTIME'] = flat_ave_exptime
        #     self.log.info("Average exposure time for flats is %f"%flat_ave_exptime)
        self.log.debug('Creating master flat frame...')
        #if there is just one, use it as flatfile or else combine all to make a master flat
        self.flat = ccdproc.CCDData.read(namelist, unit='adu', relax=True)    
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
2018-08-02 - Bias/Dark correction of darks/flats moved to stepmasterbias/dark/flat - Matt Merz
07/28/2017 - Script created by Atreyo Pal
'''
