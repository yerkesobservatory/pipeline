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
from astropy.nddata import CCDData
# import ccdproc # package for reducing optical CCD telescope data 
from astropy.io import fits #package to recognize FITS files
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepParent # pipestep stepparent object
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object

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
    
    # This function is directly lifted from CCDProc https://github.com/astropy/ccdproc/blob/master/ccdproc/core.py
    # Instead of directly calling CCDProc, we have included the function here to increase educational value
    # and to decrease reliance on external libraries.
    def subtract_bias(self, image, bias):
        """
        Subtract master bias from image.

        Parameters
        ----------
        image : `~astropy.nddata.CCDData`
            Image from which bias will be subtracted.

        bias : `~astropy.nddata.CCDData`
            Master image to be subtracted from ``ccd``.

        {log}

        Returns
        -------
        result : `~astropy.nddata.CCDData`
            CCDData object with bias subtracted.
        """
        self.log.debug('Subtracting bias...')
        result = image.copy()
        try:
            result.data = image.data - bias.data
        # we believe that we should keep this error detection in theory, the bias 
        # and image both come from seo, so their units should be the same
        except ValueError as e:
            if 'operand units' in str(e):
                raise u.UnitsError("Unit '{}' of the uncalibrated image does not "
                                "match unit '{}' of the calibration "
                                "image".format(image.unit, bias.unit))
            else:
                raise e

        self.log.debug('Subtracted bias.')
        return result

    # this code is also lifted from ccdproc https://github.com/astropy/ccdproc/blob/master/ccdproc/core.py
    # some of the code is removed from the original ccdproc because it is not relevant to how SEO currently 
    # processes data. If you are looking at this code in the future, there is more code available to draw from
    def subtract_dark(self, image, dark, scale=False, exposure_time=None, exposure_unit=None):
        """
        Subtract dark current from an image.

        Parameters
        ----------
        image : `~astropy.nddata.CCDData`
            Image from which dark will be subtracted.

        dark : `~astropy.nddata.CCDData`
            Dark image.

        exposure_time : str or `~ccdproc.Keyword` or None, optional
            Name of key in image metadata that contains exposure time.
            Default is ``None``.

        exposure_unit : `~astropy.units.Unit` or None, optional
            Unit of the exposure time if the value in the meta data does not
            include a unit.
            Default is ``None``.

        scale: bool, optional
            If True, scale the dark frame by the exposure times.
            Default is ``False``.

        {log}

        Returns
        -------
        result : `~astropy.nddata.CCDData`
            Dark-subtracted image.
        """
        
        self.log.debug('Subtracting dark...')
        result = image.copy()
        try:
            # if dark current is linear, then this first step scales the provided
            # dark to match the exposure time
            if scale:
                dark_scaled = dark.copy()

                data_exposure = image.header[exposure_time]
                dark_exposure = dark.header[exposure_time]
                # data_exposure and dark_exposure are both quantities,
                # so we can just have subtract do the scaling
                dark_scaled = dark_scaled.multiply(data_exposure /
                                                    dark_exposure)
                result.data = image.data - dark_scaled.data
            else:
                result.data = image.data - dark.data
        except (u.UnitsError, u.UnitConversionError, ValueError) as e:

            # Make the error message a little more explicit than what is returned
            # by default.
            raise u.UnitsError("Unit '{}' of the uncalibrated image does not "
                            "match unit '{}' of the calibration "
                            "image".format(image.unit, dark.unit))
        
        self.log.debug('Subtracted dark.')
        return result
    # This code is also from ccdproc. A notable removal is the option to manually choose
    # maximum and minimum flat values.
    def flat_correct(self, image, flat):
        """Correct the image for flat fielding.

        The flat field image is normalized by its mean or a user-supplied value
        before flat correcting.

        Parameters
        ----------
        ccd : `~astropy.nddata.CCDData`
            Data to be transformed.

        flat : `~astropy.nddata.CCDData`
            Flatfield to apply to the data.

        {log}

        Returns
        -------
        ccd : `~astropy.nddata.CCDData`
            CCDData object with flat corrected.
        """
        self.log.debug('Correcting flat...')
        # Use the min_value to replace any values in the flat
        flat_corrected = image.copy()
        use_flat = flat
        flat_mean_val = use_flat.data.mean()

        # Normalize the flat.
        flat_mean = flat_mean_val * use_flat.unit
        flat_normed = use_flat.divide(flat_mean)

        # divide through the flat
        flat_corrected.data = image.divide(flat_normed)

        self.log.debug('Corrected flat.')
        return flat_corrected

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
        
        image = CCDData(self.datain.image, unit='adu')
        image.header = self.datain.header
        #subtract bias from image
        image = self.subtract_bias(image, self.bias)
        #subtract dark from image
        image = self.subtract_dark(image, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second)
        #apply flat correction to image
        image = self.flat_correct(image, self.flat)
        # copy calibrated image into self.dataout - make sure self.dataout is a pipedata object
        self.dataout = DataFits(config=self.datain.config)
        self.dataout.image = image.data
        self.dataout.header = self.datain.header
        
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
        self.bias = CCDData.read(namelist, unit='adu', relax=True)
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
        self.dark = CCDData.read(namelist, unit='adu', relax=True)
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
        self.flat = CCDData.read(namelist, unit='adu', relax=True)    
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
