#!/usr/bin/env python
""" PIPE STEP MASTER DARK - Version 1.0.0

    This module creates master dark frames from inputs from the CMOS camera in HDR mode.
    Should be run with like data inputs (same binning, etc).
    Bias subtracts with matched master bias frames (from LoadAux).
    
        * Notebook for combining RAW.fits or MDARK.fits images into MDARK.fits or MMDARK.fits images, respectively. MDARK images are the mean of a group of RAW images.  MMDARK images are the mean of a group of MDARK images. At present,we assume that the images being combined have similar S/N ratios (that is, the curent algorithm does not weight the images when combining).
        * In addition to the mean dark image (in the primary HDU), each MMDARK.fits file also contain a second HDU with a noise image made by taking the standard deviation along the stack axis of the N input images and a third HDU with a table of image statistics. Notionally, these extra HDUs could be used for quality analysis or to determine weighting factors. Preserving this information may make it possible to eventually discard the RAW data to preserve disk space.

    @author: Carmen Choza, Al Harper, Marc Berthoud
"""

import os # os library
import sys # sys library
import numpy as np # numpy library
import logging # logging object library
#import astropy
from astropy.time import Time
from astropy.stats import mad_std
from astropy.io import fits
from darepype.drp import StepMIParent
from darepype.drp import DataFits

class StepMasterDarkHdr(StepMIParent):
    """ Stone Edge Pipeline Step Master Dark Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepMasterDarkHdr,self).__init__()

        self.log.debug('Init: done')
    
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
        ### Set names
        # Name of the pipeline reduction step
        self.name='masterdarkhdr'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'MDARK'
        
        ### Set parameter list
        # Clear parameter list
        self.paramlist = []
        
        # Append parameters
        self.paramlist.append(['combinemethod','median',
                               'Specifies how the files should be combined - options are median, average, sum'])
        self.paramlist.append(['outputfolder','',
                               'Output directory location - default is the folder of the input files'])
        
        
        # Set logger for this pipe step
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)  
        # Confirm end of setup.
        self.log.debug('Setup: done')

    def run(self):
        """ Runs the combining algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        
        # Create empy list for filenames and headers of loaded frames
        filelist = []
        headlist = []
        
        if len(self.datain) == 0:
            self.log.error('Dark calibration frame not found.')
            raise RuntimeError('No dark file(s) loaded.')
        elif (len(filelist) == 1):
            raise RuntimeError('Only one dark file(s) loaded.')
            
        # Build list of input darks
        for fin in self.datain:
            self.log.debug("Inputting header: %s" % fin.filename)
            filelist.append(fin.filename)
            headlist.append(fin.header)
            
        numfiles = len(filelist)
            
        # Process time and sort lists    
        date_obs = [f.get('date-obs', 'None') for f in headlist]
        t = Time(date_obs, format='isot', scale='utc')
        tsort = np.argsort(t)
        tfiles = [filelist[i] for i in tsort]
        headlist = [headlist[i] for i in tsort]
        utime = [t[i].unix for i in tsort]
        
        self.log.debug("Times logged and file lists sorted")
        
        # Make dataout
        self.dataout = DataFits(config = self.config)
        
        self.log.debug('Creating master dark frame...')
        
        ### Create master frame
        
        darkstack = []        # Build stack of darks
        median = np.zeros((numfiles))          # A group of 1D numpy arrays to hold information about each image in stack
        mean = np.zeros((numfiles))
        std = np.zeros((numfiles))             # These will be used for table construction
        mad = np.zeros((numfiles))

        for i in range(numfiles):
            dark = self.datain[i].image
            mean[i] = np.nanmean(dark)
            median[i] = np.nanmedian(dark)
            mad[i] = mad_std(dark, ignore_nan = True)
            std[i] = np.nanstd(dark)
            darkstack.append(dark)

        self.dark = np.nanmean(darkstack, axis=0)  # Currently takes a simple mean--could become a parameter
        self.noise = np.nanstd(darkstack, axis=0)        # Make a noise image
            
        # Set output header, put image into output
        self.dataout = self.datain[0].copy()
        #self.dataout.header = self.datain[0].header
        self.dataout.imageset(self.dark)
        
        # Put noise image in second HDU
        self.dataout.imageset(self.noise, 'STD')
        
        ### Extract data table information
        
        # Get elapsed time
        utime = np.asarray(utime)
        etime = utime - utime[0]        # Time elapsed from the beginning of the first exposure of the sequence. 
        index = np.arange(numfiles)   # Make a column with the file sequence numbers 

        # From headers:
        ambient = np.zeros((numfiles))
        primary = np.zeros((numfiles))
        secondar = np.zeros((numfiles))
        dewtem1 = np.zeros((numfiles))
        heatsink = np.zeros((numfiles))
        rotation = np.zeros((numfiles))
        setpoint = np.zeros((numfiles))
        drive = np.zeros((numfiles))
        for i in range(numfiles):
            ambient[i] = headlist[i]['ambient']
            primary[i] = headlist[i]['primary']
            secondar[i] = headlist[i]['secondar']
            dewtem1[i] = headlist[i]['dewtem1']
            heatsink[i] = headlist[i]['heatsink']
            rotation[i] = headlist[i]['rotation']
            setpoint[i] = headlist[i]['setpoint']
            drive[i] = headlist[i]['drive']
            
        # Compute difference between images in the stack and the mean or median of the stack
        
        dmedian = np.zeros((numfiles))
        dmean = np.zeros((numfiles))
        dmad = np.zeros((numfiles))
        dstd = np.zeros((numfiles))
        difimage = np.zeros_like(darkstack)
        for i in range(numfiles):
            difimage[i] = darkstack[i] - self.dark
            dmedian[i] = np.nanmedian(difimage[i])
            dmean[i] = np.nanmean(difimage[i])
            dmad[i] = mad_std(difimage[i],ignore_nan=True)
            dstd[i] = np.nanstd(difimage[i])
            
        ### Update output header with stack information
        
        self.dataout.setheadval('notes', 'Master dark image made from a list of dark images.')
        self.dataout.setheadval('notes1', 'Primary HDU is the mean of the input images.')
        self.dataout.setheadval('notes2', 'HDU#2 is a table with statistical and environmental data.')
        self.dataout.setheadval('filelist', 'Input files are those in date+time-range in output filename.')
        self.dataout.setheadval('imagetyp', 'MDARK', 'Image type of the output file')
        self.dataout.setheadval('bzero', 0.0)
        self.dataout.setheadval('oscnmean', 0.0)
        self.dataout.setheadval('reduceby', 'pipeline StepMasterDarkHdr', 'Reduction software')
        self.dataout.setheadval('numfiles', numfiles, 'Number of input files')
        self.dataout.header['ambient'] = np.nanmean(ambient)
        self.dataout.header['primary'] = np.nanmean(primary)
        self.dataout.header['secondar'] = np.nanmean(secondar)
        self.dataout.header['dewtem1'] = np.nanmean(dewtem1)
        self.dataout.header['heatsink'] = np.nanmean(heatsink)
        self.dataout.header['rotation'] = np.nanmean(rotation)
        self.dataout.header['setpoint'] = np.nanmean(setpoint)
        self.dataout.header['drive'] = np.nanmean(drive)

        ### Put derived and header data into a fits table and add it to the output object

        # Make file identifiers
        IDs = []
        for i in range(numfiles):
            fi = filelist[i].split('_')
            IDs.append(fi[4]+'_'+fi[5])
        fileIDs = np.asarray(IDs)

        # Make a list of fits column objects
        cols = []
        cols.append(fits.Column(name='index', format='I', array=index))
        cols.append(fits.Column(name='fileID', format='20A', array=fileIDs))
        cols.append(fits.Column(name='median', format='D', array=median, unit='ADU'))
        cols.append(fits.Column(name='mean', format='D', array=mean, unit='ADU'))
        cols.append(fits.Column(name='std', format='D', array=std, unit='ADU'))
        cols.append(fits.Column(name='mad', format='D', array=mad, unit='ADU'))
        cols.append(fits.Column(name='dmedian', format='D', array=dmedian, unit='ADU'))
        cols.append(fits.Column(name='dmean', format='D', array=dmean, unit='ADU'))
        cols.append(fits.Column(name='dstd', format='D', array=dstd, unit='ADU'))
        cols.append(fits.Column(name='dmad', format='D', array=dmad, unit='ADU'))
        cols.append(fits.Column(name='ambient', format='D', array=ambient, unit='C'))
        cols.append(fits.Column(name='primary', format='D', array=primary, unit='C'))
        cols.append(fits.Column(name='secondar', format='D', array=secondar, unit='C'))
        cols.append(fits.Column(name='dewtem1', format='D', array=dewtem1, unit='C'))
        cols.append(fits.Column(name='elapsed time', format='D', array=etime, unit='seconds'))
        cols.append(fits.Column(name='heatsink', format='D', array=heatsink, unit='C'))
        cols.append(fits.Column(name='rotation', format='D', array=rotation, unit='degrees'))
        cols.append(fits.Column(name='setpoint', format='D', array=setpoint, unit='C'))
        cols.append(fits.Column(name='drive', format='D', array=drive, unit='percent'))

        # Make table
        table_columns = fits.ColDefs(cols)
        table = fits.BinTableHDU.from_columns(table_columns)
        tabhead = table.header

        # Add table to the datafits output object
        self.dataout.tableset(table.data, tablename = 'table', tableheader=tabhead)        
        
        ### Add history
        self.dataout.setheadval('HISTORY','MasterDark: %d files used' % numfiles)
        
        self.log.debug("Filename: %s" % self.dataout.filename)

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
    StepMasterDarkHdr().execute()
        
        
""" === History ===
    2018-07-23 New step created based on StepRGB - Matt Merz
    2018-08-02 Updates to documentation, step functionality - Matt Merz
"""
