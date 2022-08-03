### BUILDING STEPMASTERFLATHDR : COPY INTO .PY FILE AND PUSH TO GITHUB TODAY
from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans # For masking/replacing

import astropy
import ccdproc # package for reducing optical CCD telescope data 
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits #package to recognize FITS files

import scipy.ndimage as nd
import numpy as np
import logging
from skimage.measure import block_reduce

class StepMasterFlatHdr(StepLoadAux, StepMIParent):
    '''Pipeline step object to produce and subtract master flats, pulls in pfits (3d image with one layer as the dark
    and one as the bias) and raw flats, requires valid config file'''
    
    def __init__(self):
            """ Constructor: Initialize data objects and variables
            """
            #Finian: taken from stepHDR ; "PFIT images, which comprise a 3D image with two planes. In the first, 
            #the pixel values are the dark current and in the second they are the bias levels. 
            #The PFIT images enable interpolation to the exact exposure time of the individual RAW flat images." 
            # - make_flat_HDR_auto , Finian: wondering why pfit + darks are loaded if the pfit contains the dark+bias?
            # in meeting Al said pfit is int.+slope of interpolated dark current vs. exptime? confused on this
            # hand-drawn diagram insinuates I only need the pfits + flats ; remove darks from code?            
            # Call superclass constructor (calls setup)
            super(StepMasterFlatHdr,self).__init__()

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
        # FINIAN: taken from stepmasterflat.py , for loadaux lines at bottom should I replace with one line
        # loading pfit?
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
        ### Set Names
        # Name of the pipeline reduction step
        self.name='masterflathdr'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'MFLATHDR'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters !!!! WHAT PARAMETERS ARE NEEDED ????? !!!!!
        self.paramlist.append(['combinemethod','median',
                               'Specifies how the files should be combined - options are median, average, sum'])
        self.paramlist.append(['outputfolder','',
                               'Output directory location - default is the folder of the input files'])
        # Get parameters for StepLoadAux, replace auxfile with biasfile
        self.loadauxsetup('bias')
        # Get parameters for StepLoadAux, replace auxfile with darkfile
        self.loadauxsetup('dark')

    def run(self):
        # FINIAN: this code from stepmasterflat.py was more legible to me than the code in stephdr.py
        # can I replace the separate dark and bias steps with a single pfit step?
        """ Runs the combining algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        # Find master dark to subtract from master dark
        biaslist = self.loadauxname('bias', multi = False)
        darklist = self.loadauxname('dark', multi = False)
        if(len(biaslist) == 0):
            self.log.error('No bias calibration frames found.')
        if(len(darklist) == 0):
            self.log.error('No bias calibration frames found.')
        self.bias = ccdproc.CCDData.read(biaslist, unit='adu', relax=True)
        self.dark = ccdproc.CCDData.read(darklist, unit='adu', relax=True)
        # Create empy list for filenames of loaded frames
        filelist=[]
        for fin in self.datain:
            self.log.debug("Input filename = %s" % fin.filename)
            filelist.append(fin.filename)
        # Make a dummy dataout
        self.dataout = DataFits(config = self.config)
        if len(self.datain) == 0:
            self.log.error('Flat calibration frame not found.')
            raise RuntimeError('No flat file(s) loaded')
        self.log.debug('Creating master flat frame...')
        # Create master frame: if there is just one file, turn it into master bias or else combine all to make master bias
        if (len(filelist) == 1):
            self.flat = ccdproc.CCDData.read(filelist[0], unit='adu', relax=True)
            self.flat = ccdproc.subtract_bias(self.flat, self.bias, add_keyword=False)
            self.flat = ccdproc.subtract_dark(self.flat, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second, add_keyword=False)
        else:
            #bias and dark correct frames
            flatlist=[]
            for i in filelist:
                flat =ccdproc.CCDData.read(i, unit='adu', relax=True)
                flatsubbias = ccdproc.subtract_bias(flat, self.bias, add_keyword=False)
                flatsubbiasdark = ccdproc.subtract_dark(flatsubbias, self.dark, scale=True, exposure_time='EXPTIME', exposure_unit=u.second, add_keyword=False)
                flatlist.append(flatsubbiasdark)
            #scale the flat component frames to have the same mean value, 10000.0
            scaling_func = lambda arr: 10000.0/numpy.ma.median(arr)
            #combine them
            self.flat = ccdproc.combine(flatlist, method=self.getarg('combinemethod'), scale=scaling_func, unit='adu', add_keyword=False)
        # set output header, put image into output
        self.dataout.header=self.datain[0].header
        self.dataout.imageset(self.flat)
        # rename output filename
        outputfolder = self.getarg('outputfolder')
        if outputfolder != '':
            outputfolder = os.path.expandvars(outputfolder)
            self.dataout.filename = os.path.join(outputfolder, os.path.split(filelist[0])[1])
        else:
            self.dataout.filename = filelist[0]
        # Add history
        self.dataout.setheadval('HISTORY','MasterFlat: %d files used' % len(filelist))

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
    StepMasterFlatHdr().execute()
