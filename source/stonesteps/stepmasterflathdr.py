'''
STEP DESCRIPTION:

Makes flat-field calibration images. Specialized for making HDR images for the 
SBIG Aluma AC4040 CMOS camera. Also makes gain ratio images by comparing signals in
the high and low gain channels when the camera is operated in 
HDR mode. In HDR mode, pixels are simultaneously read through two amplifiers with
different gains. After bias-subtraction, the total signal (sky signal plus dark signal)
of the high gain channel is divided by the total signal of the low-gain channel.

Input files are:
-   RAW flat images, taken at dusk with a script that acquires flat images of
    approximately equal exposure level. The script attempts to make the median value
    of the high-gain image approximately 3000 ADU (the high-gain channel saturates at
    4095 ADU).
-   PFIT images, which comprise a 3D image with two planes. In the first, the pixel
    values are the dark current (the slope of a linear fit to a dataset with exposure
    times between 1 and 64 seconds). The second plane contains the bias levels (the
    intercepts of the linear fits). The PFIT images enable interpolation to the exact
    exposure time of the individual RAW flat images. Separate PFIT files are
    required for high- and low-gain flat images.
    
Output is a single file with the suffix HFLAT. It has three HDUs:
-   The primary HDU is the flat-field image. It is a 3D image with two planes.
    The first plane (index = 0) is the high-gain flat and the second (index = 1)
    is the low-gain flat.
-   The second HDU is the gain ratio image. The image name is 'GAIN RATIO'.
-   The third HDU is a table with statistical information about the individual RAW
    images used to make the MFLAT image and instrumental and environmental parameters
    taken from the RAW image headers. The table name is 'TABLE'.

The actual input data to StepMasterFlatHdr are DataFits objects created by StepDataGroup.
StepDatagroup selects all RAW flat files within a directory having the same GAIN and
FILTER keywords.

The output flat image and gain image both contain NANS that identify pixels with
high dark current or anomalous gains. The threshold parameters
can be tuned by changing parameters in the config file.
'''

'''
FINIAN'S NOTES
Pipestep Master Flat HDR
Produces Master Flat images using pfit images and raw flats as inputs. 
Uses StepLoadAux to load in pfit files. 
FINIAN'S NOTES:
- Some imports are currently redundant due to being drawn from multiple sources, this can be cleaned up later.
- Potentially irrelevant code for loading in flats is in the __init__, there is code for loading them in the Run 
  funtion from Al's make_flat_HDR colab notebook. I believe this should be taken care of by the load input part
  of the config file, but am not sure exactly what to put in this notebook to load in the flats.
  
- Some functions from Al's make_flat notebook may be irrelevant and able to be removed.
- Unsure if the code to build the hotpix mask should be in here, and if so, where. 
  Al suggested moving the code that generates this mask to the code that creates the pfit, and loading it in
  using the loadaux functionality.
  
#Finian: taken from stepHDR ; "PFIT images, which comprise a 3D image with two planes. In the first, 
#the pixel values are the dark current and in the second they are the bias levels. 
#The PFIT images enable interpolation to the exact exposure time of the individual RAW flat images." 
# WILL NEED HI GAIN AND LO GAIN FLAT, LIKEWISE FOR PFIT
'''

from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class
from astropy.io import fits #package to recognize FITS files
import numpy as np
import logging
from astropy.stats import mad_std    # The median absolute deviation, a more robust estimator than std.
from astropy.time import Time
#?import sys
import os
#?import astropy
#? from astropy import units as u

class StepMasterFlatHdr(StepLoadAux, StepMIParent):
    '''
    Pipe step object to create master flats for HDR images taken with SBIG CMOS camera.
    
    Inputs are sets of high-gain (bin1H), and low-gain (bin1L) RAW flat images, all taken
    through the same spectral filter and electronic gain.
    
    Required auxiliary files include two PFIT files, one each for high and low gain data.
    PFIT files contain data required to synthesize dark images exactly matching the
    exposure times of the input flat images. PFIT files contain a 3D image with a dark
    current image in the first plane and a bias image in the second plane.
    
    Output is a single file that contains a two-plane 3D imgage in the primary HDU.
    The first plane is a high-gain flat and the second plane is a low-gain flat.
    A second HDU contains a relative gain image (the ratio of the electronic gains of
    the pixels in the high and low-gain images. A third HDU contains a table comprising
    statistical information about the input images and instrumental and environmental
    conditions during the exposures.
    '''
    
    def __init__(self):
            """ Constructor: Initialize data objects and variables
            """

            
            # Call superclass constructor (calls setup)
            super(StepMasterFlatHdr,self).__init__()

            self.hpfitloaded = False   # indicates if high-gain PFIT file has been loaded
            self.hpfit = None          # Numpy array object containing PFIT data
            self.hpfitname = ''        # name of selected PFIT file
            
            self.lpfitloaded = False   # indicates if low-gain PFIT file has been loaded
            self.lpfit = None          # Numpy array object containing PFIT data
            self.lpfitname = ''        # name of selected PFIT file

            self.flatloaded = False    # Indicates if flat file has been loaded
            self.flat = None           # Numpy array object containing flat values
            self.flatname = ''         # name of selected flat file

            # Finish up.
            self.log.debug('Init: done')

    def setup(self):
        '''
        Names and Parameters need to be set here
        Sets the internal names for the function and for saved files.
        Defines the input parameters for the current pipe step.
        Setup() is called at the end of __init__
        The parameters are stored in a list of entries containing the
        following information:
            - name:    The name for the parameter. This name is used when
                       calling the pipe step from command line or python shell.
                       It is also used to identify the parameter in the pipeline
                       configuration file.
            - default: A default value for the parameter. If nothing, set ''
                       for strings, 0 for integers and 0.0 for floats
            - help:    A short description of the parameter.
        '''
        ### Set Names
        # Name of the pipeline reduction step
        self.name='masterflathdr'
        # Shortcut for pipeline reduction step and identifier for saved file names.
        self.procname = 'HDRFLAT'
        
        ### Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['combinemethod','median',
                               'Specifies how the files should be combined - options are median, average, sum'])
        self.paramlist.append(['outputfolder','',
                               'Output directory location - default is folder with input files'])
        self.paramlist.append(['logainlim', 1.0, 'What to subtract from median gain to define low mask'])
        self.paramlist.append(['higainlim', 1.0, 'What to add to median gain to define high mask'])
        self.paramlist.append(['hotpxlim', 99.5, 'Hot pixel limit percentile'])
        self.paramlist.append(['reload', False,'Set to True to look for new pfit files for every input'])
        self.paramlist.append(['outputfolder', '','Output directory path'])
        ### Set parameters for StepLoadAux
        self.loadauxsetup('lpfit')
        self.loadauxsetup('hpfit')

        
    def timesortHDR(self, datalist, date_key = 'date-obs'):
        '''
          Sorts a list of fits files by a header keyword with a date/time value. In the case of an
          SBIG CMOS camera RAW file, the date/time is read from the second HDU.
          Arguments:
              datalist   = a list of DataFits objects
              date_key   = the header keyword containing the time/date data
          Returns:
              tfiles     = the sorted file list
              utime      = a list of the unix times of the observations
          Author(s): Al Harper, Finian Ashmead
          Modified: 210807, 210815, 220802 (modified for use in pipe step)
          Version: 1.1 (pipe)
        '''

        date_obs = []                                   # Make a list to hold the date-obs keyword strings.
        for d in datalist:
            if '_bin1L' in d.filename and '_RAW.' in d.filename:
                head = d.getheader(d.imgnames[1])       # Get the header of the second HDU (index = [1]).
                date_obs.append(head[date_key])         # Add date information to list. of string objects.
            else:
                head = d.getheader()                    # Get the header of the primary HDU (index = [0]).
                date_obs.append(head[date_key])         # Add date information to list. of string objects.
        t = Time(date_obs, format='isot', scale='utc')  # Make an astropy time object in 'isot' format.  
        tsort = np.argsort(t)                           # Make a list of indices that will sort by date_obs.
        tfiles = []
        utime = []
        for i in tsort:
            tfiles.append(datalist[i])
            utime.append(t[i].unix)
        return tfiles, utime
        
    def run(self):
        '''
        Process sets of high-gain and low-gain RAW flatfiles (datain).
        This code is based on Al's Jupyter notebook make_flat_HDR_auto_54.
        '''
        '''
        Define boolean to enable print statements for debugging. Set = True to print.
        '''
        pt = False
        
        '''
        Load PFIT files (3D images containing slopes and intercepts of dark current.
        First image plane contains pixel dark currents. Second contains intercepts of
        linear fits to dark current (a bias image).
        '''
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

        '''
        Get image data from low- and high-gain PFIT files and load it into
        a 4D image.
        '''
        hpfitimg = self.hpfit.image   # high-gain PFIT image
        lpfitimg = self.lpfit.image   # low-gain PFIT image
        if pt: print('shapes of hpfitimg and lpfitimg images =', hpfitimg.shape, lpfitimg.shape)
        polydarkimg = np.asarray([hpfitimg, lpfitimg])
        rows, cols = polydarkimg.shape[2], polydarkimg.shape[3]
        '''
        Prepare a hot pixel mask
        '''
        hotpxlim = self.getarg('hotpxlim')
        if pt: print('hotpxlim =', hotpxlim)
        img = polydarkimg[0,0]
        uppercut = np.percentile(img, hotpxlim)  # Value of img at hotpxlim percentile
        if pt: print('hot pixel cutoff =', uppercut)
        hotpix = np.where(img > uppercut)
        if pt: print('number of hot pixels =', len(hotpix[0]))
        if pt: print('')

        '''
        Make two separate lists of input DataFits objects, one for high gain
        images and one for low gain images. Then sort files by exposure time in
        order to get matched lists of high and low gain images associated with
        the same identical exposure.
        '''
        if pt: print('length of datain =', len(self.datain))
        lowgainlist = []
        highgainlist = []
        for f in self.datain:
            if '.fit' in f.filename and '_bin1L' in f.filename:
                lowgainlist.append(f)
            elif '.fit' in f.filename and '_bin1H' in f.filename:
                highgainlist.append(f)
        highgainlist, utimeH = self.timesortHDR(highgainlist, date_key = 'date-obs')
        lowgainlist, utimeL = self.timesortHDR(lowgainlist, date_key = 'date-obs')
        if pt: print('highgainlength, lowgainlength =', len(highgainlist), len(lowgainlist))

        '''
        Make a 4D image of high- and low-gain flat files.
        Dimensions are [gain, stack, row, col].
        '''
        numflats = len(highgainlist)
        flatheadlist, flatstats = [[], []], [[], []]
        flatimage = np.zeros((2, numflats, rows, cols))
        
        '''
        Fill arrays with images and statistical information from
        high-gain RAW flat images.
        '''
        median = np.zeros((numflats))     # 1D numpy array to hold array medians.
        mean = np.zeros((numflats))       # 1D numpy array to hold array medians.
        std = np.zeros((numflats))        # 1D numpy array to hold array stds.
        mad = np.zeros((numflats))        # 1D numpy array to hold array mad_stds.        
        for j in range(numflats):
            flatimage[0, j] = highgainlist[j].image[:,:4096]
            mad[j] = mad_std(flatimage[0, j],ignore_nan=True)
            median[j] = np.nanmedian(flatimage[0, j])
            mean[j] = np.nanmean(flatimage[0, j])
            std[j] = np.nanstd(flatimage[0, j])
        flatstats[0] = {'median':median, 'mean':mean, 'std':std, 'mad':mad}
        
        '''
        Fill arrays with images and statistical information from
        low-gain RAW flat images.
        '''         
        median = np.zeros((numflats))     # 1D numpy array to hold array medians.
        mean = np.zeros((numflats))       # 1D numpy array to hold array medians.
        std = np.zeros((numflats))        # 1D numpy array to hold array stds.
        mad = np.zeros((numflats))        # 1D numpy array to hold array mad_stds.
        for j in range(len(lowgainlist)):
            flatimage[1, j] = lowgainlist[j].imgdata[1][:,:4096]
            mad[j] = mad_std(flatimage[1, j],ignore_nan=True)
            median[j] = np.nanmedian(flatimage[1, j])
            mean[j] = np.nanmean(flatimage[1, j])
            std[j] = np.nanstd(flatimage[1, j])
        flatstats[1] = {'median':median, 'mean':mean, 'std':std, 'mad':mad}

        '''
        Extract header from first high-gain DataFits input object to use
        as the basis for the output object header.
        '''
        flatbaseheader = highgainlist[0].header.copy()

        '''
        Compile a list of exposure times. Set keyword variables for max and min exposure
        times for inclusion in header.
        '''
        flatexptimes = []
        imedian = flatstats[0]['median']
        for j in range(len(highgainlist)):
            exptime = highgainlist[j].header['exptime']
            flatexptimes.append(exptime)
        xtimemin = np.min(flatexptimes)
        xtimemax = np.max(flatexptimes)
        if pt: print('xtimemin =', xtimemin, '  xtimemax =', xtimemax)

        '''
        Construct a stack of dark images to match the flat image exposure times.
        '''
        darkimage = np.zeros_like(flatimage)
        for i in range(2):
            for j in range(flatimage.shape[1]):
                darkimage[i,j] = polydarkimg[i,1] + polydarkimg[i,0] * flatexptimes[j]

        '''
        Derive relative gain images for each HDR flat image in the stack. Note that
        we are measuring the total signal above bias resulting from both photons and
        dark current.
        '''
        biasimage = np.zeros((2, rows, cols))
        biasimage[0] = polydarkimg[0, 1]
        biasimage[1] = polydarkimg[1, 1]
        gainratiostack = np.zeros_like(flatimage[0])
        gainratiomedians = np.zeros((numflats))
        if pt: print('Medians of gain ratio images.')
        for i in range(len(highgainlist)):
            gainratiostack[i] = (flatimage[0][i] - biasimage[0]) / (flatimage[1][i] - biasimage[1])
            gainratiomedians[i] = np.nanmedian(gainratiostack[i])
            if pt: print(gainratiomedians[i])
        gainmean, gainmax, gainmin = np.nanmean(gainratiomedians), np.nanmax(gainratiomedians), np.nanmin(gainratiomedians)
        gainpcnt = (gainmax - gainmin)*100/gainmean
        if pt: print('Shape of gain ratio images:', gainratiostack[0].shape)
        if pt: print('')

        '''
        Create median gain ratio image and std gain ratio image from the stack of gain 
        ratio images. Calculate some statistics of the two images, ignoring 
        any nans.
        '''
        # Take the median and std along the stack axis of gainratiostack.
        mediangainratioimg = np.nanmedian(gainratiostack, axis = 0)
        stdgainratioimg = np.nanstd(gainratiostack, axis = 0)
        # Take image medians and means of the two 2D images derived above.
        grat_median = np.nanmedian(mediangainratioimg)
        grat_mean = np.nanmean(mediangainratioimg)
        stdgrat_median = np.nanmedian(stdgainratioimg)
        stdgrat_mean = np.nanmean(stdgainratioimg)
        if pt: print('median, mean of median gain ratio image =', grat_median, grat_mean)
        if pt: print('max, min =', np.nanmax(mediangainratioimg), np.nanmin(mediangainratioimg))
        if pt: print('median, mean of std gain ratio image =', stdgrat_median, stdgrat_mean)
        if pt: print('max, min =', np.nanmax(stdgainratioimg), np.nanmin(stdgainratioimg))
        if pt: print('')

        '''
        Create masked median gain ratio image and masked std image with nans replacing 
        masked pixels. Mask used here is hotpix (a hot pixel mask made from a
        dark-current image).
        '''
        maskedgainratioimg = mediangainratioimg.copy()
        maskedgainratioimg[hotpix] = np.nan
        maskedstdgainratioimg = stdgainratioimg.copy()
        maskedstdgainratioimg[hotpix] = np.nan

        '''
        Make two masks from the masked gain ratio image in order to eliminate outliers
        with anomalously high or low gain. The default low and high limits have 
        been chosen empirically to eliminate a small fraction of the total number of
        pixels while preserving the smooth (and presumably "physical") central portion of
        the histogram.
        '''
        img = maskedgainratioimg
        logainlim, higainlim = self.getarg('logainlim'), self.getarg('higainlim')
        masklow, maskhigh = grat_median - logainlim, grat_median + higainlim
        gainmasklow, gainmaskhigh = np.where(img > maskhigh), np.where(img < masklow)
        if pt: print('masklow', masklow)
        if pt: print('masklow shape', gainmasklow[0].shape)
        if pt: print('maskhigh', maskhigh)
        if pt: print('maskhigh shape', gainmaskhigh[1].shape)
        if pt: print('')

        '''
        Subtract interpolated darks from the flat images.
        '''
        if '_RAW.fit' in highgainlist[0].filename:
            flatimageDS = flatimage - darkimage

        '''
        Normalize each of the images in the high and low-gain flat image stack to its
        own median (after applying hotpix and gain masks). From this point on,
        there will be nans for every pixel in the hotpix and gain masks.
        '''
        flatimageDSN = np.zeros_like(flatimageDS)
        flatmediansDS = np.zeros((2, numflats))
        for i in range(2):
            for j in range(numflats):
                flatimageDS[i,j][hotpix] = np.nan
                flatimageDS[i,j][gainmasklow] = np.nan
                flatimageDS[i,j][gainmaskhigh] = np.nan
                flatmediansDS[i,j] = np.nanmedian(flatimageDS[i,j])
                flatimageDSN[i,j] = flatimageDS[i,j] / flatmediansDS[i,j]

        '''
        Make median high and low gain flats from the stacks and compute their
        respective median, mean, std, and mad. The result is a 3D image with
        high- and low-gain planes.
        '''
        flat = np.nanmedian(flatimageDSN, axis=1)
        flatmedian, flatmean, flatstd, flatmadstd = [0,0], [0,0], [0,0], [0,0]
        for i in range(2):  
            flatmedian[i] = np.nanmedian(flat[i])            # Median of the median flat.
            flatmean[i] = np.nanmean(flat[i])                # Mean of the median flat.
            flatstd[i] = np.nanstd(flat[i])                  # mad_std of the median flat.
            flatmadstd[i] = mad_std(flat[i],ignore_nan=True) # mad_std of the median flat.
        if pt: print('Median flat median =', flatmedian )
        if pt: print('Mean flat median =', flatmean )
        if pt: print('Median flat std =', flatstd)
        if pt: print('Median flat mad_std =', flatmadstd)
        if pt: print('Shape of flatimgage:', flatimage.shape)
        if pt: print('Shape of flat:', flat.shape)
        if pt: print('')

        '''
        if pt: print a list of the medians of the dark-subtracted high and low gain flat images.
        '''
        if pt: print('List of medians of dark-subtracted high and low gain flat images')
        for i in range(flatimage.shape[1]):
            if pt: print('{:<3}{:<10.2f}{:<10.2f}'.format(i, flatmediansDS[0,i], flatmediansDS[1,i]))
        if pt: print('')

        '''
        Find minimum and maximum values of the image medians for inclusion in output header.
        '''
        medmin, medmax = [0,0], [0,0]
        medmin[0] = np.nanmin(flatmediansDS[0])
        medmax[0] = np.nanmax(flatmediansDS[0])
        medmin[1] = np.nanmin(flatmediansDS[1])
        medmax[1] = np.nanmax(flatmediansDS[1])
        if pt: print('minimum medians =', medmin)
        if pt: print('maximum medians =', medmax)
        if pt: print('')

        '''
        Re-normalize the median flats (H and L) to their medians and mask out pixels with 
        gains greater than or less than 1.0 from the median gain. Replace the masked pixels 
        with np.nan. Hence, when a sky image is divided by the flat, those pixels will also 
        be masked (will be np.nan) in the flat-fielded sky image (in addition to any pixela
        aleady replaced with nans using the hotpix mask).
        '''
        mflat = np.zeros_like(flat)
        mflatmedian = np.zeros((2))
        mflatmadstd = np.zeros((2))
        for i in range(2):
            mflat[i] = flat[i] / flatmedian[i]
            mflatmedian[i] = np.nanmedian(mflat[i])
            mflat[i][gainmasklow] = np.nan
            mflat[i][gainmaskhigh] = np.nan
            mflatmadstd[i] = mad_std(mflat[i], ignore_nan=True)
        if pt: print('mflatmedian, mflatmadstd =', mflatmedian, mflatmadstd)
        if pt: print('')

        '''
        Compute the differences of each of the individual normalized flat images
        from the median flat of the stack. Compute the median, std, and mad_std of
        each of the difference images. Define some statistical measures of the
        uniformity of the difference images that can be reported in header keywords.
        '''
        g = 0           # Gain = 0 for high-gain, 1 for low-gain.
        vmx = mflatmadstd[g] * 2.0
        vmn = - mflatmadstd[g] * 2.0
        difimage = np.zeros_like(flatimageDS)  # Ratio of inidvidual flat images to the median flat image.
        difmadstd = np.zeros((numflats))
        difstd = np.zeros((numflats))
        difmean = np.zeros((numflats))
        difmedian = np.zeros((numflats))
        for i in range(flatimageDS.shape[1]):
            difimage[g,i] = flatimageDSN[g,i]-mflat[g]
            difmadstd[i] = mad_std(difimage[g,i],ignore_nan=True)
            difstd[i] = np.nanstd(difimage[g,i])
            difmean[i] = np.nanmean(difimage[g,i])
            difmedian[i] = np.nanmedian(difimage[g,i])            
        dstdmean, dstdmax, dstdmin = np.nanmean(difstd), np.nanmax(difstd), np.nanmin(difstd)
        dstdpcnt = (dstdmax - dstdmin)*100/dstdmean
        
        '''
        Create output DataFits object and fill with image data.
        '''
        self.dataout = DataFits(config=self.config)
        self.dataout.header = flatbaseheader.copy()
        self.dataout.image = mflat
        self.dataout.imageset(maskedgainratioimg, 'GAIN RATIO')
        
        '''
        Create output file name.
        '''
        if pt: print('length of highgainlist =', len(highgainlist))
        infolder = os.path.split(highgainlist[0].filename)[0]
        lfname = os.path.split(highgainlist[-1].filename)[1]
        ffname = os.path.split(highgainlist[0].filename)[1]
        if pt: print('lfname =', lfname)
        if pt: print('ffname =', ffname)
        lf = lfname.split('_')  # Last filename of time-sorted list.
        ff = ffname.split('_')   # First filename of time-sorted list.
        flatname = 'mflat_'+ff[1]+'_'+ff[3]+'DR'+'_'+ff[4]+'_'+ff[5]+'-'+lf[5]+'_'+ff[6]+'_'+ff[7]+'_'+'XXX.fits'
        if pt: print('flatname =',flatname)
        if pt: print('')
        
        '''
        Rename output filename
        '''
        outputfolder = self.getarg('outputfolder')
        if outputfolder != '':
            outputfolder = os.path.expandvars(outputfolder)
            self.dataout.filename = os.path.join(outputfolder, flatname)
        else:
            self.dataout.filename = os.path.join(infolder, flatname)
        
        
        '''
        Create numpy arrays with information about input data from headers
        and timesortDF functions.
        '''
        # From statistical calculations on input data:
        imedianH, imadH, imeanH, istdH = flatstats[0]['median'], flatstats[0]['mad']\
                                        , flatstats[0]['mean'], flatstats[0]['std']
        imedianL, imadL, imeanL, istdL = flatstats[1]['median'], flatstats[1]['mad']\
                                        , flatstats[1]['mean'], flatstats[1]['std']                          
        # From time_sortDF:
        utime = np.asarray(utimeH)
        etime = utime - utime[0]   # Time elepased from beginning of first exposure of sequence.
        # From header:
        ambient = np.zeros((numflats))
        primary = np.zeros((numflats))
        secondar = np.zeros((numflats))
        dewtem1 = np.zeros((numflats))
        for i in range(numflats):
            ambient[i] = highgainlist[i].header['ambient']
            primary[i] = highgainlist[i].header['primary']
            secondar[i] = highgainlist[i].header['secondar']
            dewtem1[i] = highgainlist[i].header['dewtem1']
        # Make a column with the file sequence numbers.
        index = np.arange(numflats)
        
        '''
        Now put derived and header data into a fits table and add it to the output object.
        '''
        
        '''
        Make file identifiers. The file identifiers are the date and time fields
        of the input filenames.
        '''
        IDs = []
        for i in range(numflats):
            fi = highgainlist[i].filename.split('_')
            IDs.append(fi[4]+'_'+fi[5])
        fileIDs = np.asarray(IDs)

        '''
        Make a list of fits column objects.
        '''
        tcols = []
        tcols.append(fits.Column(name='index', format='I', array=index))
        tcols.append(fits.Column(name='fileID', format='20A', array=fileIDs))
        tcols.append(fits.Column(name='median', format='D', array=imedianH, unit='ADU'))
        tcols.append(fits.Column(name='mean', format='D', array=imeanH, unit='ADU'))
        tcols.append(fits.Column(name='std', format='D', array=istdH, unit='ADU'))
        tcols.append(fits.Column(name='mad', format='D', array=imadH, unit='ADU'))
        tcols.append(fits.Column(name='dmedian', format='D', array=difmedian, unit='ADU'))
        tcols.append(fits.Column(name='dmean', format='D', array=difmean, unit='ADU'))
        tcols.append(fits.Column(name='dstd', format='D', array=difstd, unit='ADU'))
        tcols.append(fits.Column(name='dmad', format='D', array=difmadstd, unit='ADU'))
        tcols.append(fits.Column(name='ambient', format='D', array=ambient, unit='C'))
        tcols.append(fits.Column(name='primary', format='D', array=primary, unit='C'))
        tcols.append(fits.Column(name='secondar', format='D', array=secondar, unit='C'))
        tcols.append(fits.Column(name='dewtem1', format='D', array=dewtem1, unit='C'))
        tcols.append(fits.Column(name='elapsed time', format='D', array=etime, unit='seconds'))
        tcols.append(fits.Column(name='gain ratio', format='D', array=gainratiomedians, unit='ADU'))

        '''
        Make table and add it to the output object.
        '''
        c = fits.ColDefs(tcols)
        table = fits.BinTableHDU.from_columns(c)
        tabhead = table.header
        self.dataout.tableset(table.data, tablename = 'table', tableheader=tabhead)

        '''
        Populate output header with new keyword data.
        '''
        self.dataout.header['notes'] = '1st HDU: 3D mflat image (img[0]=higain, img[1]=logain)'
        self.dataout.header['notes2'] = '2nd HDU: gain ratio image (higain/logain)'
        self.dataout.header['notes3'] = 'Table HDU: statistical and environmental data'
        self.dataout.header['imagetyp'] = 'HDR MFLAT'
        self.dataout.header['bzero'] = 0.0
        self.dataout.header['ambient'] = np.nanmean(ambient)
        self.dataout.header['primary'] = np.nanmean(primary)
        self.dataout.header['secondar'] = np.nanmean(secondar)
        self.dataout.header['dewtem1'] = np.nanmean(dewtem1)        
        
        self.dataout.setheadval('xtimemin', xtimemin, 'Minimum exposure in the set of flats')
        self.dataout.setheadval('xtimemax', xtimemax, 'Maximum exposure in the set of flats')
        self.dataout.setheadval('medminH', medmin[0], 'Minimum high gain median in the set of flats')
        self.dataout.setheadval('medmaxH', medmax[0], 'Maximum high gain median in the set of flats')
        self.dataout.setheadval('medminL', medmin[1], 'Minimum low gain median in the set of flats')
        self.dataout.setheadval('medmaxL', medmax[1], 'Maximum low gain median in the set of flats')
        self.dataout.setheadval('dstdmean', dstdmean, 'Mean std in the set of difference images')
        self.dataout.setheadval('dstdpcnt', dstdpcnt, '(dstdmax-dstdmin)*100.0/dstdmean')
        self.dataout.setheadval('gainmean', grat_mean, 'Mean gain ratio')
        self.dataout.setheadval('gainpcnt', gainpcnt, '(gainmax-gainmin)*100.0/gainmean')
        self.dataout.setheadval('numfiles', numflats, 'Number of RAW exposures in input datasets')
        self.dataout.setheadval('lgainlim', logainlim, 'Lower limit for unmasked pixel gains')
        self.dataout.setheadval('hgainlim', higainlim, 'Upper limit for unmasked pixel gains')
        self.dataout.setheadval('hotpxlim', hotpxlim, 'Upper limit percentile for unmasked dark current')
        
        '''
        Add a history keyword
        '''
        self.dataout.setheadval('HISTORY','HDR Master Flat: made from %d x 2 files' % numflats)
        if pt: print(self.dataout.header)
        
        
        
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
