### BUILDING STEPMASTERFLATHDR
'''
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
STEP DESCRIPTION:
Makes flat-field calibration images. Specialized for making HDR images for the SBIG CMOS camera.
Also makes gain ratio images by comparing signals in the high and low gain channels when the camera is operated in 
HDR mode. After bias-subtraction, the total signal (sky signal plus dark signal) of the high gain channel is 
divided by the total signal of the low-gain channel.
Input files are:
-   RAW flat images, taken at dusk with a script that acquires flat images of approximately equal exposure level. 
    The script attempts to make the median value of the image approximately 3000 ADU (the high-gain channel 
    saturates at 4095 ADU).
-   PFIT images, which comprise a 3D image with two planes. In the first, the pixel values are the dark current and 
    in the second they are the bias levels. The PFIT images enable interpolation to the exact exposure time of the 
    individual RAW flat images.
    
Output is a file with the suffix MFLAT. It has three HDUs:
-   The primary HDU is the flat-field image. It is a 3D image with two planes. The first (index = 0) is the high-gain 
    flat and the second (index = 1) is the low-gain flat.
-   The second HDU is the gain ratio image. The image name is 'GAIN RATIO'.
-   The third HDU is a table with information about the individual RAW images used to make the MFLAT image.
The table name is 'TABLE'.
The notebook batch-processes files within a directory according to keywords in user-configurable lists, 
one for gain and one for filter.
The output flat image and gain image both contain NANS that identify pixels determined using histograms to be 
outliers. The threshold criteria can be tuned by changing parameters in the code.
'''

from darepype.drp import DataFits # pipeline data object class
from darepype.drp.stepmiparent import StepMIParent # pipestep Multi-Input parent
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object class
###from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans # For masking/replacing

import astropy
#? from astropy import units as u
from astropy.io import fits #package to recognize FITS files

###import scipy.ndimage as nd
import numpy as np
import logging
###from skimage.measure import block_reduce

## BELOW: imports from Al's make_flat colab notebook 

 ## Variables used to control execution of the code.
#? camera = 'SBIG'    # Identifies camera being used. Currently the choices are SBIG and FLI (CMOS and CCD, respectively).
#? reduceby = 'make_flat_HDR_auto_51'  # Reduction software (this notebook)
 ## Core imports.
import sys
import os
print(sys.executable) ## Check to see if Jupyter is using the correct path.
 ### Generally useful mports and system path additions. May contain some imports not required for this notebook. 
from astropy.stats import mad_std           # The median absolute deviation, a more robust estimator than std.
from astropy.time import Time


############################################################################

class StepMasterFlatHdr(StepLoadAux, StepMIParent):
    '''Pipeline step object to produce and subtract master flats, pulls in pfits (3d image with one layer as the dark
    and one as the bias) and raw flats, requires valid config file'''
    
    def __init__(self):
            """ Constructor: Initialize data objects and variables
            """
            #Finian: taken from stepHDR ; "PFIT images, which comprise a 3D image with two planes. In the first, 
            #the pixel values are the dark current and in the second they are the bias levels. 
            #The PFIT images enable interpolation to the exact exposure time of the individual RAW flat images." 
            # WILL NEED HI GAIN AND LO GAIN FLAT, LIKEWISE FOR PFIT
            # Call superclass constructor (calls setup)
            
            super(StepMasterFlatHdr,self).__init__()

            # Pfit values
            self.hpfitloaded = False # indicates if bias has been loaded
            self.hpfit = None # Numpy array object containing bias values
            self.hpfitname = '' # name of selected bias file

            self.lpfitloaded = False # indicates if bias has been loaded
            self.lpfit = None # Numpy array object containing bias values
            self.lpfitname = '' # name of selected bias file

            # Flat values
            self.flatloaded = False # indicates if flat has been loaded
            self.flat = None # Numpy array object containing flat values
            self.flatname = '' # name of selected flat file

            # Finish up.
            self.log.debug('Init: done')

    def setup(self):
        # FINIAN: taken from stepmasterflat.py 
        ## replace load aux steps to load in pfit, that will be the only aux file
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
        
        # FINIAN: Al says I will need to change these parameters here
        self.paramlist.append(['combinemethod','median',
                               'Specifies how the files should be combined - options are median, average, sum'])
        self.paramlist.append(['outputfolder','',
                               'Output directory location - default is the folder of the input files'])
        # Get parameters for StepLoadAux, replace auxfile with pfit
        self.loadauxsetup('lpfit')
        self.loadauxsetup('hpfit')

        
    def timesortHDR(self, datalist, date_key = 'date-obs'):
        '''
          Sorts a list of fits files by a header keyword with a date/time value. In the case of an
          SBIG CMOS camera RAW file, the date/time is read from the second HDU.
          Arguments:
              filelist   = a list of fits files
              datapath   = the path to the files
              date_key   = the header keyword containing the time/date data
              print_list = if True, print the sorted file list
          Returns:
              tfiles     = the sorted file list
              utime      = a list of the unix times of the observations
          Author(s): Al Harper
          Modified: 210807, 210815
          Version: 1.1
        '''

        date_obs = []                                                # Make a list to hold the date-obs keyword strings.
        for d in datalist:
            if '_bin1L' in d.filename and '_RAW.' in d.filename:
                head = d.getheader(d.imgnames[1])                  # Get the header of the second HDU (index = [1]).
                date_obs.append(head[date_key])                      # Add date information to list. of string objects.
            else:
                head = d.getheader()                                # Get the header of the primary HDU (index = [0]).
                date_obs.append(head[date_key])                      # Add date information to list. of string objects.
        t = Time(date_obs, format='isot', scale='utc')               # Make an astropy time object in 'isot' format.  
        tsort = np.argsort(t)                                        # Make a list of indices that will sort by date_obs.
        tfiles = []
        utime = []
        for i in tsort:
            tfiles.append(datalist[i])
            utime.append(t[i].unix)
        return tfiles, utime
        
    def run(self):
        ## this bit taken from Al's colab notebook make_flat_HDR_auto_52_220803
        '''Loop to input and process RAW flatfiles'''
        
        '''Load and process a set of high and low gain flatfiles.'''
        ### Load pfit files
        self.hpfitname = self.loadauxname('hpfit', multi = False)
        self.hpfit = DataFits(config = self.config)
        self.hpfit.load(self.hpfitname)

        self.lpfitname = self.loadauxname('lpfit', multi = False)
        self.lpfit = DataFits(config = self.config)
        self.lpfit.load(self.lpfitname)
        
        #1st index low vs high, then intercepts and slopes of pfits
        
        polydarkimg = np.zeros((2,2,4096,4096))
        ## could use np.zeros_like(flatimage)
        print('Begin reduction.')

        '''Prepare a hot pixel mask'''
        print('')
        hotpixlim = 99.5                           # Input parameter for hot pixel limit
        print('hotpixlim =', hotpixlim)
        
        img = polydarkimg[0,0]
        rows, cols = polydarkimg.shape[2], polydarkimg.shape[3]
        uppercut = np.percentile(img, hotpixlim)
        print('uppercut =', uppercut)
        hotpix = np.where(img > uppercut)
        print('len hotpix =', len(hotpix[0]))
        print('')



        '''Use list comprehensions to select subsets of files to load. Then sort files by exposure time.'''

        ### HERE CREATE TWO LISTS, ONE FOR HIGH GAIN AND ONE FOR LOW GAIN, THEN TIME SORT EACH OF THESE AND PUT INTO 4D NUMPY STRUCTURE
        lowgainlist = []
        highgainlist = []
        for f in self.datain:
            if '.fit' in f.filename and '_bin1L' in f.filename:
                lowgainlist.append(f)
            elif '.fit' in f.filename and '_bin1H' in f.filename:
                highgainlist.append(f)
          
        highgainlist, utimeH = self.timesortHDR(highgainlist, date_key = 'date-obs')
        lowgainlist, utimeL = self.timesortHDR(lowgainlist, date_key = 'date-obs')

        '''Create output file name (assuming high and low gain flats are to be stored as
        a 2D image).'''
        
        lf = highgainlist[-1].filename.split('_')
        ff = highgainlist[0].filename.split('_')
        flatname = 'mflat_'+ff[1]+'_'+ff[3]+'DR'+'_'+ff[4]+'_'+ff[5]+'-'+lf[5]+'_'+ff[6]+'_'+ff[7]+'_'+'MFLAT.fits'
        print('flatname =',flatname)
        print('')

        '''Make a 4D image of high and low gain flat files. Dimensions are [gain, stack, row, col].'''

        flatnum = len(highgainlist)
        flatheadlist, flatstats = [[], []], [[], []]
        flatimage = np.zeros((2, flatnum, rows, cols))
        
        ## NEED TO CREATE NP ARRAYS TO HOLD STATISTICAL INFORMATION
        ## FILL THE ARRAYS WITHIN THESE LOOPS , LOOK AT TABLE TO MAKE SURE THEY'LL FILL IT CORRECTLY
        ## REPLACE ALL PRINTS WITH SELF.LOG.DEBUG
        for j in range(len(highgainlist)):
            flatimage[0, j] = highgainlist[j].image
            # Calculate some statistical information.
            mad[j] = mad_std(image[j],ignore_nan=True)
            median[j] = np.nanmedian(image[j])
            mean[j] = np.nanmean(image[j])
            std[j] = np.nanstd(image[j])
            
        
        for j in range(len(lowgainlist)):
            flatimage[1, j] = lowgainlist[j]
            # Calculate some statistical information.
            mad[j] = mad_std(image[j],ignore_nan=True)
            median[j] = np.nanmedian(image[j])
            mean[j] = np.nanmean(image[j])
            std[j] = np.nanstd(image[j])
           
                
  ### EXTRACT HEADER FROM HIGH GAIN LIST TO USE FOR THIS
        flatbaseheader = highgainlist[0].header.copy()

        ### NEED TO FIX INDENTATION BELOW THIS (minus four spaces)
        '''
            Print list of the files and median signal values.
            Compile a list of exposure times. Set variables for max and min exposure
            time for inclusion in header.
        '''

        print('Print filenames and image medians,')
        print('')
        for i in range(2):
            flatexptimes = []
            imedian = flatstats[i]['median']
            for j in range(len(flatfiles[i])):
                exptime = flatheadlist[i][j]['exptime']
                flatexptimes.append(exptime)
                print ('{:<4}{:<70}{:12.2f}'.format(j,flatfiles[i][j],imedian[j]))
            print('')
        xtimemin = np.min(flatexptimes)
        xtimemax = np.max(flatexptimes)
        print('xtimemin =', xtimemin, '  xtimemax =', xtimemax)
        print('')


        '''Construct a stack of dark images to match the flat image exposure times.'''

        darkimage = np.zeros_like(flatimage)
        for i in range(2):
            for j in range(flatimage.shape[1]):
                darkimage[i,j] = polydarkimg[i,1] + polydarkimg[i,0] * flatexptimes[j]


        '''Derive relative gain images for each HDR flat image in the stack. Note that
        we are measuring the total signal above bias resulting from both photons and
        dark current.'''

        biasimage = np.zeros((2, rows, cols))
        biasimage[0] = polydarkimg[0, 1]
        biasimage[1] = polydarkimg[1, 1]

        gainratiostack = np.zeros_like(flatimage[0])
        print('Medians of gain ratio images.')
        for i in range(len(flatfiles[0])):
            gainratiostack[i] = (flatimage[0][i] - biasimage[0]) / (flatimage[1][i] - biasimage[1])
            print(np.nanmedian(gainratiostack[i]))
        print('')

        '''
        Create median gain ratio image and std gain ratio image from the stack of gain 
        ratio images. Calculate and print some statistics of the two images, ignoring 
        masked values.
        '''

        # Take the median and std along the stack axis of gainratiostack.
        mediangainratioimg = np.nanmedian(gainratiostack, axis = 0)
        stdgainratioimg = np.nanstd(gainratiostack, axis = 0)

        # Take image medians and means of the two 2D images derived above.
        grat_median = np.nanmedian(mediangainratioimg)
        grat_mean = np.nanmean(mediangainratioimg)
        stdgrat_median = np.nanmedian(stdgainratioimg)
        stdgrat_mean = np.nanmean(stdgainratioimg)
        print('median, mean of median gain ratio image =', grat_median, grat_mean)
        print('max, min =', np.nanmax(mediangainratioimg), np.nanmin(mediangainratioimg))
        print('')
        print('median, mean of std gain ratio image =', stdgrat_median, stdgrat_mean)
        print('max, min =', np.nanmax(stdgainratioimg), np.nanmin(stdgainratioimg))
        print('')


        '''
        Create masked median gain ratio image and masked std image with nans replacing 
        masked pixels. Mask used here is hotpix (a hot pixel mask made from a dark exposure 
        in a preceding cell). Print some statistics, ignoring masked values.
        '''

        maskedgainratioimg = mediangainratioimg.copy()
        maskedgainratioimg[hotpix] = np.nan

        maskedstdgainratioimg = stdgainratioimg.copy()
        maskedstdgainratioimg[hotpix] = np.nan

        print('median, mean of masked gain ratio image =', np.nanmedian(maskedgainratioimg), 
            np.nanmean(maskedgainratioimg))
        print('max, min = of masked gain ratio image', np.nanmax(maskedgainratioimg), 
            np.nanmin(maskedgainratioimg))
        print('')
        print('median, mean of masked std gain ratio image =', np.nanmedian(maskedstdgainratioimg), 
            np.nanmean(maskedstdgainratioimg))
        print('max, min of masked std gain ratio image =', np.nanmax(maskedstdgainratioimg), 
            np.nanmin(maskedstdgainratioimg))
        print('')

######################################################################################################

#   PART THAT REPLACES HISTOGRAM CODE FOR GAIN MASK GOES HERE
#   WILL NEED TO REPLACE HOTPIX STUFF WITH HI AND LO GAIN LIMITS , DEFINED IN DCONF_FA
        img = maskedgainratioimg
        titlestring = ''
        masklow, maskhigh = grat_median - logainlim, grat_median + higainlim
        gainmask = np.where((img > masklow) & (img < maskhigh))
        
        #  AL: NOT SURE WHAT END PRODUCT OF THIS BLOCK SHOULD BE / IF I IMPLEMENTED THIS PROPERLY
      


            #print('')
            #hotpixlim = 99.5                           # Input parameter for hot pixel limit
            #print('hotpixlim =', hotpixlim)
            #img = polydarkimg[0,0]
            #uppercut = np.percentile(img, hotpixlim)
            #print('uppercut =', uppercut)
            #hotpix = np.where(img > uppercut)
            #print('len hotpix =', len(hotpix[0]))
            #print('')

#######################################################################################################


        '''Subtract interpolated darks from the flat images.'''

        if '_RAW.fit' in flatfiles[0][0]:
            flatimageDS = flatimage - darkimage
        #     print(flatimageDS.shape)

        '''
        Normalize each of the images in the high and low-gain flat image stack to its own median,
        and apply the gain-outlier mask derived above. From this point on, there will be nans
        for every pixel in both the hotpix and mgrmask masks.'''

        flatimageDSN = np.zeros_like(flatimageDS)
        flatmediansDS = np.zeros((2, flatnum))
        flatmadstdsDSN = np.zeros((2, flatnum))
        for i in range(2):
            for j in range(flatnum):
                flatimageDS[i,j][mgrmasklow] = np.nan
                flatimageDS[i,j][mgrmaskhigh] = np.nan
                flatmediansDS[i,j] = np.nanmedian(flatimageDS[i,j])
                flatimageDSN[i,j] = flatimageDS[i,j] / flatmediansDS[i,j]
                flatmadstdsDSN[i,j] = mad_std(flatimageDSN[i,j], ignore_nan=True)
        #     print(flatimageDSN.shape, flatmediansDS.shape)


        '''Make median high and low gain flats from the stacks and compute their respective
        median, mean, std, and mad. The result is a 3D image with high and low-gain planes.'''

        flat = np.nanmedian(flatimageDSN, axis=1)

        flatmedian, flatmean, flatstd, flatmadstd = [0,0], [0,0], [0,0], [0,0]

        for i in range(2):  
            flatmedian[i] = np.nanmedian(flat[i])                   # Median of the median flat.
            flatmean[i] = np.nanmean(flat[i])                       # Mean of the median flat.
            flatstd[i] = np.nanstd(flat[i])                         # mad_std of the median flat.
            flatmadstd[i] = mad_std(flat[i],ignore_nan=True)        # mad_std of the median flat.
        #     print('')
        print('Median flat median =', flatmedian )
        print('Mean flat median =', flatmean )
        print('Median flat std =', flatstd)
        print('Median flat mad_std =', flatmadstd)
        print('')


        '''Print a list of the medians of the dark-subtracted high and low gain flat images.'''
        print('List of medians of dark-subtracted high and low gain flat images')
        for i in range(flatimage.shape[1]):
            print('{:<3}{:<10.2f}{:<10.2f}'.format(i, flatmediansDS[0,i], flatmediansDS[1,i]))
        print('')

        '''Find minimum and maximum values of the image medians for inclusion in output header.'''

        medmin, medmax = [0,0], [0,0]
        medmin[0] = np.nanmin(flatmediansDS[0])
        medmax[0] = np.nanmax(flatmediansDS[0])
        medmin[1] = np.nanmin(flatmediansDS[1])
        medmax[1] = np.nanmax(flatmediansDS[1])
        print('minimum medians =', medmin, '   ', 'maximum medians =', medmax )
        print('')

        '''
        Re-normalize the median flats (H and L) to their medians and mask out pixels with 
        gains greater than or less than 1.0 from the median gain. Replace the masked pixels 
        with np.nan. Hence, when a sky image is divided by the flat, those pixels will also 
        be masked (will be np.nan) in the flat-fielded sky image (in addition to the pixela
        aleady replaced with nans using the hotpix mask).
        '''

        mflat = np.zeros_like(flat)
        mflatmedian = np.zeros((2))
        mflatmadstd = np.zeros((2))
        for i in range(2):
            mflat[i] = flat[i] / flatmedian[i]
            mflatmedian[i] = np.nanmedian(mflat[i])
            mflat[i][mgrmasklow] = np.nan
            mflat[i][mgrmaskhigh] = np.nan
            mflatmadstd[i] = mad_std(mflat[i], ignore_nan=True)
        print('mflatmedian, mflatmadstd =', mflatmedian, mflatmadstd)
        print('')

        '''Free up memory being used by in-line images.'''



        '''
        Compute the differences of each of the individual normalized flat images from the
        median flat of the stack. Compute and print the median, std, and mad_std of
        each of the difference images. If display = True, also show the images.
        '''

        g = 0           # Gain = 0 for high-gain, 1 for low-gain.

        display = show_difimages  # If True, display the individual difference images.

        vmx = mflatmadstd[g] * 2.0
        vmn = - mflatmadstd[g] * 2.0

        difimage = np.zeros_like(flatimageDS)  # Ratio of inidvidual flat images to the median flat image.
        difmadstd = np.zeros((flatnum))
        difstd = np.zeros((flatnum))
        difmean = np.zeros((flatnum))
        difmedian = np.zeros((flatnum))


        print('List difimage median, mean, std, and mad_std')
        for i in range(flatimageDS.shape[1]):
            difimage[g,i] = flatimageDSN[g,i]-mflat[g]
            difmadstd[i] = mad_std(difimage[g,i],ignore_nan=True)
            difstd[i] = np.nanstd(difimage[g,i])
            difmean[i] = np.nanmean(difimage[g,i])
            difmedian[i] = np.nanmedian(difimage[g,i])
            print(i, difmedian[i], difmean[i], difstd[i], difmadstd[i])  


        '''Create output DataFits object and fill with image data.'''

        dataout = DataFits(config=config)
        dataout.header = flatbaseheader.copy()
        dataout.image = mflat
        dataout.imageset(maskedgainratioimg, 'GAIN RATIO')


        '''
        Create numpy arrays with information from headers, make_stack, and timesortDF functions.
        '''

        imedianH, imadH, imeanH, istdH = flatstats[0]['median'], flatstats[0]['mad']\
                                        , flatstats[0]['mean'], flatstats[0]['std']     # From make_stackDF
        imedianL, imadL, imeanL, istdL = flatstats[1]['median'], flatstats[1]['mad']\
                                        , flatstats[1]['mean'], flatstats[1]['std']     # From make_stackDF

        utime = np.asarray(utimeH)  # From time_sortDF.
        etime = utime - utime[0]    # Time elepased from the beginning of the first exposure of the sequence.
        index = np.arange(flatnum)  # Make a column with the file sequence numbers 

        # From header:
        ambient = np.zeros((flatnum))
        primary = np.zeros((flatnum))
        secondar = np.zeros((flatnum))
        dewtem1 = np.zeros((flatnum))
        for i in range(flatnum):
            ambient[i] = flatheadlist[0][i]['ambient']
            primary[i] = flatheadlist[0][i]['primary']
            secondar[i] = flatheadlist[0][i]['secondar']
            dewtem1[i] = flatheadlist[0][i]['dewtem1']

        '''Now put put derived and header data into a fits table and add it to the output object.'''

        # Make file identifiers:
        IDs = []
        for i in range(flatnum):
            fi = flatfiles[0][i].split('_')
            IDs.append(fi[4]+'_'+fi[5])
        fileIDs = np.asarray(IDs)

        # Make a list of fits column objects.
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

        # Make table
        c = fits.ColDefs(tcols)
        table = fits.BinTableHDU.from_columns(c)
        tabhead = table.header

        # Add table to the datafits output object.
        dataout.tableset(table.data, tablename = 'table', tableheader=tabhead)

        '''Populate output header with new keyword data.'''

        dataout.header['notes'] = 'High and low gain master flats (3D image)'
        dataout.header['notes2'] = '2nd HDU is a gain ratio image. Contains nans'
        dataout.header['notes3'] = 'Table contains statistical and environmental data'
        dataout.header['filelist'] = filestring
        dataout.header['oscnmean'] = 0.0
        dataout.header['imagetyp'] = 'HDR MFLAT'
        dataout.setheadval('xtimemin', xtimemin, 'Minimum exposure in the set of flats')
        dataout.setheadval('xtimemax', xtimemax, 'Maximum exposure in the set of flats')
        dataout.setheadval('medminH', medmin[0], 'Minimum high gain median in the set of flats')
        dataout.setheadval('medmaxH', medmax[0], 'Maximum high gain median in the set of flats')
        dataout.setheadval('medminL', medmin[1], 'Minimum low gain median in the set of flats')
        dataout.setheadval('medmaxL', medmax[1], 'Maximum low gain median in the set of flats')
        dataout.setheadval('reduceby', reduceby, 'Reduction software')
        dataout.header['bzero'] = 0.0
        dataout.header['numfiles'] = flatnum
        dataout.header['ambient'] = np.nanmean(ambient)
        dataout.header['primary'] = np.nanmean(primary)
        dataout.header['secondar'] = np.nanmean(secondar)
        dataout.header['dewtem1'] = np.nanmean(dewtem1)
        dataout.header['meandstd'] = np.nanmean(dstd)
            ## might also be useful to have a max dstd column, if I use setheadval I can add a comment line
            ## .header and .setheadval are somewhat interchangeable


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
