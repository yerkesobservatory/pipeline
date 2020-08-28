#!/usr/bin/env python
""" PIPE RGB - Version 1.1.0

    Code for StepRGB in pipeline: combines filters to colored image & adds a label
    
    @author: M. Berthoud and A. Pal
"""

import os # os library
import numpy # numpy library
import logging # logging object library
import pylab # pylab library for creating rgb image
import time # time library for dated folders
from astropy.visualization import simple_norm
from PIL import Image # image library for saving rgb file as JPEG
from PIL import ImageFont # Libraries for adding a label to the color image
from PIL import ImageDraw
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepMOParent # pipe step parent object

class StepRGB(StepMOParent):
    """ Stone Edge Pipeline Step RGB Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version

    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
	# call superclass constructor (calls setup)
        super(StepRGB,self).__init__()
	# list of data
        self.datalist = [] # used in run() for every new input data file
	# set configuration
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
        ### Set Names
        # Name of the pipeline reduction step
        self.name='makergb'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'rgb'
        # Set Logger for this pipe step
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['minpercent', 0.05, 
		'Specifies the percentile for the minimum scaling'])
        self.paramlist.append(['maxpercent', 0.999,
		'Specifies the percentile for the maximum scaling'])
        self.paramlist.append(['filterorder', 'sii|i-band|h-alpha|r-band|clear|oiii|g-band',
        'Specifies the filters to be used in this step from reddest to bluest. ' + 
        'Filters should be separated by the pipe (|) character. ' + 
        'Example: "i-band|r-band|g-band"'])
        self.paramlist.append(['filterprefs', 'i-band|r-band|g-band',
        'Specifies the preferred filters for the R, G, and B channels. ' +
        'Filters should be separated by the pipe (|) character. ' + 
        'Example: "i-band|r-band|g-band"'])
        self.paramlist.append(['folderpaths', '',
        'Specifies optional folders to copy the finished .jpeg file to. ' +
        'Paths should be separated by the colon (:) character and can be ' +
        'time-formatted. ' +
        r'Example: "~/myfiles:/data/images/StoneEdge/0.5meter/%Y/%Y-%m-%d'])
        self.paramlist.append(['createfolders', False,
        'Specifies whether nonexistent paths given as save folders should be created'])

    def run(self):
        """ Runs the combining algorithm. The self.datain is run
            through the code, the result is in jpeg_dataout.
        """
        ''' Select 3 input dataset to use, store in datause '''
        #Store number of inputs
        num_inputs = len(self.datain)
        # Create variable to hold input files
        # Copy input to output header and filename
        datause = [None, None, None]
        self.log.debug('Number of input files = %d' % num_inputs)
        
        if num_inputs == 0:   # Raise exception for no input
            raise ValueError('No input')
        elif num_inputs == 1:
            datause = [self.datain[0], self.datain[0], self.datain[0]]
        elif num_inputs == 2:
            datause = [self.datain[0], self.datain[0], self.datain[1]]
        else:
            filterorder_list = self.getarg('filterorder').split('|')
            filterprefs_list = self.getarg('filterprefs').split('|')

            datain_filter_list = [element.getheadval('filter') for element in self.datain]
            used_filter_flags = [False] * len(self.datain)


            if len(filterprefs_list) != 3:
                self.log.error('Invalid number of preferred filters provided (should be 3): ' + 
                               self.getarg('filterprefs'))
            else:
                # Locate data matching the filters specified in filterprefs
                for i, preferred_filter in enumerate(filterprefs_list):
                    for j, element in enumerate(self.datain):
                        if element.getheadval('filter') == preferred_filter:
                            datause[i] = element
                            used_filter_flags[j] = True
                            break

            filterorder_walker = 0
            for i, channel in enumerate(datause):
                if channel == None:
                    for ordered_filter in filterorder_list[filterorder_walker:]:
                        filterorder_walker = filterorder_walker + 1
                        if ordered_filter in datain_filter_list:
                            datain_index = datain_filter_list.index(ordered_filter)
                            if not used_filter_flags[datain_index]:
                                datause[i] = self.datain[datain_index]
                                used_filter_flags[datain_index] = True
                                break
                elif channel.getheadval('filter') in filterorder_list:
                    filterorder_walker = filterorder_list.index(channel.getheadval('filter'))
                    
            for i, channel in enumerate(datause):
                if channel == None:
                    for j, datain_filter in enumerate(datain_filter_list):
                        if not used_filter_flags[j]:
                            datause[i] = self.datain[j]
                            used_filter_flags[j] = True
                            break
                    
        self.log.debug('Files used: R = %s  G = %s  B = %s' % (datause[0].filename, datause[1].filename, datause[2].filename) )
        jpeg_dataout = DataFits(config = self.config)
        jpeg_dataout.header = datause[0].header
        jpeg_dataout.filename = datause[0].filename
        img = datause[0].image
        img1 = datause[1].image
        img2 = datause[2].image
        
        ''' Finding Min/Max scaling values '''
        # Create a Data Cube with floats
        datacube = numpy.zeros((img.shape[0], img.shape[1], 3), dtype=float)
        # Enter the image data into the cube so an absolute max can be found
        datacube[:,:,0] = img
        datacube[:,:,1] = img1
        datacube[:,:,2] = img2
        # Find how many data points are in the data cube
        datalength = img.shape[0] * img.shape[1] * 3
        # Create a 1-dimensional array with all the data, then sort it	
        datacube.shape=(datalength,)
        datacube.sort()
        # Now use arrays for each filter to find separate min values
        rarray = img.copy()
        garray = img1.copy()
        barray = img2.copy()
        # Shape and sort the arrays
        arrlength = img.shape[0] * img.shape[1]
        rarray.shape=(arrlength,)
        rarray.sort()
        garray.shape=(arrlength,)
        garray.sort()
        barray.shape=(arrlength,)
        barray.sort()
        # Find the min/max percentile values in the data for scaling
        # Values are determined by parameters in the pipe configuration file
        minpercent = int(arrlength * self.getarg('minpercent'))
        maxpercent = int(datalength * self.getarg('maxpercent'))
        # Find the final data values to use for scaling from the image data
        rminsv = rarray[minpercent]  #sv stands for "scalevalue"
        gminsv = garray[minpercent]
        bminsv = barray[minpercent]
        maxsv = datacube[maxpercent]
        self.log.info(' Scale min r/g/b: %f/%f/%f' % (rminsv,gminsv,bminsv))
        self.log.info(' Scale max: %f' % maxsv)
        # The same min/max values will be used to scale all filters
        ''' Finished Finding scaling values	'''

        ''' Combining Function '''
        # Make new cube with the proper data type for color images (uint8)
        # Use square root (sqrt) scaling for each filter
        # log or asinh scaling is also available
        #astropy.vidualizations.SqrtStretch()
        imgcube = numpy.zeros((img.shape[0], img.shape[1], 3), dtype='uint8')
        minsv = [rminsv, gminsv, bminsv]
        for i in range(3):
            # Make normalization function
            norm = simple_norm(datause[i].image, 'sqrt', min_cut = minsv[i], max_cut = maxsv)
            # Apply it
            imgcube[:,:,i] = norm(datause[i].image) * 255.
        jpeg_dataout.image = imgcube
        # Create variable containing all the scaled image data
        imgcolor = Image.fromarray(jpeg_dataout.image, mode='RGB')
        # Save colored image as a .tif file (without the labels)
        imgcolortif = imgcube.copy()
        imgcolortif.astype('uint16')
        ### tiff.imsave('%s.tif' % jpeg_dataout.filenamebase, imgcolortif)
        ''' End of combining function '''
        
        ''' Add a Label to the Image '''
        draw = ImageDraw.Draw(imgcolor)
        # Use a variable to make the positions and size of text relative
        imgwidth = img.shape[1]
        imgheight = img.shape[0]
        # Open Sans-Serif Font with a size relative to the picture size
        try:
            # This should work on Linux
            font = ImageFont.truetype('/usr/share/fonts/liberation/LiberationSans-Regular.ttf',imgheight//41)
        except:
            try:
                # This should work on Mac
                font = ImageFont.truetype('/Library/Fonts/Arial Unicode.ttf',imgheight//41)
            except:
                try:
                    # This should work on Windows
                    font = ImageFont.truetype('C:\\Windows\\Fonts\\arial.ttf',imgheight//41)
                except:
                    # This should work in Colab
                    font = ImageFont.truetype('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf',imgheight//41)
                    # If this still doesn't work - then add more code to make it run on YOUR system
        # Use the beginning of the FITS filename as the object name
        filename = os.path.split(jpeg_dataout.filename)[-1]
        try:
            objectname = filename.split('_')[0]
            objectname = objectname[0].upper()+objectname[1:]
        except Exception:
            objectname = 'Unknown.'
        objectname = 'Object:  %s' % objectname
        # Read labels at their respective position (kept relative to image size)
        # Left corner: object, observer, observatory
        # Right corner: Filters used for red, green, and blue colors
        draw.text((imgwidth/100,imgheight/1.114), objectname, (255,255,255), font=font)
        # Read FITS keywords for the observer, observatory, and filters
        if 'OBSERVER' in jpeg_dataout.header:
            observer = 'Observer:  %s' % jpeg_dataout.getheadval('OBSERVER')
            draw.text((imgwidth/100,imgheight/1.073), observer, (255,255,255), font=font)
        if 'OBSERVAT' in jpeg_dataout.header:
            observatory = 'Observatory:  %s' % jpeg_dataout.getheadval('OBSERVAT')
            draw.text((imgwidth/100,imgheight/1.035), observatory, (255,255,255), font=font)
        if 'FILTER' in datause[0].header:
            red = 'R:  %s' % datause[0].getheadval('FILTER')
            draw.text((imgwidth/1.15,imgheight/1.114),red, (255,255,255), font=font)
        if 'FILTER' in datause[1].header:
            green = 'G:  %s' % datause[1].getheadval('FILTER')
            draw.text((imgwidth/1.15,imgheight/1.073),green, (255,255,255), font=font)
        if 'FILTER' in datause[2].header:
            blue = 'B:  %s' % datause[2].getheadval('FILTER')
            draw.text((imgwidth/1.15,imgheight/1.035),blue, (255,255,255), font=font)

        # Make image name
        imgname = jpeg_dataout.filenamebegin
        if imgname[-1] in '_-,.': imgname=imgname[:-1]
        imgname +='.jpg'
        # Save the completed image
        imgcolor.save(imgname)
        self.log.info('Saving file %sjpg' %jpeg_dataout.filenamebegin)

        # Optional folder output setup
        baseimgname = os.path.basename(imgname)
        folderpaths_list = self.getarg('folderpaths').split(':')
        for path in folderpaths_list:
            path = time.strftime(path, time.localtime())
            if not os.path.exists(path):
                if self.getarg('createfolders'):
                    os.makedirs(path)
                    self.log.info('Creating directory %s' %path)
                else:
                    self.log.info('Invalid folder path %s' %path)
            try:
                imgcolor.save(os.path.join(path, baseimgname))
            except:
                self.log.exception('Could not save image to directory %s' %path)

        
        ''' End of Label Code '''
        # Set complete flag
        jpeg_dataout.setheadval('COMPLETE',1,
                                    'Data Reduction Pipe: Complete Data Flag')

        ### Make output data
        self.dataout = self.datain.copy()
        self.dataout.append(jpeg_dataout)
        
        
    
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.log.debug('Reset: done')
        
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step rgb')

        # log message
        self.log.info('Testing pipe step rgb - Done')
    
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
    StepRGB().execute()

#!/usr/bin/env python


''' This is a version of a script that executes the Data Reduction Pipeline
    from a folder one level above the images. The final goal is to have a script 
    connected to Yerkes' File Manager server that will automatically execute 
    whenever new images are uploaded, or on a set time interval (whichever is 
    more practical).

    This version of the pipeline is optimized for use with i/r/g-band images.
    It will still work if they are not there, but the order of images run through
    the pipeline may be incorrect.  If a different combination of filters are
    being used, the code can be edited appropriately (this version contains
    'dummy' code that is currently commented out that can be easily
    edited/uncommented to work with any set of desired filters).
'''



""" === History ===
    2014-06-30 New file created by Neil Stilin from template file by Nicolas Chapman
    2014-07-09 Main code for creating RBG image added by Neil Stilin
    2014-07-29 Code has been improved by adding better scaling and image labels
    2014-08-06 Added 'if' functions to the label printing so that if keywords do not exist in the header(s), they are skipped rather than raising an error --NS
    2014-08-11 This file was essentially just renamed. The file called steprgb.py now uses raw inputs to determine the scaling values.  --NS
"""
