#!/usr/bin/env python
""" PIPE CLEAN - Version 1.0.0

    This module defines the HAWC pipeline step parent object. Pipe steps are
    the modules responsible for all HAWC data reduction. They are called by
    the pipeline and work with pipedata objects. All pipe step objects are
    descendants from this one. Pipe steps are callable objects that return
    the reduced data product (as pipedata object).
    
    @author: berthoud
"""

import os # os library
import numpy # numpy library
import logging # logging object library
import pyfits # pyfits library (for accessing header data) @UnresolvedImport
import scipy # scipy library
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepParent # pipe step parent object

class StepDark(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepDark,self).__init__()
        # list of data and darks
        self.datalist = [] # used in run() for every new input data file
        # dark values
        self.darkloaded = 0 # indicates if dark has been loaded
        self.darks = [] # list containing arrays with dark values
        self.goodpixmap = numpy.zeros([1,1], dtype=numpy.int16) # array with good pixels
        # dark file info and header keywords to fit
        self.darkfile = '' # name of selected dark file
        self.darkhead = pyfits.PrimaryHDU(numpy.array(1)) # header of dark file
        self.fitkeys = [] # FITS keywords that have to fit
        self.keyvalues = [] # values of the keywords (from the first data file)
        # set configuration
        self.log.debug('Init: done')
    
    def setup(self):
        """ ### Names and Prameters need to be Set Here ###
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
        self.name='dark'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'drk'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['darkfile', 'search',
            'Filename for dark file or "search" for searching ' +
            'file in darkfolder (default = search)'])
        self.paramlist.append(['fitkeys', ['DETSIZE', 'PASSBAND', 'PUPIL'],
            'List of keys that need to match dark and data file ' +
            '(default = DETSIZE PASSBAND PUPIL) ' +
            '- only used if darkfile=search'])
        self.paramlist.append(['darkfolder', '.',
            'Folder for searching dark files ' +
            '(default = . ) - only used if darkfile=search'])
        self.paramlist.append(['l0method', 'NO',
            'Method to normalize data: NOne, REal, IMag and ABSolute ' +
            '(default = NO)'])
        self.paramlist.append(['datalist', [],
            'List of data sets to subtract dark in input file ' +
            '(default = [] i.e. only subtract dark from image cube in first HDU)'])

    def run(self):
        """ Runs the dark subtraction step. The self.datain is run
            through the code, the result is in self.dataout.
        """
	''' Preparation (finding dark file) '''
        # Load dark file if necessary
        if not self.darkloaded:
            self.loaddark()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.fitkeys)):
                if self.keyvalues[keyind] != self.datain.getheadval(self.fitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.fitkeys[keyind])
	# Copy incoming image data
	self.dataout = self.datain.copy()
	img = self.datain.image
	''' Dark Subtraction code '''
	darkimg = self.darks[0]  # Not sure if this is correct
	img = img - darkimg
	''' Dark Subtraction code (END)'''
	self.dataout.image = img
        # Set complete flag
        self.dataout.setheadval('COMPLETE',1,
                                'Data Reduction Pipe: Complete Data Flag')

    def loaddark(self):
        """ Loads the dark information for the instrument settings
            described in the header of self.datain.
            
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        ### identify dark file to load, search if requested
        darkfile = self.getarg('darkfile')
        if darkfile == 'search' :
            # get list of keywords to fit
            fitkeys = self.getarg('fitkeys')
            # check format (make first element uppercase)
            try:
                _ = fitkeys[0].upper()
            except AttributeError:
                # AttributeError if it's not a string
                self.log.error('LoadDark: fitkeys config parameter is ' +
                               'incorrect format')
                raise TypeError('fitkeys config parameter is incorrect format')
            # get keywords from data
            datakeys=[]
            for fitkey in fitkeys:
                datakeys.append(self.datain.getheadval(fitkey))
            # get dark files from darkdir folder
            darkfolder = self.getarg('darkfolder')
            filelist=[name for name in os.listdir(darkfolder)
                      if name[0] != '.' and name.find('.fit') > -1 ]
            if len(filelist) < 1:
                self.log.error('LoadDark: no dark files found in folder ' +
                               darkfolder)
                raise ValueError('no dark files found in folder ' +
                                 darkfolder)
            # match dark files, return best dark file
            bestind = 0 # index of file with best match in filelist
            bestfitn = 0 # number of keywords that match in best match
            fileind = 0 # index for going through the list
            while fileind < len(filelist) and bestfitn < len(fitkeys):
                # load keys of dark file
                filehead = pyfits.getheader(darkfolder+'/'+filelist[fileind])
                filekeys=[]
                for fitkey in fitkeys:
                    try:
                        filekeys.append(filehead[fitkey])
                    except KeyError:
                        self.log.warning('LoadDark: missing key [%s] in dark <%s>'
                                       % (fitkey, filelist[fileind] ) )
                        filekeys.append('')
                # determine number of fitting keywords
                keyfitn=0
                while ( keyfitn < len(fitkeys) and 
                        datakeys[keyfitn] == filekeys[keyfitn] ):
                    keyfitn = keyfitn + 1
                # compare with previous best find
                if keyfitn > bestfitn:
                    bestind = fileind
                    bestfitn = keyfitn
                fileind=fileind+1
            darkfile = darkfolder+'/'+filelist[ bestind ]
            if bestfitn < len(fitkeys):
                self.log.warn('Could not find perfect dark file match')
                self.log.warn('Best dark file found is <%s>'
                              % filelist[bestind] )
            else:
                self.log.info('Best dark file found is <%s>'
                              % filelist[bestind] )
            self.fitkeys = fitkeys
            self.keyvalues = datakeys            
        ### load dark data into a DataFits object
        self.darkfile = darkfile
        darkdata = DataFits(config = self.config)
        darkdata.load(self.darkfile)
        ### find dark image data arrays and store them
        # get sizes of input data
        datalist = self.getarg('datalist')
        if len(datalist) > 0:
            # There are items in datalist -> loop over items
            self.darks = []
            # Check if necessary number of images in darkdata
            if len(darkdata.imgdata) < len(datalist): 
                msg = 'Number of images in dark file < '
                msg += 'number of entries in datalist'
                self.log.error('LoadDark: %s' % msg)
                raise ValueError(msg)
            # Loop through datalist items
            for dataind in range(len(datalist)):
                dataitem = datalist[dataind]
                # Search for dataitem in self.datain images
                if dataitem.upper() in self.datain.imgnames:
                    dataimg = self.datain.imageget(dataitem)
                    self.log.debug('LoadDark: Found image <%s> to subtract dark'
                                   % dataitem)
                # Search dataitem in self.table columns
                else:
                    try:
                        dataimg = self.datain.table[dataitem]
                        self.log.debug('LoadDark: Found column <%s> to subtract dark'
                                       % dataitem)
                    except:
                        msg = 'No data found for <%s>' % dataitem
                        self.log.error('LoadDark: %s' % msg)
                        raise ValueError(msg)
                # Get dimensions - append dark to list
                datasiz = dataimg.shape
                if self.getarg('l0method').upper() != 'NO':
                    datasiz = datasiz[1:]
                darksiz = darkdata.imgdata[dataind].shape
                self.darks.append(darkdata.imgdata[dataind])
                # Check dimension with dark data
                print(datasiz,darksiz)
                if len(datasiz) >= len(darksiz):
                    # Data has >= dimensions than dark -> compare
                    begind = len(datasiz)-len(darksiz)
                    if datasiz[begind:] != darksiz:
                        msg = 'Dark "%s" does not fit data - A' % dataitem
                        self.log.error('LoadDark: %s' % msg)
                        raise ValueError(msg)
                else:
                    # More dimensions in dark data -> report error
                    msg = 'Dark "%s" does not fit data - B' % dataitem
                    self.log.error('LoadDark: %s' % dataitem)
                    raise ValueError(msg)
        else:
            # Empty datalist -> Subtract dark from first image in data with first dark
            datasiz = self.datain.image.shape
            if self.getarg('l0method').upper() != 'NO':
                datasiz = datasiz[1:]
            darksiz = darkdata.image.shape
            if len(datasiz) >= len(darksiz):
                # Data has >= dimensions than dark -> compare
                begind = len(datasiz)-len(darksiz)
                if datasiz[begind:] != darksiz:
                    self.log.error('LoadDark: Dark does not fit data - A')
                    raise ValueError('Dark does not fit data - A')
            else:
                # More dimensions in dark data -> report error
                self.log.error('LoadDark: Dark does not fit data - B')
                raise ValueError('Dark does not fit data - B')
            self.log.debug('LoadDark: Subtracting Dark from first data image with first dark')
            self.darks=[darkdata.image]
        ### make good pixel map for each detector and add to 
        #darktemp = numpy.abs(data[0,...])+numpy.abs(data[1,...])
        #self.goodpixmap = numpy.ones(data.shape[1:])
        #self.goodpixmap [ numpy.where(darktemp == 0.0)] = 0.0
        # Finish up
        self.darkloaded = 1
        self.log.debug('LoadDark: done')
    
    def undo(self):
        """ Reverts the pipeline step to its status before the last call
        """
        self.log.debug('Undo: done')
    
    def reset(self):
        """ Resets the step to the same condition it was when it was
            created. Stored dark image data is erased and the configuration
            information is cleared.
        """
        self.darkloaded = 0
        self.darkvalue = numpy.zeros([1,1])
        self.darkphase = numpy.zeros([1,1])
        self.darkheader = pyfits.PrimaryHDU(numpy.array(1))
        self.darkfile = ''
        
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)

        # log message
        self.log.info('Testing pipe step %s - Done' %self.name)
    
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
    StepDark().execute()

""" === History ===
	2014-12-29 -- First basic version created by Neil Stilin
	2014-12-30 -- Copied/Edited code for loading dark files from StepFlat (made by Marc Berthoud)
"""
