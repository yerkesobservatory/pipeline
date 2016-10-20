""" PIPE FLAT - Version 0.1.0

    This module defines the HAWC pipeline flatfielding step object.
    The object is a child of stepparent. Flat loads the appropriate
    flat file.
    
    In the input data, Flat looks for data objects to be flatted named in 
    datalist and flats them with the images in the flatfile. The objects
    in datalist are flatted with the first images in the flat file. If
    datalist is empty (None), the first image of the input data is flatted
    using the first image of the flatfile. Alternatively if no images are
    found in the data the Flat looks for the data as columns of the
    first table. In that case the data is moved from the table to images.
    
    Flat also accepts input images that contain real and
    imaginary components of the demodulated data. In that case l0method
    is used to normalize to real values.
    
    The selection of the flat file is done my matching specific keywords in
    the header of the input data. If no such file can be found, the best fit
    is selected. For each subsequent data set the keywords in the header
    are checked against the ones from the flat file.
    
    Task List (changes to do)
    * loadflat: still make good pixel map (now list with one for each
      detector) - WILL BE DONE LATER (probably need badpixel map and
      several codes for bad pixels)
      * Make badpixmap for bad pixels
"""

import os # os library
import numpy # numpy library
import logging # logging object library
from astropy.io import fits # pyfits library (for accessing header data) @UnresolvedImport
from drp.pipedata import PipeData # pipeline data object
from drp.stepparent import StepParent # pipe step parent object

class StepFlat(StepParent):
    """ HAWC Pipeline Flatfielding Step Object
    """
    
    stepver = '0.1' # pipe step version

    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepFlat,self).__init__()
        # list of data and flats
        self.datalist = [] # used in run() for every new input data file
        # flat values
        self.flatloaded = 0 # indicates if flat has been loaded
        self.flats = [] # list containing arrays with flat values
        self.goodpixmap = numpy.zeros([1,1], dtype=numpy.int16) # array with good pixels
        # flat file info and header keywords to fit
        self.flatfile = '' # name of selected flat file
        self.flathead = fits.PrimaryHDU(numpy.array(1)) # header of flat file
        self.fitkeys = [] # FITS keywords that have to fit
        self.keyvalues = [] # values of the keywords (from the first data file)
        # set configuration
        self.log.debug('Init: done')
        
    def setup(self):
        """ ### Names and Prameters need to be Set Here ###
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
        self.name='flat'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'fla'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['flatfile', 'search',
            'Filename for flat file or "search" for searching ' +
            'file in flatfolder (default = search)'])
        self.paramlist.append(['fitkeys', ['DETSIZE', 'PASSBAND', 'PUPIL'],
            'List of keys that need to match flat and data file ' +
            '(default = DETSIZE PASSBAND PUPIL) ' +
            '- only used if flatfile=search'])
        self.paramlist.append(['flatfolder', '.',
            'Folder for searching flat files ' +
            '(default = . ) - only used if flatfile=search'])
        self.paramlist.append(['l0method', 'NO',
            'Method to normalize data: NOne, REal, IMag and ABSolute ' +
            '(default = NO)'])
        self.paramlist.append(['datalist', [],
            'List of data sets to flatten in intput file ' +
            '(default = [] i.e. only flatten image cube in first HDU)'])

    def run(self):
        """ Runs the flatfielding algorithm. The flatfielded data is
            returned in self.dataout
        """
        ### Preparation
        # Load flat files if necessary
        if not self.flatloaded:
            self.loadflat()
        # Else: check data for correct instrument configuration
        else:
            for keyind in range(len(self.fitkeys)):
                if self.keyvalues[keyind] != self.datain.getheadval(self.fitkeys[keyind]):
                    self.log.warn('New data has different FITS key value for keyword %s' %
                                  self.fitkeys[keyind])
        ### Copy datain to dataout
        self.dataout = self.datain.copy()
        ### Apply Flatfield
        # Only one data set -> flatfield it
        datalist = self.getarg('datalist')
        if len(datalist) == 0:
            # Get Image
            image = self.datain.image.copy()
            # Flatfield it
            self.dataout.image = self.flatfield(image,self.flats[0])
        # Loop through data sets
        else:
            for dataind in range(len(datalist)):
                dataitem = datalist[dataind]
                # Search for dataitem in self.datain images -> flatfield it
                if dataitem.upper() in self.datain.imgnames:
                    image = self.datain.imageget(dataitem)
                    image = self.flatfield(image, self.flats[dataind])
                    self.dataout.imageset(image, dataitem)
                    continue # go to next dataitem (skip end of loop)
                # Search for dataitem in self.table columns
                try:
                    image = self.datain.table[dataitem]
                except:
                    msg = 'No data found for %s in file %s' % (dataitem, self.datain.filename)
                    self.log.error('Run: %s' % msg)
                    raise ValueError(msg)
                # Flatfield
                image = self.flatfield(image, self.flats[dataind])
                # Store image in dataout
                self.dataout.imageset(image, dataitem)
                # Remove table column from dataout
                self.dataout.tabledelcol(dataitem)
        # Add goodpixmap to output data
        self.dataout.imageset(self.goodpixmap, 'goodpixmap')
        ### Finish - cleanup
        # Set complete flag
        self.dataout.setheadval('COMPLETE',1,
                           'Data Reduction Pipe: Complete Data Flag')
    
    def loadflat(self):
        """ Loads the flat information for the instrument settings
            described in the header of self.datain.
            
            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        ### identify flat file to load, search if requested
        flatfile = self.getarg('flatfile')
        if flatfile == 'search' :
            # get list of keywords to fit
            fitkeys = self.getarg('fitkeys')
            # check format (make first element uppercase)
            try:
                _ = fitkeys[0].upper()
            except AttributeError:
                # AttributeError if it's not a string
                self.log.error('LoadFlat: fitkeys config parameter is ' +
                               'incorrect format')
                raise TypeError('fitkeys config parameter is incorrect format')
            # get keywords from data
            datakeys=[]
            for fitkey in fitkeys:
                datakeys.append(self.datain.getheadval(fitkey))
            # get flat files from flatdir folder
            flatfolder = self.getarg('flatfolder')
            filelist=[name for name in os.listdir(flatfolder)
                      if name[0] != '.' and name.find('.fit') > -1 ]
            if len(filelist) < 1:
                self.log.error('LoadFlat: no flat files found in folder ' +
                               flatfolder)
                raise ValueError('no flat files found in folder ' +
                                 flatfolder)
            # match flat files, return best flat file
            bestind = 0 # index of file with best match in filelist
            bestfitn = 0 # number of keywords that match in best match
            fileind = 0 # index for going through the list
            while fileind < len(filelist) and bestfitn < len(fitkeys):
                # load keys of flat file
                filehead = fits.getheader(flatfolder+'/'+filelist[fileind])
                filekeys=[]
                for fitkey in fitkeys:
                    try:
                        filekeys.append(filehead[fitkey])
                    except KeyError:
                        self.log.warning('LoadFlat: missing key [%s] in flatfile <%s>'
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
            flatfile = flatfolder+'/'+filelist[ bestind ]
            if bestfitn < len(fitkeys):
                self.log.warn('Could not find perfect flat file match')
                self.log.warn('Best flat file found is <%s>'
                              % filelist[bestind] )
            else:
                self.log.info('Best flat file found is <%s>'
                              % filelist[bestind] )
            self.fitkeys = fitkeys
            self.keyvalues = datakeys            
        ### load flat data into a PipeData object
        self.flatfile = flatfile
        flatdata = PipeData(config = self.config)
        flatdata.load(self.flatfile)
        ### find flatfields data arrays and store them
        # get sizes of input data
        datalist = self.getarg('datalist')
        if len(datalist) > 0:
            # There are items in datalist -> loop over items
            self.flats = []
            # Check if necessary number of images in flatdata
            if len(flatdata.imgdata) < len(datalist): 
                msg = 'Number of images in flatfield file < '
                msg += 'number of entries in datalist'
                self.log.error('LoadFlat: %s' % msg)
                raise ValueError(msg)
            # Loop through datalist items
            for dataind in range(len(datalist)):
                dataitem = datalist[dataind]
                # Search for dataitem in self.datain images
                if dataitem.upper() in self.datain.imgnames:
                    dataimg = self.datain.imageget(dataitem)
                    self.log.debug('LoadFlat: Found image <%s> to flat'
                                   % dataitem)
                # Search dataitem in self.table columns
                else:
                    try:
                        dataimg = self.datain.table[dataitem]
                        self.log.debug('LoadFlat: Found column <%s> to flat'
                                       % dataitem)
                    except:
                        msg = 'No data found for <%s>' % dataitem
                        self.log.error('LoadFlat: %s' % msg)
                        raise ValueError(msg)
                # Get dimensions - append flat to list
                datasiz = dataimg.shape
                if self.getarg('l0method').upper() != 'NO':
                    datasiz = datasiz[1:]
                flatsiz = flatdata.imgdata[dataind].shape
                self.flats.append(flatdata.imgdata[dataind])
                # Check dimension with flat data
                print(datasiz,flatsiz)
                if len(datasiz) >= len(flatsiz):
                    # Data has >= dimensions than flat -> compare
                    begind = len(datasiz)-len(flatsiz)
                    if datasiz[begind:] != flatsiz:
                        msg = 'Flat "%s" does not fit data - A' % dataitem
                        self.log.error('LoadFlat: %s' % msg)
                        raise ValueError(msg)
                else:
                    # More dimensions in flat data -> report error
                    msg = 'Flat "%s" does not fit data - B' % dataitem
                    self.log.error('LoadFlat: %s' % dataitem)
                    raise ValueError(msg)
        else:
            # Empty datalist -> Flat first image in data with first flat
            datasiz = self.datain.image.shape
            if self.getarg('l0method').upper() != 'NO':
                datasiz = datasiz[1:]
            flatsiz = flatdata.image.shape
            if len(datasiz) >= len(flatsiz):
                # Data has >= dimensions than flat -> compare
                begind = len(datasiz)-len(flatsiz)
                if datasiz[begind:] != flatsiz:
                    self.log.error('LoadFlat: Flat does not fit data - A')
                    raise ValueError('Flat does not fit data - A')
            else:
                # More dimensions in flat data -> report error
                self.log.error('LoadFlat: Flat does not fit data - B')
                raise ValueError('Flat does not fit data - B')
            self.log.debug('LoadFlat: Flatfielding first data image with first flat')
            self.flats=[flatdata.image]
        ### make good pixel map for each detector and add to 
        #flattemp = numpy.abs(data[0,...])+numpy.abs(data[1,...])
        #self.goodpixmap = numpy.ones(data.shape[1:])
        #self.goodpixmap [ numpy.where(flattemp == 0.0)] = 0.0
        # Finish up
        self.flatloaded = 1
        self.log.debug('LoadFlat: done')
    
    def flatfield(self,imgin,flat):
        """ Flatfields an array using flat.
            If r0method != NO then the real image is computed
            This method checks that imagein and flat are compatible
        """
        # Do r0method correction if necessary
        l0method = self.getarg('l0method')
        if l0method != 'NO':
            # Return L0 for chosen method
            if l0method == 'ABS':
                data=numpy.zeros(imgin.shape[1:], dtype=complex)
                data.real=imgin[0,...].copy()
                data.imag=imgin[1,...].copy()
                imgin = abs(data)
            elif l0method =='IM':
                imgin = imgin[1,...].copy()
            else:
                imgin = imgin[0,...].copy()
        # Check flatfield dimension
        datasiz = imgin.shape
        flatsiz = flat.shape
        if len(datasiz) >= len(flatsiz):
            # Data has >= dimensions than flat -> compare
            begind = len(datasiz)-len(flatsiz)
            if datasiz[begind:] != flatsiz:
                msg='Flat does not fit data in file %s'%self.datain.filename
                self.log.error('FlatField: %s' % msg)
                raise ValueError(msg)
        else:
            # More dimensions in flat data -> report error
            msg = 'Flat does not fit data in file %s' % self.datain.filename
            self.log.error('LoadFlat: %s' % msg)
            raise ValueError(msg)
        # Apply flatfield
        imgout = imgin / flat
        return imgout
    
    def reset(self):
        """ Resets the step to the same condition it was when it was
            created. Stored flatfield data is erased and the configuration
            information is cleared.
        """
        self.flatloaded = 0
        self.flatvalue = numpy.zeros([1,1])
        self.flatphase = numpy.zeros([1,1])
        self.flatheader = fits.PrimaryHDU(numpy.array(1))
        self.flatfile = ''
    
    def test(self):
        """ Test Pipe Step Flat Object: Runs basic tests
        """
        # initial log message
        self.log.info('Testing pipe step flat')
        # get testin and a configuration
        if len(self.config) > 2: # i.e. if real config is loaded
            testin = PipeData(config=self.config)
        else:
            testin = PipeData(config=self.testconf)
        # load sample data
        datain=PipeData(config=testin.config)
        #infile = 'mode_chop/120207_000_00HA012.chop.dmd.fits'
        infile = 'mode_chop/120306_000_00HA006.chop.dmd.fits'
        #infile = 'mode_chop/120402_000_00HA035.chop.dmd.fits'
        #infile = 'sharp/sharc2-048485.dmdsqr.fits'
        testfile = os.path.join(datain.config['testing']['testpath'],
                                infile)
        #testfile = '/Users/berthoud/testfit.fits'
        datain.load(testfile)
        if False:
            # change data (make complex number array with
            #              Re=0,1,2,3,4,5,6 . . . in time Im=0)
            dataval=numpy.ones(datain.image.shape)
            dataval[...,1]=0.0
            inclist=numpy.arange(dataval.shape[0])
            incshape=[1+i-i for i in dataval.shape[0:-1] ]
            incshape[0]=dataval.shape[0]
            inclist.shape=incshape
            dataval[...,0]=dataval[...,0]*inclist
        #datain.image=dataval
        # run first flat
        dataout=self(datain)
        #print dataout.image[100,...] # print 100th image
        #print dataout.image[range(0,dataval.shape[0],1000),0,0] # 1 val per img
        dataout.save()
        # final log message
        self.log.info('Testing pipe step flat - Done')

if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
    """
    StepFlat().execute()

""" === History ===
    2013-9-16 Marc Berthoud: Upgraded to new step format
    2010-12-20 Marc Berthoud: Added creation of the GOODPIXMAP image frame
    2010-10-26 Marc Berthoud: Upgraded to new pipedata
    2009-12-01 Marc Berthoud, Ver 0.1.0: Written and tested
"""
