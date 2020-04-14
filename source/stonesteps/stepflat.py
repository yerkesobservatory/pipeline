#!/usr/bin/env python
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

    A specific flat file can be specified in the configuration. Alternatively,
    a flat file is selected by matching specific keywords in the header of
    the input data to files in a specified folder. If no perfect match can be
    found, the best fit is selected. For each subsequent data set the keywords
    in the header are checked against the ones from the initial data file.

    In addition, Flat can add images (HDUs) from the flat file to the data
    for later use (for instance bad pixel maps or distortion information).
    The addfromfile keyword is used for that, it should contain a list of HDU
    names to add.

"""

import os # os library
import numpy # numpy library
import logging # logging object library
from astropy.io import fits
from darepype.drp import DataFits # pipeline data object
from darepype.drp import StepParent # pipe step parent object
from darepype.tools.steploadaux import StepLoadAux # pipestep steploadaux object

class StepFlat(StepLoadAux, StepParent):
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
        self.flatdata = DataFits() # Pipedata object containing the flat file
        # flat file info and header keywords to fit
        self.flatfile = '' # name of selected flat file
        self.fitkeys = [] # FITS keywords that have to fit
        self.keyvalues = [] # values of the keywords (from the first data file)
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
        self.paramlist.append(['reload', 'False',
            'Set to True to look for new flat files for every input'])
        self.paramlist.append(['l0method', 'NO',
            'Method to normalize data: NOne, REal, IMag and ABSolute ' +
            '(default = NO)'])
        self.paramlist.append(['datalist', [],
            'List of data sets to flatten in intput file ' +
            '(default = [] i.e. only flatten image cube in first HDU)'])
        self.paramlist.append(['addfromfile', [],
            'List of data sets from the flat file to add to the output data' +
            '(default = [] i.e. no data to add)'])
        # Get parameters for StepLoadAux, replace auxfile with flatfile
        self.loadauxsetup('flatfile')
        

    def run(self):
        """ Runs the flatfielding algorithm. The flatfielded data is
            returned in self.dataout
        """
        ### Preparation
        # Load flat files if necessary
        if not self.flatloaded or self.getarg('reload'):
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
        ### Add additional image frames from flat file
        for dataitem in self.getarg('addfromfile'):
            ind = self.flatdata.imageindex(dataitem.upper())
            if ind > -1:
                self.dataout.imageset(self.flatdata.imageget(dataitem),
                                      imagename = dataitem,
                                      imageheader = self.flatdata.getheader(dataitem))
                continue
            ind = self.flatdata.tableindex(dataitem.upper())
            if ind > -1:
                self.dataout.tableset(self.flatdata.tableget(dataitem),
                                      tablename = dataitem,
                                      tableheader = self.flatdata.getheader(dataitem))
                continue
            msg = 'No data to add found for <%s>' % dataitem
            self.log.error('Run: %s' % msg)
            raise ValueError(msg)
        # Remove the instrumental configuration HDU
        if 'CONFIGURATION' in self.dataout.imgnames:
            self.dataout.imagedel('CONFIGURATION')    

        ### Finish - cleanup
        # Update DATATYPE
        self.dataout.setheadval('DATATYPE','IMAGE')
        # Add flat file to History
        self.dataout.setheadval('HISTORY','FLAT: %s' % self.flatdata.filename)
        # Update PROCSTAT to level 2
        self.dataout.setheadval('PROCSTAT','LEVEL_2')

    def loadflat(self):
        """ Loads the flat information for the instrument settings
            described in the header of self.datain.

            If an appropriate file can not be found or the file is invalid
            various warnings and errors are returned.
        """
        ### Search for flat and load it into data object
        self.flatdata = self.loadauxfile()
        ### find flatfields data arrays and store them
        # get sizes of input data
        datalist = self.getarg('datalist')
        if len(datalist) == 0:
            # Empty datalist -> Flat first image in data with first flat
            self.checksize(self.datain.image.shape, self.flatdata.image.shape)
            self.log.debug('LoadFlat: Flatfielding first data image with first flat')
            self.flats=[self.flatdata.image]
        else:
            # There are items in datalist -> loop over items
            self.flats = []
            # Check if necessary number of images in flatdata
            if len(self.flatdata.imgdata) < len(datalist):
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
                self.checksize(dataimg.shape, self.flatdata.imgdata[dataind].shape)
                self.flats.append(self.flatdata.imgdata[dataind])
        ### Ensure that data listed in addfromfile is present in flatfile
        for dataitem in self.getarg('addfromfile'):
            if dataitem.upper() in self.flatdata.imgnames:
                pass
            elif dataitem.upper() in self.flatdata.tabnames:
                pass
            else:
                msg = 'No data to add found for <%s>' % dataitem
                self.log.error('LoadFlat: %s' % msg)
                raise ValueError(msg)
        # Finish up
        self.flatloaded = 1
        self.log.debug('LoadFlat: done')

    def flatfield(self,imgin,flat):
        """ Flatfields an array using flat.
            If r0method != NO then the real image is computed
            This method checks that imagein and flat are compatible
        """
        # Check flatfield dimension
        self.checksize(imgin.shape, flat.shape)
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
        # Apply flatfield
        imgout = imgin * flat
        return imgout

    def checksize(self, datashape, flatshape):
        """ Checks that the shape of the flat is comptatible to be used
            with image data of datashape. Raises exceptions otherwise.
        """
        if self.getarg('l0method').upper() != 'NO':
            datashape = datashape[1:]
        if len(datashape) >= len(flatshape):
            # Data has >= dimensions than flat -> compare
            begind = len(datashape)-len(flatshape)
            if datashape[begind:] != flatshape:
                msg='Flat does not fit data in file %s'%self.datain.filename
                self.log.error('FlatField: %s' % msg)
                raise ValueError(msg)
        else:
            # More dimensions in flat data -> report error
            msg = 'Flat does not fit data in file %s' % self.datain.filename
            self.log.error('LoadFlat: %s' % msg)
            raise ValueError(msg)

    def reset(self):
        """ Resets the step to the same condition it was when it was
            created. Stored flatfield data is erased and the configuration
            information is cleared.
        """
        # initialize input and output
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
        if self.config != None and len(self.config) > 2: # i.e. if real config is loaded
            testin = DataFits(config=self.config)
        else:
            testin = DataFits(config=self.testconf)
        # load sample data
        datain=DataFits(config=testin.config)
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
    2016-4-13 Marc Berthoud: Flat is now MULTIPLIED with data not divided
    2013-9-16 Marc Berthoud: Upgraded to new step format
    2010-12-20 Marc Berthoud: Added creation of the GOODPIXMAP image frame
    2010-10-26 Marc Berthoud: Upgraded to new pipedata
    2009-12-01 Marc Berthoud, Ver 0.1.0: Written and tested
"""
