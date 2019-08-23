""" PIPE DATA FITS - Version 1.0.0

    This module defines the pipeline FITS data object. The object
    is used by the pipeline to store raw, intermediate and reduced data.

"""

import numpy # numpy library
import os # os library
import logging # logging library
import time # time library
import configobj # config object library
import validate  # to check input config file is correct
import re # regexp
import gc # garbage collect
from astropy.io import fits # fits libary
from dataparent import DataParent # Pipe Data Parent
from __builtin__ import True

class DataFits(DataParent):
    """ Pipeline Data FITS Object
        The data is stored as a list of (multi-dimensional) images and a
        similar list of tables. Fits  file name and all headers are
        stored as well.
    """

    # File Name Fit: Regexp expression that fits valid filenames
    filenamefit = r'\.(fits|fts|fit)\Z' # will fit .fits or .fts at end of name

    def __init__(self, filename='', config=None):
        """ Constructor: Initialize data object and variables
        """
        # set up internal variables
        #     Remark: to preserve order of images no dictionary is used
        self.filename=filename # Name of loaded file
        # Image Variables:
        #   At all times len(imgdata) == len(imgnames) == len(imgheads)
        #   - If only loadhead() was called or first HDU has no data:
        #     imghead[0]=header, imgdata[0]=None, imgnames[0]='Primary Header'
        self.imgdata=[] # List with image data - data is in array of int/float
        self.imgnames=[] # List with image names (all uppercase)
        self.imgheads=[] # List with image headers
        # Table Variables:
        #   At all times len(tabdata) == len(tabnames) == len(tabheads)
        self.tabdata=[] # List with tables each item is a record array
        self.tabnames=[] # List with table names (all uppercase)
        self.tabheads=[] # List with table headers
        # set up logging and configuration
        self.log=logging.getLogger('pipe.datafits') # Logger
        self.config = None # Configuration
        self.setconfig(config) # Set the Configuration
        # Count to decrease number of frequent log messages
        self.logcount=0
        # if file exists, load it
        if os.path.exists(filename):
            self.load(filename)
        elif filename != '': # user specified a non-existant filename
            msg = 'No such file %s' %filename
            self.log.error(msg)
            raise IOError(msg)
        self.log.debug('Init: done')

    def __getattr__(self, name):
        """ Get Attribute function that allows to access these attributes:
            - header: header of primary HDU
            - image: the first image
            - table: the first table
            - filenamebase: filename start - depreciated
            From DataParent.__getattr__:
            - filenamebegin: filename start
            - filenameend: filename end
        """
        # return data
        if name == 'data':
            # return image if available
            if len(self.imgdata) > 0 and self.imgdata[0] != None:
                return self.imgdata[0]
            # else return a table
            elif len(self.tabdata) > 0:
                return self.tabdata[0]
            else:
                msg="'%s' object has no data" % (type(self).__name__)
                self.log.warn(msg)
                return None
        # return image if it's requested - raise error if no images present
        if name == 'image':
            if len(self.imgdata) == 0:
                msg="'%s' object has no image data" % (type(self).__name__)
                self.log.warn(msg)
                return None
            elif self.imgdata[0] is None:
                msg="'%s' object has no image data" % (type(self).__name__)
                self.log.warn(msg)
            return self.imgdata[0]
        # return table if it's requested - raise error if no tables present
        if name == 'table':
            if len(self.tabdata) == 0:
                msg="'%s' object has no table data" % (type(self).__name__)
                self.log.warn(msg)
                return None
            return self.tabdata[0]
        # return header if it's requested - raise error if no images present
        if name == 'header':
            if len(self.imgheads) == 0:
                msg="'%s' object has no header data" % (type(self).__name__)
                self.log.warn(msg)
                return None
            return self.imgheads[0]
        # return filenamebase if it's requested
        if name == 'filenamebase':
            self.log.warn('Filenamebase is depreciated - use filenamebegin')
            return self.filenamebegin
        # run parent function (filenamebegin and filenameend)
        return super(DataFits, self).__getattr__(name)
        # raise error if attribute is unknown
        msg="'%s' object has no attribute '%s'" % (type(self).__name__,name)
        raise AttributeError(msg)

    def __setattr__(self, name, value):
        """ Set Attribute function that allows to access the first
            image as PipeData.image
        """
        # set the data to image
        if name == 'data':
            # if it's a table put it in table
            if issubclass(value.__class__, numpy.recarray):
                self.table = value
            # else assume it's an image
            else:
                self.image = value
        # set the image if it's requested
        if name == 'image':
            if len(self.imgdata) > 0:
                self.imgdata[0] = value
            else:
                self.imgdata=[value]
                self.imgnames=['PRIMARY IMAGE']
                self.imgheads=[fits.Header()]
        # set the table if it's requested
        elif name == 'table':
            if len(self.tabdata) > 0:
                self.tabdata[0] = value
            else:
                self.tabdata=[value]
                self.tabnames=['PRIMARY TABLE']
                self.tabheads=[fits.Header()]
        # set the header if it's requested
        elif name == 'header':
            if len(self.imgdata) > 0:
                self.imgheads[0] = value
            else:
                self.imgdata=[None]
                self.imgnames=['PRIMARY HEADER']
                self.imgheads=[value]
        # else pass the command to the parent function
        else:
            object.__setattr__(self, name, value)

    def loadhead(self,filename='', dataname=''):
        """ Loads and returns the primary header of the FITS file given.
            This also checks for file existence and type. The pipedata object
            is filled with only that header.
            - filename: a string with the name of the file to load
              (if omitted, self.filename is used)
            - dataname specified the EXTNAME value of the header to be loaded
              If such a header is not found, or dataname=='', the first header
              is loaded. This option should be used if the main file
              information is not in the primary header.
        """
        ### check for file existence, type and get primary header
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # read fits header, checks for existing valid fits file
        try:
            hdus = fits.open(filename)
        except IOError, error:
            # IOError: 2 No such file or directory (missing file)
            self.log.error('LoadHead: no such file or head read error, '
                           + filename)
            raise error
        except IndexError, error:
            # IndexError: list index out of range (invalid format for fits)
            self.log.error('LoadHead: '+filename+' is not a valid FITS file')
            raise error
        except TypeError, error:
            # TypeError: invalid filename (not a string)
            self.log.error('LoadHead: filename is invalid type')
            raise error
        # No dataname -> return primary header
        if dataname == '':
            header = hdus[0].header
        # Look for correct dataname
        else:
            hdui = 0
            found = 0
            while not found:
                if 'EXTNAME' in hdus[hdui].header.keys():
                    if hdus[hdui].header['EXTNAME'].strip() == dataname:
                        found = 1
                        header = hdus[hdui].header
                hdui += 1
            if not found:
                msg = "loadhead: HDU with EXTNAME=%s not found" % dataname
                self.log.error(msg)
                raise ValueError(msg)
        # Fill in the header, if necessary fill in data
        self.imgheads = [header]
        self.imgdata = [None]
        try:
            self.imgnames = [header['EXTNAME'].upper()]
        except:
            self.imgnames = ['PRIMARY HEADER']
        # Fill the filename
        self.filename = filename
        self.log.debug('LoadHead: done')
        hdus.close()

    def load(self,filename=''):
        """ Load a file into the data object. The file can be a raw HAWC
            file or a file saved by pipedata.
            - filename: a string with the name of the file to load
              (if omitted, self.filename is used)
            - raises various errors if file can not be read
        """
        ### Clear file data
        self.imgdata, self.imgheads, self.imgnames = [], [], []
        self.tabdata, self.tabheads, self.tabnames = [], [], []
        ### Get filename and file checks
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # check for file existence, type and get primary header
        #   fills first entry into imgdata, imgnames
        self.loadhead(filename)
        ### Read File
        hdus=fits.open(filename,memmap=False)
        ### Collect images / Load them
        # Search for ImageHDUs (does not include PrimaryHDU)
        imgind=[i for i in range(len(hdus))
                if isinstance(hdus[i], fits.ImageHDU) ]
        imgn = len(imgind) # store number of images
        # get naxis and naxis1 from primary header, check if keywords exist
        try:
            naxis=self.getheadval('NAXIS')
            if naxis > 0:
                naxis1=self.getheadval('NAXIS1',errmsg=False)
            else:
                naxis1=0
        except KeyError:
            # KeyError: keyword not found (no keywords available)
            self.log.warn('Load: missing naxis keywords in fits file '
                           + filename)
        # Check if file has no image (i.e. if naxis==1 and naxis1==0)
        # -> No change, primary image stays as it was from loadhead()
        if naxis * naxis1 == 0:
            self.log.info('Load: No image data in first HDU')
        # else: Load first HDU data if there is image data in it
        else:
            imgn += 1
            # Load data
            self.imgdata[0] = hdus[0].data
            # Set image name
            if 'EXTNAME' in self.imgheads[0]:
                self.imgnames[0] = self.imgheads[0]['EXTNAME'].upper()
            else:
                self.imgnames[0] = 'PRIMARY IMAGE'
        # Message with number of images
        if imgn == 0:
            self.log.debug('Load: No image in file %s' % filename)
        elif imgn == 1:
            self.log.debug('Load: 1 image in file %s' % filename)
        else:
            self.log.debug('Load: %d images in file %s' %(imgn, filename))
        # Load subsequent images
        for ind in imgind:
            # get name
            if 'EXTNAME' in hdus[ind].header.keys():
                self.imgnames.append(hdus[ind].header['EXTNAME'].upper())
            else:
                self.imgnames.append('SECONDARY IMAGE %d' % ind)
            # get data and header
            self.imgdata.append(hdus[ind].data)
            self.imgheads.append(hdus[ind].header)
        ### Look for Tables / Read them
        # Search for BinTableHDUs
        tabind=[i for i in range(len(hdus))
                if isinstance(hdus[i], fits.BinTableHDU) ]
        tabn = len(tabind)
        # Messages on number of tables
        if tabn == 0:
            self.log.debug('Load: No table in file %s' % filename)
        elif tabn == 1:
            self.log.debug('Load: %d table in file %s' %(tabn, filename))
        else:
            self.log.debug('Load: %d tables in file %s' % (tabn, filename))
        # Load all tables
        for ind in tabind:
            # if table is empty -> warning and skip it
            try:
                if hdus[ind].data is None:
                    msg = 'Load: Table in HDU number %d has no data' % ind
                    msg += ' -> Ignoring this HDU'
                    self.log.warn(msg)
                    continue
            except:
                msg = 'Load: Problem loading table in HDU number %d' % ind
                msg += ' -> Ignoring this HDU'
                self.log.warn(msg)
                continue

            # get name
            if 'EXTNAME' in hdus[ind].header.keys():
                self.tabnames.append(hdus[ind].header['EXTNAME'].upper())
            else:
                if ind > 0:
                    self.tabnames.append('SECONDARY TABLE')
                else:
                    self.tabnames.append('PRIMARY TABLE')
            # get data, header
            self.tabdata.append(numpy.rec.array(hdus[ind].data))
            self.tabheads.append(hdus[ind].header)
        # Close the file
        hdus.close()
        gc.collect()
        ### Message
        self.log.debug('Load: loaded fits file')

    def save(self,filename = None):
        """ Save the data in the object to the specified file. Existing files are
            overwritten.
            - filename: a string with the file name to save to
        """
        # get file name
        if filename == None:
            filename = self.filename
        # update pipeline keywords
        self.setheadval('PIPEVERS', DataParent.pipever.replace('.','_'),
                        'Pipeline Version')
        self.setheadval('FILENAME',os.path.split(filename)[-1])
        self.setheadval('DATE',time.strftime('%Y-%m-%dT%H:%M:%S'))
        ### save file
        # make the data primary HDU -> List
        hdul = [] # hdu list
        for i in range(len(self.imgnames)):
            if i == 0:
                hdui = fits.PrimaryHDU(self.imgdata[i],self.imgheads[i])
            else:
                hdui = fits.ImageHDU(self.imgdata[i],self.imgheads[i])
            hdui.header['EXTNAME'] = (self.imgnames[i].upper(),'ID of the HDU')
            hdul.append(hdui)
        if len(hdul) == 0: # no IMAGE HDUs added
            hdul.append(fits.PrimaryHDU(None,self.imgheads[0]))
        # make hdus for tables
        for i in range(len(self.tabnames)):
            hdut = fits.BinTableHDU(self.tabdata[i],self.tabheads[i])
            hdut.header['EXTNAME'] = (self.tabnames[i].upper(),'ID of the HDU')
            hdul.append(hdut)

        # make an HDU list
        hdulist=fits.HDUList(hdul)
        # save the file (produce errors if not successful)
        try:
            hdulist.writeto(filename, output_verify='fix', clobber=True)
        except TypeError, error:
            self.log.error('Save: filename is invalid type')
            raise error
        except IOError, error:
            self.log.error('Save: Failed to write fits file to '+filename)
            raise error
        self.log.debug('Save: wrote FITS file %s' % filename)

    def copy(self):
        """ Returns a copy of self
        """
        out=DataFits(config=self.config) # create new object
        # copy all images
        out.imgnames=self.imgnames[:]
        out.imgdata=[]
        out.imgheads=[]
        for imgi in range(len(self.imgdata)):
            if self.imgdata[imgi] is not None:
                out.imgdata.append(self.imgdata[imgi].copy())
            else:
                out.imgdata.append(None)
            out.imgheads.append(self.imgheads[imgi].copy())
        # copy tables
        out.tabnames=self.tabnames[:]
        out.tabdata=[]
        out.tabheads=[]
        for tabi in range(len(self.tabdata)):
            out.tabdata.append(self.tabdata[tabi].copy())
            out.tabheads.append(self.tabheads[tabi].copy())
        # copy filename
        out.filename=self.filename
        # return message and new object
        self.log.debug('Copy: done')
        return out

    def mergedata(self,other):
        """ Appends images and table from the new object to the existing
            data. The keywords from the primary headers are merged, in doubt
            the keyword from the local object is kept. The extension headers
            are left unchanged. If self has no data, the data from other is
            copied into the object.
            - other: a pipedata object with the same kind of data
            - raises ValueErrors if other has invalid format
        """
        ### check if other is correct data type
        if not isinstance(other,DataFits):
            self.log.error('Mergedata: other data is not PipeData type')
            raise TypeError('other data is not PipeData type')
        ### If self has no images and no tables: Copy All data from other
        if self.imgnames==[] and self.tabnames==[]:
            self=other.copy()
            return
        ### Check data compatibility
        # check if image number and names agree
        if self.imgnames != other.imgnames and self.imgnames != []:
            msg = 'other data has different image list'
            self.log.error('Mergedata: ' + msg)
            raise ValueError(msg)
        # check if table number and names agree
        if self.tabnames != other.tabnames:
            msg = 'other data has different table list'
            self.log.error('Mergedata: ' + msg)
            raise ValueError(msg)
        ### append the images (make new, copy, move)
        for imgi in range(len(self.imgdata)):
            # get new shape
            ssel=self.imgdata[imgi].shape
            soth=other.imgdata[imgi].shape
            # check dimensions, set new shape
            msg = 'dimension mismatch at image ' + self.imgnames[imgi]
            # if number of dimensions are equal
            if len(ssel) == len(soth):
                # check shape
                if ssel[1:] != soth[1:]:
                    self.log.error('MergeData: ' + msg)
                    raise ValueError(msg)
                nsel=ssel[0]
                noth=soth[0]
            # if number of other dimensions = number of self - 1
            elif len(ssel) == len(soth)+1:
                # check shape
                if ssel[1:] != soth:
                    self.log.error('MergeData: ' + msg)
                    raise ValueError(msg)
                nsel=ssel[0]
                noth=1
            # if dimensions are unmatched
            else:
                self.log.error('MergeData: ' + msg)
                raise ValueError(msg)
            nnew = nsel + noth
            snew = list(ssel)
            snew[0] = nnew
            # allocate new image
            newimg = numpy.empty(snew, dtype=self.imgdata[imgi].dtype)
            # copy images to new
            newimg[0:nsel, ...] = self.imgdata[imgi]
            newimg[nsel:nnew, ...] = other.imgdata[imgi]
            self.imgdata[imgi] = newimg
        ### append tables
        for tabi in range(len(self.tabdata)):
            # Check if table formats agree
            if self.tabdata[tabi].dtype != other.tabdata[tabi].dtype :
                msg='other table columns are incompatible'
                self.log.error('MergeData: '+msg)
                raise TypeError(msg)
            # make new shape
            nsel=len(self.tabdata[tabi])
            noth=len(other.tabdata[tabi])
            nnew=nsel+noth
            newshape=(nnew,)
            newtable=numpy.empty(newshape,dtype=self.tabdata[tabi].dtype)
            # copy existing table and put to self.data
            newtable[0:nsel]=self.tabdata[tabi]
            newtable[nsel:nnew]=other.tabdata[tabi]
            self.tabdata[tabi]=newtable
        ### combine headers
        self.mergehead(other)
        self.log.debug('MargeData: done')

    def mergehead(self, other):
        """ Merges the header of another data object to the existing header.
            Most of the header of the new data is ignored, COMMENT,
            HISTORY keywords are copied.
        """
        # get self and other cards and parents
        othercards = other.header.cards
        selfcards = self.header.cards
        # add history keywords (no duplicates)
        otherhist=[card.value for card in othercards
                  if card.keyword=='HISTORY'] # hist values for other
        selfhist=[card.value for card in selfcards
                  if card.keyword=='HISTORY'] # hist values for self
        for hist in otherhist: # loop through other history values
            if not hist in selfhist: # only add if not duplicate
                self.header.add_history(hist)
        # add comment keywords (no duplicates)
        othercomm=[card.value for card in othercards
                  if card.keyword=='COMMENT'] # comm values for other
        selfcomm=[card.value for card in selfcards
                  if card.keyword=='COMMENT'] # comm values for self
        for comm in othercomm: # loop through other comment values
            if not comm in selfcomm: # only add if not duplicate
                self.header.add_comment(comm)
        # Go through keywords listed in headmerge: assume self is first
        headmerge = self.config['headmerge']
        for key in headmerge.keys():
            if key in self.header and key in other.header:
                selfval = self.header[key]
                otherval = other.header[key]
                operation = headmerge[key].upper()
                if operation == 'LAST': selfval = otherval
                elif operation == 'MIN': selfval = min(selfval,otherval)
                elif operation == 'MAX': selfval = max(selfval,otherval)
                elif operation == 'SUM': selfval += otherval
                elif operation == 'OR': selfval = selfval | otherval
                elif operation == 'AND': selfval = selfval & otherval
                elif operation == 'CONCATENATE':
                    if ',' in selfval:
                        vlist = selfval.split(',')
                    else:
                        vlist = [selfval]
                    if otherval not in vlist:
                        vlist.append(otherval)
                        selfval = ','.join(sorted(vlist))
                elif operation == 'DEFAULT':
                    if type(selfval) is str:
                        selfval = 'UNKNOWN'
                    elif type(selfval) is int or type(selfval) is long:
                        selfval = -9999
                    elif type(selfval) is float:
                        selfval = -9999.0
                self.header[key] = selfval
        self.log.debug('MergeHead: done')

    def copyhead(self,other,name=None,overwrite=True):
        """ Similar to mergehead, except copy all header keywords, comments,
            and history from other.  Will overwrite existing cards, unless
            overwrite flag is set to False.  Exceptions: HISTORY and COMMENT
            cards are always appended to the end of the list of such keywords
            present in self.

            If name is None, use the header from the first HDU (self.header)
        """

        if name is None:
            head1 = other.header
            head2 = self.header
        else:
            head1 = other.getheader(name)
            head2 = self.getheader(name)
        try:
            nhist = len(head1['history'])
        except KeyError:
            nhist = 0
        try:
            ncomm = len(head1['comment'])
        except KeyError:
            ncomm = 0

        for k in head1.iterkeys():
            if k in ('COMMENT','HISTORY'): # handle later
                pass
            elif k in head2.keys():
                if overwrite:
                    head2[k] = (head1[k],head1.comments[k])
            else:
                head2[k] = (head1[k],head1.comments[k])
        for i in xrange(nhist):
            head2.add_history(head1['history'][i])
        for i in xrange(ncomm):
            head2.add_comment(head1['comment'][i])
        self.log.debug("Copydata: done")

    def adddata(self,images,tablerow):
        """ Appends images and a row to primary table
            - images: a dictionary with images to add
              (images are indexed by correct image names)
            - tablerow: one row of table values for the images (can be list)
                        This row is added to the first table
            - raises Errors if images or tablerow have wrong
              format or dimensions
        """
        # Comment: This function makes too many checks that would just
        #          raise exceptions when the error happens
        ### check input
        # check type of images
        if not isinstance(images, dict):
            self.log.error('Adddata: images needs to be a dictionary')
            raise TypeError('images needs to be a dictionary')
        # check keywords and data dimensions
        for key in images.keys():
            # check keyword / get index
            try:
                ind = self.imgnames.index(key.upper())
            except ValueError, error:
                msg = 'invalid image name (%s) in images' % key
                self.log.error('Adddata: ' + msg)
                raise error
            # check data type / get shape
            try:
                news = images[key].shape
            except AttributeError, error:
                msg = 'invalid image type for image %s' % key
                self.log.error('Adddata: ' + msg)
                raise error
            olds = self.imgdata[ind].shape
            # compare shapes (new can be part of old)
            if len(news) == len(olds):
                if news[1:] != olds[1:]:
                    msg = 'incompatible image shape in %s' % key
                    self.log.error('Adddata: ' + msg)
                    raise ValueError(msg)
            else:
                if news != olds[1:]:
                    msg = 'incompatible image shape in %s' % key
                    self.log.error('Adddata: ' + msg)
                    raise ValueError(msg)
        # check if table format is correct (by converting it)
        if len(tablerow) != len(self.table[0]) :
            self.log.error('AddData: table row has invalid length')
            raise ValueError, 'table row has invalid length'
        try:
            tablerow=numpy.rec.array(tablerow,dtype=self.table.dtype)
        except (ValueError,TypeError), error:
            self.log.error('AddData: table row has invalid element type')
            raise ValueError, 'table row has invalid element type'
        ### append data
        for key in images.keys():
            # get index
            ind = self.imgnames.index(key.upper())
            # get new size
            news = images[key].shape
            olds = self.imgdata[ind].shape
            oldn = olds[0]
            if len(news) == len(olds):
                newn = news[0]
            else:
                newn = 1
            totn = oldn + newn
            tots = list(olds)
            tots[0]=totn
            # allocate new data
            newimg = numpy. empty(tots, dtype=self.imgdata[ind].dtype)
            # copy images and set new data
            newimg[0:oldn, ...] = self.imgdata[ind]
            newimg[oldn:totn, ...] = images[key]
            self.imgdata[ind] = newimg
        ### append table
        # make new shape
        totn = len(self.table) + 1
        newshape=(totn,)
        newtable=numpy.empty(newshape,dtype=self.table.dtype)
        # copy existing table and put to self.data
        newtable[0:totn-1]=self.table
        newtable[totn-1]=tablerow
        self.table=newtable
        self.log.debug('AddData: done')

    def copydata(self, other, dataname):
        """ Copies data (image or table) form another pipedata object to the
            current object. If an object of that name already exists it is
            overwritten. Both the data and the header are copied.
        """
        # Check type and presence of data in other object
        isimg = 0
        if dataname.upper() in other.imgnames:
            isimg = 1
            othind = other.imageindex(dataname)
        elif dataname.upper() in other.tabnames:
            othind = other.tableindex(dataname)
        else:
            self.log.error('CopyData: dataname not found')
            raise ValueError('dataname not found')
        # Check if name exists in local data, delete if it's different type
        if dataname.upper() in self.imgnames:
            if isimg < 1:
                self.imagedel(dataname)
        elif dataname.upper() in self.tabnames:
            if isimg > 0:
                self.tabledel(dataname)
        # Copy / Overwrite object
        if isimg>0:
            self.imageset(other.imgdata[othind].copy(), dataname,
                          other.imgheads[othind].copy())
        else:
            self.tableset(other.tabdata[othind].copy(), dataname,
                          other.tabheads[othind].copy())
        self.log.debug('CopyData: done')

    def imageindex(self,imagename=None):
        """ Returns the index of an image, such that the
            images can be accessed pipedata.imgdata[index]
            - imagename: The name of the requested image (None for first image)
        """
        # empty name -> return first image
        if imagename == None:
            self.log.debug('Imageindex: return 0 - done')
            return 0
        # check for valid name
        try:
            ind = self.imgnames.index(imagename.upper())
        except ValueError, error:
            msg = 'invalid image name (%s)' % imagename
            self.log.error('Imageindex: '+msg)
            raise error
        # return index
        self.log.debug('Imageindex: done')
        return ind

    def imageget(self,imagename=None):
        """ Returns an image
            - imgname: The name of the requested image (None for first image)
        """
        # get index
        ind = self.imageindex(imagename)
        # return image
        self.log.debug('Imageget: done')
        return self.imgdata[ind]

    def imageset(self, imagedata, imagename=None, imageheader=None, index=-1):
        """ Sets an image. This should be used to add a new image or if the
            image data should be replaced. The index flag allows to determine
            the position of the image in the image list.
            - imagedata: A multi dimensional array containing the image data
            - imagename: The name of the image to set (None for first image)
            - imageheader: The header for the image to set
            - index: Indicates the position of the image in the image list
                     Ignored if == -1.
        """
        totalindex = len(self.imgdata)
        # If index is not set, get a valid index
        if index < 0:
            if imagename == None:
                index = 0
                if totalindex > 0:
                    imagename = self.imgnames[0]
                else:
                    imagename = 'PRIMARY'
            elif imagename.upper() in self.imgnames:
                index = self.imgnames.index(imagename.upper())
            else: # name given, but not an existing name
                if totalindex > 0 and self.imgdata[0] is None:
                    index = 0 # add to primary HDU, if primary HDU empty
                else:
                    index = totalindex # add to end otherwise
        # re-use header, if none specified
        if imageheader is None:
            if index < totalindex:
                imageheader = self.imgheads[index]
        # Set image
        if index == 0:
            hdu = fits.PrimaryHDU(imagedata,header=imageheader)
        else:
            hdu = fits.ImageHDU(imagedata,header=imageheader)
        if index < totalindex: # overwriting an existing image
            self.imgnames[index] = imagename.upper()
            self.imgdata[index]  = imagedata
            self.imgheads[index] = hdu.header
        else:
            self.imgnames.append(imagename.upper())
            self.imgdata.append(imagedata)
            self.imgheads.append(hdu.header)
        self.log.debug('Imageset: done')

    def imagedel(self, imagename=None):
        """ Removes the image specified by imagename.
            - imagename: The name of the image to delete (None for first image)
        """
        # get index
        ind = self.imageindex(imagename)
        # delete image
        del self.imgnames[ind]
        del self.imgdata[ind]
        del self.imgheads[ind]
        self.log.debug('Imagedel: done')

    def tableindex(self, tablename=None):
        """ Returns the index of a table, such that the tables can be
            accessed pipedata.tabdata[index]
            - tablename: The name of the requested table (None for first table)
        """
        # Check if tables are present
        if len(self.tabnames) == 0:
            msg = 'no tables in data'
            self.log.error('Imageindex: ' + msg)
            raise RuntimeError(msg)
        # empty name -> return first table
        if tablename == None:
            self.log.debug('Tableindex: return 0 - done')
            return 0
        # check for valid name
        try:
            ind = self.tabnames.index(tablename.upper())
        except ValueError, error:
            msg = 'invalid table name (%s)' % tablename
            self.log.error('Tableindex: '+msg)
            raise error
        # return index
        self.log.debug('Tableindex: done')
        return ind

    def tableget(self, tablename = None):
        """ Returns a table
            - tablename: The name of the requested table (None for first table)
        """
        # get index
        ind = self.tableindex(tablename)
        # return table
        self.log.debug('Tableget: done')
        return self.tabdata[ind]

    def tableset(self, tabledata, tablename = None, tableheader = None,
                 index = -1):
        """ Sets a table. This should be used to add a new image or if the
            table data should be replaced. The index flag allows to determine
            where (in the order of tables) the table should be.
            - tabledata: A record array with table contents
            - tablename: The name of the table to set (None for first table)
            - index: Indicates the position of the table in the table list
                     Ignored if == -1.
        """
        # If index is not set, get valid index
        if index < 0:
            if tablename == None:
                index = 0
            # if image exists - replace
            elif tablename.upper() in self.tabnames:
                index = self.tabnames.index(tablename.upper())
        # Set table
        if index > len(self.tabnames)-1 or index < 0:
            self.tabnames.append(tablename.upper())
            self.tabdata.append(tabledata)
            # make sure tableheader is valid
            if tableheader == None:
                tableheader = fits.Header()
            self.tabheads.append(tableheader)
        else:
            self.tabnames[index] = tablename.upper()
            self.tabdata[index] = tabledata
            if tableheader != None:
                self.tabheads[index] = tableheader
        self.log.debug('Tableset: done')

    def tabledel(self, tablename = None):
        """ Removes the table specified by tablename.
            - tablename: The name of the table to be deleted        # Remove Imag columns
        #for tabcol in ['R array Imag','T array Imag','Chop Offset Imag']:
        #    self.dataout.table.columns.del_col(tabcol)

        # Remove Imag columns
        #for tabcol in ['R array Imag','T array Imag','Chop Offset Imag']:
        #    self.dataout.table.columns.del_col(tabcol)

        # Remove Imag columns
        #for tabcol in ['R array Imag','T array Imag','Chop Offset Imag']:
        #    self.dataout.table.columns.del_col(tabcol)

        # Remove Imag columns
        #for tabcol in ['R array Imag','T array Imag','Chop Offset Imag']:
        #    self.dataout.table.columns.del_col(tabcol)

                         (None for first table)
        """
        # get index
        ind = self.tableindex(tablename)
        # delete table
        del self.tabnames[ind]
        del self.tabdata[ind]
        del self.tabheads[ind]
        self.log.debug('Tabledel: done')

    def tableaddcol(self, colname, array, tablename = None, dtype=None):
        """ Adds a column to the table. If the table under tablename
            doesn't exist, it is created.
            - colname: a string specifying the name of the new column
            - array: a list containing the values for the new column
            - tablename: The name of the affected table (None for first table)
            - dtype: (optional) data type for the new column
            WARNING: This may not work if the elements in the table have more
                     than one dimension.
        """
        ### Make new data type
        # make a basic array
        array=numpy.asarray(array)
        # set correct data type
        if dtype is None:
            dtype=array.dtype
        # get additional dimension
        if len(array.shape) > 1:
            newtype = ( colname, dtype, array.shape[1] )
        else:
            newtype = ( colname, dtype )
        ### Get table index (-1 if table needs to be created)
        tabind = -1
        if tablename is None:
            if len(self.tabnames) > 0:
                tabind = 0
            else:
                tablename = "Table"
        else:
            if tablename.upper() in self.tabnames:
                tabind = self.tableindex(tablename)
        ### Make new table if necessary
        if tabind < 0:
        # If there is no existing table:
            # Get new table data type
            newdtype = numpy.dtype([newtype])
            # Make new table
            newtable = numpy.empty(array.shape[0],dtype=newdtype)
            # Add to list of tables
            tabind = len(self.tabnames)
            self.tabdata.append(newtable)
            self.tabnames.append(tablename.upper())
            self.tabheads.append(fits.Header())
        else:
        ### Else (table exists) base new table on old table
            table = self.tabdata[tabind]
            # check if dimension of new value array is correct
            if len(table) != len(array):
                msg = 'column array len (%d) != table len (%d)' % (len(array), len(table))
                self.log.error('TableAddCol: ' + msg)
                raise ValueError(msg)
            # get new data type for table
            newdtype = numpy.dtype(table.dtype.descr+[newtype])
            # make new table
            newtable = numpy.empty(table.shape,dtype=newdtype)
            # fill old table values
            for field in table.dtype.fields:
                newtable[field] = table[field]
        # fill new table values
        newtable[colname] = array
        # copy new table to self.table
        self.tabdata[tabind] = newtable
        self.log.debug('TableAddCol: done')

    def tableaddrow(self, tablerow, tablename = None):
        """ Adds a row to the data table
            - tablerow: A list containing the elements of the row to be added
            - tablename: The name of the affected table (None for first table)
        """
        # Get table
        tabind = self.tableindex(tablename)
        table = self.tabdata[tabind]
        # check if tablerow format is correct (by converting it)
        if len(tablerow) != len(table[0]) :
            self.log.error('TableAddRow: table row has invalid length')
            raise ValueError, 'table row has invalid length'
        try:
            tablerow=numpy.rec.array(tablerow,dtype=table.dtype)
        except (ValueError,TypeError):
            self.log.error('TableAddRow: table row has invalid element type')
            raise ValueError, 'table row has invalid element type'
        # Add to table
        self.tabdata[tabind] = numpy.insert(table, len(table), tablerow)
        self.log.debug('TableAddRow: done')

    def tabledelcol(self, colname, tablename = None):
        """ Deletes a column of the data table.
            - colname: A string or list of strings containing the
                       name(s) of the column(s) to delete
            - tablename: The name of the affected table
                         (None for first table)
        """
        # Get table
        tabind = self.tableindex(tablename)
        table = self.tabdata[tabind]
        # Check if colname is valid, make sure it's a list
        olddt = table.dtype
        if isinstance(colname,str):
            colname = [colname]
        if isinstance(colname,(list,tuple)):
            for c in colname:
                if c not in olddt.names:
                    msg = 'Invalid column name %s' % c
                    self.log.error('TableDelCol: ' + msg)
                    raise ValueError(msg)
        else:
            msg = "Invalid colname '%s'.  Must be string or list/tuple" %colname
            self.log.error('TableDelCol: ' + msg)
            raise ValueError(msg)
        # if all columns deleted, clear table
        if len(olddt) - len(colname) <= 0:
            del self.tabdata[tabind]
            del self.tabnames[tabind]
            del self.tabheads[tabind]
        # else -> remove column(s) from table
        else:
            # if it's a pyfits.FITS_rec -> Treat specially
            # (this was pyfits.fitsrec.FITS_rec)
            # test for FITS_rec.  ndarray and recarray don't support names method
            try:
                names   = table.names
                names   = table.names # list of fields
                formats = table.columns.formats
                dims    = table.columns.dims
                units   = table.columns.units

                cols = []
                for n,f,d,u in zip(names,formats,dims,units):
                    if n not in colname:
                        cols.append(fits.Column(name=n,format=f,dim=d,unit=u,
                            array=table.field(n)))
                tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
                self.tabdata[tabind] = tbhdu.data
                self.tabheads[tabind] = tbhdu.header
            except AttributeError: # assume table is a regular record array
                newnames = [n for n in olddt.names if n not in colname]
                newtable = table[newnames]
                self.tabdata[tabind] = newtable

        self.log.debug('TableDelCol: done')

    def tabledelrow(self, index, tablename = None):
        """ Deletes a row of the data table.
            - index: The index of the row to delete
            - tablename: The name of the affected table (None for first table)
        """
        # Get table
        tabind = self.tableindex(tablename)
        table = self.tabdata[tabind]
        # Check if index is valid
        if index >= len(table):
            self.log.error('TableDelRow: Invalid row index %d' % index)
            raise ValueError('Invalid row index %d' % index)
        # if table has one row -> clear table
        if len(table) < 2:
            del self.tabdata[tabind]
            del self.tabnames[tabind]
        # else -> remove row from table
        else:
            # make index list
            #indlist = [i for i in range(len(table)) if i != index]
            # fill in the values
            #newtable = table[indlist]
            # copy new table
            #self.tabdata[tabind] = newtable.copy()
            self.tabdata[tabind] = numpy.delete(table,index).copy()
        self.log.debug('TableDelRow: done')

    def tablemergerows(self, rows):
        """ Merges several table rows into a single row that is returned.
            Each column is merged according to the rules defined in the
            [table] section of the configuration file.

            Usage: - Get lines to be merged from input table
                     Ex: inrows = indata.table[10:20]
                   - Call tablemergerows
                     Ex: outrow = indata.tablemergerows(inrows)
                   - Add / Insert row to output table
                     Ex: outdata.tableaddrow(outrow)
        """
        # check if rows has same format than table
        try:
            _newdtype = rows.dtype
        except AttributeError, error:
            msg = 'input rows are incorrect data type'
            self.log.error('TableMergeRows: %s' % msg)
            raise error
        #if newdtype != table.dtype :
        #    msg = 'input rows have different format than table'
        #    self.log.error('TableMergeRows: %s' % msg)
        #    raise TypeError(msg)
        # make output table row (copy from first row)
        outrow = (rows[0:1].copy())[0]
        # run through columns and merge values
        for colname in rows.dtype.names:
            # get merge function
            try:
                funct = self.config['table'][colname.lower()].lower()
            except :
                # if keyword is not available
                self.log.warn('TableMergeRows: Missing table merge entry for '
                    + 'column -%s- returning first row value' % colname)
                funct = 'first'
            # Try to run the function
            # Comment: float() is necessary, otherwise values get messed up
            try:
                if funct == 'first':
                    pass
                elif funct == 'last':
                    outrow[colname] = rows[-1][colname]
                elif funct == 'min':
                    outrow[colname] = float(numpy.nanmin(rows[colname]))
                elif funct == 'max':
                    outrow[colname] = float(numpy.nanmax(rows[colname]))
                elif funct == 'med':
                    outrow[colname] = float(numpy.nanmedian(rows[colname]))
                elif funct == 'avg':
                    outrow[colname] = float(numpy.nanmean(rows[colname] ))
                elif funct == 'sum':
                    outrow[colname] = float(numpy.nansum(rows[colname]))
                elif funct == 'wtavg':
                    tmp = float(numpy.nansum(rows[colname]*rows['Samples']))
                    tmpsum = float(numpy.nansum(rows['Samples']))
                    if tmpsum > 0:
                        outrow[colname] = tmp/tmpsum
                    else:
                        outrow[colname] = tmp*0.0
                else:
                    self.log.warn('TableMergeRows: Unknown operation -'+funct
                     +'- for column -'+colname+'- returning first row value')
            except NameError:
                # if unsuccessful return error
                self.log.warn('TableMergeRows: Error in %s( %s ) - '
                    % (funct, colname) + ' returning first row value')
        # return output row
        self.logcount += 1
        if self.logcount%10 == 0:
            self.log.debug('TableMergeRows: done')
        return outrow

    def tablemergetables(self,tables):
        """Returns a new table containing data merged from the input table(s).
           Columns are merged according to the rules define in the [table]
           section of the configuration file.

           tables is a list of tables, previously obtained with something
           like tableget().  Note that each table is assumed to have a single
           row of data.  If you need to merge rows of data, use something
           like tablemergerows().

        """

        ntable = len(tables)
        self.log.debug("ntable = %d" %ntable)
        names   = tables[0].names # list of fields
        formats = tables[0].columns.formats
        dims    = tables[0].columns.dims
        units   = tables[0].columns.units
        # loop through all tables, make sure names, formats, dims, and units
        # are the same among them.
        for i in xrange(1,ntable):
            if cmp(names,tables[i].names): # equals zero if identical
                self.log.error('column names differ for merging')
                raise ValueError
            if cmp(formats,tables[i].columns.formats): # equals zero if identical
                self.log.error('column formats differ for merging')
                raise ValueError
            if cmp(dims,tables[i].columns.dims): # equals zero if identical
                self.log.error('column names differ for merging')
                raise ValueError
            if cmp(units,tables[i].columns.units): # equals zero if identical
                self.log.error('column names differ for merging')
                raise ValueError

        cols = []
        for n,f,d,u in zip(names,formats,dims,units):
            # get merge function
            try:
                funct = self.config['table'][n.lower()].lower()
            except :
                # if keyword is not available
                self.log.warn('TableMergeTables: Missing table merge entry for '
                    + 'column -%s- returning first row value' % n)
                funct = 'first'
            if funct == 'first':
                tmp = tables[0][n]
            elif funct == 'last':
                tmp = tables[-1][n]
            elif funct == 'min':
                tmp = numpy.nanmin([a[n] for a in tables])
            elif funct == 'max':
                tmp = numpy.nanmax([a[n] for a in tables])
            elif funct == 'med': # median
                tmp = numpy.nanmedian([a[n] for a in tables])
            elif funct == 'avg': # average
                tmp = numpy.nanmean([a[n] for a in tables])
            elif funct == 'sum':
                tmp = numpy.nansum([a[n] for a in tables])
            elif funct == 'wtavg':
                tmp = numpy.nansum([a[n]*a['Samples'] for a in tables])
                tmpsum = numpy.nansum([a['Samples'] for a in tables])
                tmp = tmp/float(tmpsum)
            else:
                self.log.warn('TableMergeTables: Unknown operation -'+funct
                    +'- for column -'+n+'- returning first row value')
                tmp = tables[0][n]
            cols.append(fits.Column(name=n, format=f, dim=d, unit=u,
                array=[tmp]))

        tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        return tbhdu

    def getheader(self, dataname = ''):
        """ Returns the stored header. If dataname is specified, the
            header of the named image/table is returned, otherwise the
            primary header.
        """
        # Return primay header
        if dataname == '':
            header = self.header
        # Return image header
        elif dataname.upper() in self.imgnames:
            index = self.imageindex(dataname)
            header = self.imgheads[index]
        # Return table header
        elif dataname.upper() in self.tabnames:
            index = self.tableindex(dataname)
            header = self.tabheads[index]
        # Dataname not found -> return error message
        else:
            msg = 'Invalid data name (%s)' % dataname
            self.log.error('GetHeader: ' + msg)
            raise ValueError(msg)
        self.log.debug('GetHeader: done')
        return header

    def setheader(self,header, dataname = ''):
        """ Overwrites the stored header. If dataname is specified, the
            header of the named image/table is set, otherwise the
            primary header.
        """
        # Set primay header
        if dataname == '':
            self.header = header
        # Set image header
        elif dataname.upper() in self.imgnames:
            index = self.imageindex(dataname)
            self.imgheads[index] = header
        # Set table header
        elif dataname.upper() in self.tabnames:
            index = self.tableindex(dataname)
            self.tabheads[index] = header
        # Dataname not found -> return error message
        else:
            msg = 'Invalid data name (%s)' % dataname
            self.log.error('SetHeader: ' + msg)
            raise ValueError(msg)
        self.log.debug('SetHeader: done')


    def getheadval(self, key, dataname = '', errmsg=True):
        """ Get Header Value: Returns the value of the requested key from
            the header. If the keyword is present in the [Header] section
            of the configuration that value is returned instead. In case that
            value from the configuration file is itself a header key, the value
            stored under that key is returned. If the key can not be found an
            KeyError is produced and a warning is issued.
        """
        val = None
        inkey = key # retain key which was input in case key changes
        # Look in the config
        try:
            # get the value
            val = self.config['header'][key.upper()]
            # Check if it's optional header replacement i.e. starts with '?_'
            if val[:2] in ['?_', '? ', '?-']:
                # if key is not in the header -> use key name under value instead
                if key not in self.header:
                    self.log.info('Getheadval: Using %s keyword for %s' %
                                  (val[2:].upper(), key))
                    key = val[2:].upper()
                val = None
            # Check if it's a Header replacement (but not T/F)
            elif val[0].isupper() and val[:2] not in ['T ', 'F ']:
                self.log.info('Getheadval: Using %s value for %s' %
                              (val.upper(), key ))
                key = val.upper()
                val = None
            else:
                # make it a pyfits.Card then get value and comment
                card = fits.Card()
                card = card.fromstring(key.upper() + ' = ' + val)
                # Following code is unused because non-string values
                #   don't get converted correctly. However the new code
                #   may not work with all version of pyfits - tbd - nlc 2013
                #tmp = val.split('/')
                #if len(tmp) > 1:
                #    card=pyfits.Card(key,tmp[0].strip(),tmp[1].strip())
                #else:
                #    card=pyfits.Card(key,tmp[0].strip())
                self.log.info('Getheadval: Setting %s to %s' % (key, val))
                val=card.value
                # update value in header
                self.setheadval(key, val, card.comment)
        except KeyError :
            # if key is not in config - continue
            pass
        except TypeError :
            # if config is not yet loaded - return error
            self.log.warn('GetHeadVal: Missing Configuration')
            # The following line is commented out such that pipedata can
            #     still be used without configuration
            #raise RuntimeError('Missing Pipe Configuration')
        # Look in the header
        if val is None:
            # get the header from dataname
            header = self.getheader(dataname)
            # get value from header
            try:
                val = header[key]
            except KeyError:
                # if keyword is not found
                msg='Missing %s keyword in header %s' % (key, dataname)
                if errmsg:
                    self.log.error('GetHeadVal: %s' % msg)
                raise KeyError(msg)
        # if Value is a pyfits.core.Undefined i.e. no keyword
        if isinstance(val, fits.Undefined):
            msg = 'Missing value for key = %s - returning empty string' % key
            self.log.warn('GetHeadVal: %s' % msg)
            val = ''
        # Debug message (only first line if multi-line value)
        if '\n' in str(val):
            self.log.debug('GetHeadval: done (%s=%s ...)' % (inkey, str(val[0])))
        else:
            self.log.debug('GetHeadVal: done (%s=%s)' % (inkey,str(val)))
        return val

    def setheadval(self, key, value, comment=None, dataname = ''):
        """ Sets a FITS keyword in the current header
        """
        # If no header exitst, make a first empty image
        if len(self.imgheads) == 0:
            self.imgheads = [fits.PrimaryHDU().header]
            self.imgdata = [None]
            self.imgnames = ['PRIMARY']
        # Set primay header
        if dataname == '':
            hdr = self.header
        # Set image header
        elif dataname.upper() in self.imgnames:
            index = self.imageindex(dataname)
            hdr = self.imgheads[index]
        # Set table header
        elif dataname.upper() in self.tabnames:
            index = self.tableindex(dataname)
            hdr = self.tabheads[index]
        # Dataname not found -> return error message
        else:
            msg = 'Invalid data name (%s)---' % dataname
            self.log.error('SetHeadVal: ' + msg)
            raise ValueError(msg)
        # Set value into hdr
        if comment != None: hdr[key] = (value,comment)
        else: hdr[key] = value
        self.log.debug('SetHeadVal: done')

    def delheadval(self, key, dataname = None):
        """ Delete one or more FITS keyword in specified header, which defaults
            to the first header.
        """

        if isinstance(key,(list,tuple)):
            for k in key:
                self.delheadval(k,dataname)
        elif isinstance(key,str):
            if dataname is None:
                try:
                    del(self.header[key])
                except KeyError: # don't care if it doesn't exist
                    pass
            elif dataname.upper() in self.imgnames:
                index = self.imageindex(dataname)
                try:
                    del(self.imgheads[index][key])
                except KeyError: # don't care if it doesn't exist
                    pass
            elif dataname.upper() in self.tabnames:
                index = self.tableindex(dataname)
                try:
                    del(self.tabheads[index][key])
                except KeyError: # don't care if it doesn't exist
                    pass
            else: # Dataname not found -> return error message
                msg = 'Invalid data name (%s)---' % dataname
                self.log.error('DelHeadVal: ' + msg)
                raise ValueError(msg)
        else:
            msg = 'Invalid key (%s).  Must be a str, list, or tuple' %repr(key)
            self.log.error('DelHeadVal: ' + msg)
            raise ValueError(msg)
        self.log.debug('DelHeadVal: done')

    def test(self):
        """ Test Pipe Data Object: Runs a set of tests on a pipedata object.
            This function needs to start in the folder contining the test
            files.
        """
        self.log.info("Testing Data Object")
        ### Setup
        # test loading config
        #self.setconfig(1) # should raise TypeError
        self.setconfig(None)
        conf = self.testconf
        self.setconfig(conf)
        conf=configobj.ConfigObj(conf)
        self.setconfig(conf)
        ### Loading and Saving DAta
        # self.load(1) # should create TypeError
        # self.load('nothere.fits') # should create IOError
        # self.load('makebasic.pro') # should create IndexError
        testpath = self.config['testing']['testpath']
        print(testpath+'/makedata/testimgtabimg.fits')
        self.load(testpath+'/makedata/testimgtabimg.fits')
        self.log.debug('Images = ' + repr(self.imgnames))
        self.save(testpath+'/testsave.fits')
        print(self.filenamebegin, self.filenameend)
        #self.save(1) # creates TypeError
        ### Copy, Merge and Add Data
        other=self.copy()
        other.mergedata(self)
        other.load(testpath+'/makedata/testimgtabimg.fits')
        other.imgdata[1]=numpy.zeros((10,20))
        #other.table = numpy.rec.array([(1)], names='nr', formats='i2')
        self.mergedata(other)
        self.adddata({'First Image':numpy.zeros((10,10))},['Delta',4])
        for i in range(len(self.imgdata)):
            s='%s %s' % (self.imgnames[i],self.imgdata[i].shape)
            self.log.debug(s)
        self.log.debug(repr(self.table))
        ### Images / Table / Header
        # Images
        img=self.imageget('SECOND IMAGE')
        self.imageset(numpy.zeros((5,5,5)),'First Image Copy')
        self.imagedel('SECOND IMAGE')
        self.imageset(img,'Second Image',index=-1)
        other.imgnames[0]='WILDCARD'
        self.copydata(other,'WILDCARD')
        # Image variable
        self.image=numpy.zeros((5,5))
        self.image[0,:]=range(5)
        # Table
        i=len(self.table)
        #print(self.table,self.table.dtype,self.table.shape)
        self.tableaddcol('index',numpy.arange(i)-2)
        self.tableaddrow(['Zeta',0,100])
        self.tabledelcol('NR')
        self.tabledelrow(3)
        #print(self.table,self.table.dtype,self.table.shape)
        self.table['index']=numpy.arange(len(self.table))+2
        self.table[0]=('first',100)
        self.tableaddcol('arr', numpy.ones((7,2)))
        self.log.debug(repr(self.table))
        self.log.debug(repr(self.table.dtype))
        self.log.debug('TabMrg='+repr(self.tablemergerows(self.table)))
        # Second table
        tab2 = self.tableget()
        self.tableset(tab2,'Second Table')
        print(self.tableget('Second Table').dtype,self.tabdata[1])
        self.tabledelcol('arr','Second Table')
        self.tabledelrow(6,'Second Table')
        self.tabledelrow(5,'Second Table')
        self.tableaddrow(['velo',19], 'Second Table')
        self.tableaddcol('NewNr', range(6), 'Second Table')
        print(self.tableget('Second Table').dtype,self.tabdata[1])
        # Delete first table and save
        self.tabledel('TABLE')
        self.save(testpath+'/testsave.fits')
        # Header
        self.header['Reduced'] = 0
        self.setheadval('TST','Test Value','') #,'Second Table')
        print(self.header)
        self.getheadval('TST') #,'Second Table')
        self.setheadval('TST',15)
        print(self.header)
        self.setheader(self.getheader())
        # Additional Load and Save Tests
        if os.path.exists(testpath+'/makedata/outtest.fits'):
            os.remove(testpath+'/makedata/outtest.fits')
        self.save(testpath+'/makedata/outtest.fits')
        #self.save('pipedata/outtest.fits')
        #self.load('pipedata/outtest.fits')
        #self.load('sample.chop_nod.dog.raw.fits')
        #self.load('sample.noise.fft.raw.fits')
        self.log.info("Test Finished")
        print('Done')
        return 0

if __name__ == '__main__':
    """ Main Function for testing:
        Changes to the folder with the test files and runs pipedata.test()
    """
    logging.basicConfig(level=logging.DEBUG)
    data=DataFits()
    data.test()

""" === History ===
    2016-9-26 Marc Berthoud: Split PipeData (fits only) into
        DataParent and DataFits
        Main Changes:
        - Added data attribute which can go to image or table
    2015-9-25 Marc Berthoud: Edited save() to update keywords for
        FILENAME and DATE
    2014-7-31 Marc Berthoud: Added output_verify=fix to write()
        Fixed tablemergerows to work for FITS_rec
        Now use numpy.insert for tableaddrow, numpy.delete in tabledelrow
    2014-2-17 Marc Berthoud: Fix error with load if primary header empty
    2012-11-15 Marc Berthoud: Loadhead now loads from other HDUs
    2012-9-27 Marc Berthoud: Add functionality for header keyword replacement
    2012-9-13 Marc Berthoud: Allow files to have empty HDUs
        i.e. image.shape=(0,)
    2012-5-17 Marc Berthoud: Added storate and management of headers for all
    2012-2-17 Marc Berthoud: Added multiple tables option
        Many new commands (tableindex/get/set/del) and upgrades to table
        commands.
        COMMENT: adddata now only affects first table (may have to change)
    2011-2-4 Marc Berthoud: Added filenamebase to __getattr__
    2010-12-20 Marc Berthoud: Made sure that all image names would always be
        uppercase.
    2010-12-3 Marc Berthoud: Added code to tableaddcol such that arrays
        can be stored in tables.
    2010-10-1  Marc Berthoud, Ver0.1.1: New structure for images:
        - Now there can be multiple images
        - The size of the primary image doesn't have to agree with the
          table size
    2009-11-17 Marc Berthoud, Added remap to sky index
    2008-10-28 Marc Berthoud, Ver0.1: Wrote and tested
"""
