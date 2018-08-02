""" PIPE DATA - Version 1.0.0

    This module defines the pipeline data parent object. The object
    has limited use to store, load and save data, it's child objects
    should be used instead.

    If the configuration file is properly populated, calling DataParent.load() or
    DataParent.loadhead() will return the proper DataChild object with the
    loaded file.

    The module also contains crucial variables for the entire pipeline,
    such as pipever and testconf.

"""

import numpy # numpy library
import os # os library
import logging # logging library
import time # time library
import configobj # config object library
import validate  # to check input config file is correct
import re # regexp
import gc # garbage collect
from __builtin__ import True
from fileinput import filename

class DataParent(object):
    """ HAWC Pipeline Data Object
        This object stores a header and the config file. Potential data is stored
        as a single image (numpy array).
    """
    # General variables: These are valid for all pipeline and all pipesteps
    pipever = '1.3.0beta3' # Pipeline version

    #testconf = 'config/pipeconf_master.txt' # Test configuration
    testconf = 'config/pipeconf_mgb.txt' # Test configuration

    # File Name Fit: Regexp expression that fits filenames
    filenamefit = '\A\Z' # Empty filename, only fits empty string

    def __init__(self, config=None):
        """ Constructor: Initialize data object and variables
        """
        # set up internal variables
        #     Remark: to preserve order of images no dictionary is used
        self.filename='' # Name of loaded file
        # Data Variable:
        self.data = None # Empty data array
        # Header: A dictionary. Lists for HISTORY and COMMENT entries.
        self.header = {}
        # set up logging and configuration
        self.log=logging.getLogger('pipe.data') # Logger
        self.config = None # Configuration
        self.setconfig(config) # Set the Configuration
        # Log message
        self.log.debug('Init: done')

    def __getattr__(self, name):
        """ Get Attribute function that allows to access these attributes:
            - filenamebegin: filename start. based on data:filenamebegin in
                             the configuration file, defaults to all before
                             .something.something.
            - filenameend: filename end. based on data:filenameend in the
                           configuration file, defaults to .something at the
                           end of the file.
            - filenum: file number. based on data:flilenum in the configuration
                       file. Returns None if there is no file number.
        """
        # return filenamebegin if it's requested
        if name == 'filenamebegin':
            (fpath,fname) = os.path.split(self.filename)
            # if filenamebegin is specified, use it
            try:
                filenamebegin = self.config['data']['filenamebegin']
                match = re.search(filenamebegin,fname)
            except (KeyError, TypeError):
                match = False
                filenamebegin = "None"
            if match :
                # found file name beginning -> return it
                return os.path.join(fpath,match.group())
            else:
                # assume filename format is name.filestep.fits or ....fits.gz
                msg = "Filename=%s doesn't match pattern=%s" % (fname, filenamebegin)
                self.log.warn(msg)
                extloc = fname.rfind('.f')
                #print extloc,fname
                if ( extloc < 0 or
                     fname[extloc:] not in ['.fts','.fits','.fits.gz'] ) :
                    self.log.warn('Filename has non-fits extension')
                    extloc = fname.rfind('.')
                    if extloc < 0: extloc = len(fname)
                    else: extloc += 1 # to add the '.'
                else: extloc += 1 # to add the '.'
                typeloc = fname[0:extloc-1].rfind('.')
                if typeloc < 0:
                    typeloc = extloc
                else:
                    typeloc += 1
                return os.path.join(fpath,fname[:typeloc])
        # return filenameend if it's requested
        if name == 'filenameend':
            (fpath,fname) = os.path.split(self.filename)
            # if filenameend is specified, use it
            try:
                filenameend = self.config['data']['filenameend']
                match = re.search(filenameend, fname)
            except (KeyError, TypeError):
                match = False
                filenameend = "None"
            if match :
                # found file name end -> return it
                return match.group()
            else:
                msg = "Filename=%s doesn't match pattern=%s" % (fname, filenameend)
                extloc = fname.rfind('.f')
                if ( extloc < 0 or
                     fname[extloc:] not in ['.fts','.fits','.fits.gz'] ) :
                    self.log.warn('Filename has non-fits extension')
                    extloc = fname.rfind('.')
                    if extloc < 0: extloc = len(fname)
                if extloc < 0:
                    return ''
                else:
                    return fname[extloc:]
        # return file number if it's requested
        if name=='filenum':
            (fpath,fname) = os.path.split(self.filename)
            # if filenum is specified in pipeconf, use it
            try:
                filenum = self.config['data']['filenum']
                match = re.search(filenum, fname)
            except (KeyError, TypeError):
                match = False
                filenum = None
            if match :
                # found file num -> return it
                for m in match.groups():
                    if m is not None:
                        return m
                return match.groups()
            else:
                return None
        # raise error if attribute is unknown
        msg="'%s' object has no attribute '%s'" % (type(self).__name__,name)
        raise AttributeError(msg)

    def setconfig(self,config):
        """ Sets configuration for the pipe data: The configuration
            object is returned. The config variable can be one of these
            - A ConfigObj object
            - A string to the filename of a valid config file
            - A list with strings to valid config files. In this case
              the first file is loaded, the others are merged with it.
              If all files have the same folder, only the first
              filepathname needs the full path.
        """
        if isinstance(config,configobj.ConfigObj):
            # if config is a ConfObj -> set it
            self.config=config
            self.log.debug('SetConfig: skipping configuration file validation')
            retmsg='received ConfigObj'
        elif isinstance(config,basestring):
            # if config is a string - check for file existence -> load it
            if os.path.isfile(config):
                try:
                    valpath = os.path.dirname(__file__)
                    valpath = os.path.join(valpath,'conf_validate.txt')
                    valpath = ''
                    self.config=configobj.ConfigObj(config,configspec=valpath)
                    #print self.config.configspec
                    #self.validateconfig() # taken out for now
                    # configuration validation couldn't keep up with
                    # new configuration parameters
                except configobj.ConfigObjError, error:
                    msg = 'Error while loading configuration file'
                    self.log.error('SetConfig: ' + msg)
                    raise error
                retmsg='loaded config file=' + config
            else:
                msg='<%s> is invalid file name for configuration' % config
                self.log.error('SetConfig: '+msg)
                raise IOError(msg)
        elif config is None:
            # no configuration -> No change
            # (( The following code would have been to use default config
            #    from the config validation file ))
            #self.config = None
            #path = os.path.dirname(__file__)
            #path = os.path.join(path,'conf_validate.txt')
            #self.config = configobj.ConfigObj(configspec=path)
            #self.config.validate(validate.Validator(),copy=True)
            #write out default
            #fp = open("pipe_test_conf.txt",'w')
            #self.config.write(fp)
            #fp.close()
            retmsg='no config given'
        else: # Assume it's a list-like object of config files
            # Get first object and check if it's really a list
            try:
                baseconf = config[0]
            except: # If it's not list-ish 
                # Invalid configuration - error
                self.log.error('SetConfig: Invalid configuration variable = %s' % repr(config))
                raise TypeError('Invalid configuration variable = %s' % repr(config))
            # Load first config
            self.setconfig(baseconf)
            basefolder = os.path.split(baseconf)[0]
            # Load other configs on top
            for conf in config[1:]:
                # Check if file exists
                if os.path.isfile(conf):
                    self.mergeconfig(conf)
                else:
                    cnf = os.path.join(basefolder, os.path.split(conf)[1])
                    if os.path.isfile(cnf):
                        self.mergeconfig(cnf)
                    else:
                        self.log.error('SetConfig: Invalid configuration file in list = %s' % conf)
                        raise ValueError('Invalid configuration file in list = %s' % conf)
            retmsg  = 'config list with %d files' % len(config)
        # Set environment variables from [envars]
        if not self.config is None:
            if self.config.has_key('envars'):
                for var in self.config['envars']:
                    os.environ[var] = str( self.config['envars'][var] )
        # Return
        self.log.debug('SetConfig: done ('+retmsg+')')
        return self.config

    def mergeconfig(self, newconfig):
        """ Merges a configuration object (or file) into the existing
            configuration. All values from the new configuration are used,
            overwriting old values if they are already in the old
            configuration.
        """
        # If there is no existing config: Just load the new config
        if not self.config:
            self.setconfig(newconfig)
            return
        # Store the existing config and load new config
        oldconf = self.config
        newconf = self.setconfig(newconfig)
        # Merge new to old config and store it
        oldconf.merge(newconf)
        self.config = oldconf
        
    def validateconfig(self):
        """ Test config against configspec and print errors if it doesn't
            conform.
            
            This code is currently not used

            Things that need to be fixed:
            * Config validation should not overwrite default value from
              paramlist in steps
            * If a keyword is missing in the [HEADER] section config validation
              shouldn't put it in.
        """

        errFlag = False # set to true if errors encountered
        results = self.config.validate(validate.Validator(),copy=True)
        extra   = configobj.get_extra_values(self.config) # keywords not in spec
        # List configuration keywords that failed validation
        if results != True:
            for (section_list, key,_) in configobj.flatten_errors(self.config,
                                                                  results):
                blah = ', '.join(section_list)
                if key is not None:
                    msg = "ValidateConfig: key '%s' in section '%s' failed"
                    self.log.error( msg % (key, blah))
                    errFlag = True
                else:
                    self.log.error("ValidateConfig: Section '%s' failed" %blah)
                    errFlag = True
        # Warn for keywords not found in validation
        for s,k in extra:
            if len(s) == 0:
                msg = 'ValidateConfig: Skippping unknown global keyword "%s"'
                self.log.warning(msg % k)
            else:
                msg = 'ValidateConfig: Skippping unknown keyword "%s" in section "%s"'
                self.log.warning( msg %(k,s[0]))
        if errFlag:
            raise validate.ValidateError
        self.log.debug('ValidateConfig: done')

    def getobject(self,objname):
        """ Looks for a module with the listed object. This function
            searches all package entries in config[general][steppacks] for
                package.objname
            If the module loads, it makes an objname object and returns it.
        """
        # get list of packages
        try:
            steppacks = self.config['general']['steppacks']
            if isinstance(steppacks,str): # ensure we have a list if only
                steppacks = [steppacks]    # 1 item, steppacks is str
        except KeyError, error:
            self.log.error('Setup: Missing steppacks item in configuration')
            raise error

        # Load object module
        stepmodule = None
        for pack in steppacks:
            if pack == '.': # look in local directory
                mod = objname.lower()
            else:
                mod = '%s.%s' % (pack,objname.lower())
            self.log.debug('looking for %s' % mod)
            # import the module
            try:
                stepmodule = __import__(mod, globals(), locals(),
                                        [objname])
                self.log.debug('Pipe step %s found in %s' %
                              (objname,pack))
                break
            except ImportError,msg:
                tmp = 'No module named %s' % objname.lower()
                if str(msg).startswith(tmp): # module not present in directory
                    self.log.debug('Pipe object %s not found in %s' %
                                   (objname, pack))
                else: # module tries to import a missing package
                    raise
            except:   # print out import errors not due to missing
                raise # modules (e.g., typos in code)
        # If not found -> Raise error
        if stepmodule == None:
            msg = 'Could not find object %s' % objname
            self.log.error('GetObject: ' + msg)
            raise ImportError(msg)
        # Make an object instance
        try:
            retobj = stepmodule.__dict__[objname]()
        except KeyError, error:
            msg = 'Pipe object %s' % objname
            msg+= ' not found in module %s' % mod
            self.log.error('GetObject: %s' % msg)
            raise error
        # Return the object
        return retobj

    def datamatch(self,filename=''):
        """ Returns an instance of the first data object in
            config['data']['dataobjects'] which matches the file under
            filename. This function returns an error if the file doesn't
            exist. The data is not loaded into the returned object.
        """
        # Get list of data objects
        try:
            dataobjects = self.config['data']['dataobjects']
            if isinstance(dataobjects,str): # ensure we have a list if only
                dataobjects = [dataobjects]    # 1 item, steppacks is str
        except KeyError, error:
            self.log.error('Setup: Missing dataobjects item in configuration')
            raise error
        # Loop through dataobjects
        found = False
        for dataname in dataobjects:
            # Get an object
            dataobj = self.getobject(dataname)
            dataobj.setconfig(self.config)
            # See if file fits
            if re.search(dataobj.filenamefit, filename) is not None:
                found = True
                break
        if not found:
            msg = 'DataMatch: No matching data object in %s for filename %s' % (repr(dataobjects), filename)
            self.log.error(msg)
            raise ValueError(msg)
        else:
            msg = 'DataMatch: Found matching data object %s for filename %s' % (dataname, filename)
            self.log.debug(msg)
        # Return the object
        return dataobj

    def loadhead(self,filename=''):
        """ Searches for a matching dataobject for the file given, then
            calls the loadhead() function of that object on the file.
            The new dataobject is returned.
            !! Only DataParent.loadhead returns such an object,     !!
            !! no return is expected from loadhead of child objects !!
        """
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # Check for file existance
        if not os.path.exists(filename):
            msg = 'LoadHead: %s is invalid filename' % filename
            self.log.error(msg)
            raise ValueError(msg)
        # Create a data object for the file
        dataobj = self.datamatch(filename)
        # Load the header for the parent object
        dataobj.loadhead(filename)
        # Return the parent object
        return dataobj

    def load(self,filename=''):
        """ Searches for a matching dataobject for the file given, then
            calls the load() function of that object on the file.
            The new dataobject is returned.
            !! Only DataParent.load returns such an object,     !!
            !! no return is expected from load of child objects !!
        """
        # set self.filename and filename
        if len(filename) > 0:
            self.filename=filename
        else:
            filename=self.filename
        # Check for file existance
        if not os.path.exists(filename):
            msg = 'Load: %s is invalid filename' % filename
            self.log.error(msg)
            raise ValueError(msg)
        # Create a data object for the file
        dataobj = self.datamatch(filename)
        # Load the header for the parent object
        dataobj.load(filename)
        # Return the parent object
        return dataobj

    def save(self,filename = ''):
        """ Save the data in the object to the specified file. Existing files are
            overwritten. This issues a warning (dataparent should not be used for
            storing data) then searches for a matching dataobject to the given
            filename. Data and header are copied to that dataobject and the data
            is saved.
            - filename: a string with the file name to save to
        """
        # Issue warning
        self.log.warn('DataParent.save should not be used as DataParent is not ' +
                      'intended to store data.')
        # get file name
        if len(filename) > 0:
            self.filename = filename
        else:
            filename = self.filename
        # Create data object for the file
        dataobj = self.datamatch(filename)
        # Set filename
        dataobj.filename = filename
        # Copy header
        for key in self.header.keys():
            if key != 'COMMENT' and key != 'HISTORY':
                dataobj.setheadval(key,self.header[key])
        # Copy Comments
        if 'COMMENT' in self.header:
            for com in self.header['COMMENT']:
                dataobj.setheadval('COMMENT',com)
        # Copy History
        if 'HISTORY' in self.header:
            for hist in self.header['HISTORY']:
                dataobj.setheadval('HISTORY',hist)
        # Copy data
        dataobj.data = self.data
        # Save
        dataobj.save()

    def copy(self):
        """ Returns a copy of self
        """
        out=DataParent(config=self.config) # create new object
        # copy filename and header
        out.filename = self.filename
        out.header = self.header.copy()
        # Copy data - backup if no copy() available
        try:
            out.data = out.data.copy()
        except:
            out.data = out.data
        # return message and new object
        self.log.debug('Copy: done')
        return out

    def mergehead(self, other):
        """ Merges the header of another data object to the existing header.
            Most of the header of the new data is ignored, PARENT,
            HISTORY keywords are copied.
        """
        # get selfhist and otherhist lists
        if 'HISTORY' in self.header:
            selfhist = self.header['HISTORY']
        else: selfhist = []
        if 'HISTORY' in other.header:
            otherhist = other.header['HISTORY']
        else: otherhist = []
        # add history keywords (no duplicates)
        selfhist += [hist for hist in otherhist if not hist in selfhist]
        # if there is something add write back to header
        if len(selfhist) : self.header['HISTORY'] = selfhist
        # get selfcomm and othercomm lists
        if 'COMMENT' in self.header:
            selfcomm = self.header['COMMENT']
        else: selfcomm = []
        if 'COMMENT' in other.header:
            othercomm = other.header['COMMENT']
        else: othercomm = []
        # add comment keywords (no duplicates)
        selfcomm += [comm for comm in othercomm if not comm in selfcomm]
        # if there is something add write back to header
        if len(selfcomm) : self.header['COMMENT'] = selfcomm

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

    def getheadval(self, key, errmsg=True):
        """ Get Header Value: Returns the value of the requested key from
            the header. If the key is present in the [Header] section
            of the configuration that value is returned instead, the following
            entries are possible in the configuration file:
              * KEY = VALUE * VALUE is returned, the system checks if value is an
                              int or a float, else a string is returned.
              * KEY = NEWKEY * The VALUE under header[NEWKEY] is returned.
              * KEY = ?_ALTKEY * If the keyword KEY is present, header[KEY] is
                                 returned, else header[ALTKEY] is returned.
            If the key can not be found an KeyError is produced and a warning is
            issued (unless key is present in the [Header] section of the
            configuration).

            errmsg: Flag indicating if a log error message should be
                    issued if the keyword is not found. This can be disabled
                    (set it to False) if getheadval is used to probe a dataset.
        """
        val = None
        inkey = key # retain key which was input in case key changes
        # Look in the config
        try:
            # get the value
            val = self.config['header'][key]
            # Check if it's optional header replacement i.e. starts with '?_'
            if val[:2] in ['?_', '? ', '?-']:
                # if key is not in the header -> use key name under value instead
                if not key in self.header:
                    key = val[2:].upper()
                val = None
            # Check if it's a Header replacement (but not T/F)
            elif val[0].isalpha() and val[:2] not in ['T ', 'F '] and val not in ['T', 'F']:
                self.log.info('Getheadval: Using %s value for %s' %
                              (val.upper(), key ))
                key = val.upper()
                val = None
            # Else: read value
            else:
                # Try as T / F
                found = True
                if val == 'T' or val[:2] == 'T ':
                    val = True
                elif val == 'F' or val[:2] == 'F ':
                    val = False
                else:
                    found = False
                # Try as int
                if not found:
                    try:
                        val = int(val)
                        found = True
                    except:
                        pass
                # Try as float
                if not found:
                    try:
                        val = float(val)
                        found = True
                    except:
                        pass
                # If not found - just leave value as string
                # update value in header
                self.setheadval(key, val)
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
            # get value from header
            try:
                val = self.header[key]
            except KeyError:
                # if keyword is not found
                msg='Missing %s keyword in header' % key
                if errmsg:
                    self.log.error('GetHeadVal: %s' % msg)
                raise KeyError(msg)

        # Debug message (only first line if multi-line value)
        if '\n' in str(val):
            self.log.debug('GetHeadval: done (%s=%s ...)' % (inkey, str(val[0])))
        else:
            self.log.debug('GetHeadVal: done (%s=%s)' % (inkey,str(val)))
        return val

    def setheadval(self, key, value, comment = ''):
        """ Sets a FITS keyword in the header
        """
        # If key==HISTORY or COMMENT: add to history list
        if key=='HISTORY' or key=='COMMENT':
            if key in self.header:
                self.header[key].append(value)
            else:
                self.header[key] = [value,]
        # Else: add as normal keyword
        else:
            self.header[key] = value
            if len(comment) > 0:
                self.setheadval('COMMENT','%s, %s' % (key,comment))
        self.log.debug('SetHeadVal: done')

    def delheadval(self, key):
        """ Delete one or more FITS keyword in specified header, which defaults
            to the first header.

            WARNING: This will delete all COMMENT or HISTORY entries
        """
        # If key is a list, remove all entries
        if isinstance(key,(list,tuple)):
            for k in key:
                self.delheadval(k)
        # Else if it's a string delete the key - this may rise KeyError
        else:
            del(self.header[key])
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
    data=DataParent()
    data.test()

""" === History ===
    2016-9-26 Marc Berthoud: Split PipeData (fits only) into
        DataParent and DataFits
        Main Changes:
        - init no longer accepts a file name (since it can't return
          the correct object.
        - each object has a filenamefit variable (is fitted with regexp)
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
