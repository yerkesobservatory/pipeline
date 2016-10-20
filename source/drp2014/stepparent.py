""" PIPE STEP PARENT - Version 1.0.0

    This module defines the HAWC pipeline step parent object. Pipe steps are
    the modules responsible for all HAWC data reduction. They are called by
    the pipeline and work with pipedata objects. All pipe step objects are
    descendants from this one. Pipe steps are callable objects that return
    the reduced data product (as pipedata object).
    
    @author: berthoud
"""

"""
    2DO:
    - Convert all steps to new format and test them

    - later:
      - rename callstart/end -> runstart/end
      - remove parameter from callstart/end
"""

import time # time library
import logging # logging object library
import argparse # Argument parsing library
from configobj import ConfigObj # configuration object
from drp.pipedata import PipeData # pipeline data object

class StepParent(object):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    
    stepver = '0.1.1' # pipe step version
    #testconf = 'config/pipeconf_master.txt' # Test configuration
    testconf = 'config/pipeconf_mgb.txt' # Default test configuration
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
            calls the setup function.
        """
        # initialize input and output
        self.datain = PipeData()
        self.dataout = PipeData()
        # set names
        self.name = None
        self.procname = None
        # set parameters / logger
        self.arglist={} # Dictionary with current arguments
        self.paramlist=[] # List with possible parameters
        # set configuration / logger
        self.config = None
        self.log = None
        # specify whether this step runs on a single PipeData object with a 
        # single output PipeData object (SISO), multiple input PipeData objects
        # with multiple output PipeData objects (MIMO), or multiple input Pipefile
        # objects with a single output PipeData object (MISO).
        self.iomode = 'SISO'
        # do local setup
        self.setup()
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
        self.name='parent'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['sampar', 1.0,
            'Sample Parameter - parent only - no practical use'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        # Log the value of sample parameter
        self.log.debug("Sample Parameter = %.2f" % self.getarg('sampar'))
        # Copy datain to dataout
        self.dataout = self.datain
        # Set complete flag
        self.dataout.setheadval('COMPLETE',1,
                                'Data Reduction Pipe: Complete Data Flag')


    def __call__(self,datain, **arglist):
        """ Object Call: returns reduced input data
        """
        # Get input data
        self.datain = datain
        # Start Setup
        self.runstart(self.datain, arglist)
        # Call the run function
        self.run()
        # Finish - call end
        self.runend(self.dataout)
        # return result
        return self.dataout
    
    def runstart(self, data, arglist):
        """ Method to call at the beginning of the pipe step call.
            - Sends initial log message
            - Checks the validity of the input data
            - Gets configuration from input data and checks type
        """
        # Start Message
        self.log.info('Start Reduction: Pipe Step %s' % self.name)
        # Set input arguments
        for k in arglist.keys():
            self.arglist[k.lower()] = arglist[k]
        # Check input data type and set data config
        if isinstance(data,PipeData):
            self.config = data.config
        else:
            msg = 'Invalid input data type: PipeData object is required'
            self.log.error(msg)
            raise TypeError('Runstart: '+msg)
        # Set Configuration
        if self.config is None: # no config specified, make an empty one
            self.config = ConfigObj()
            self.config[self.name] = {}
        # Check configuration
        if not isinstance(self.config,ConfigObj):
            msg='Invalid configuration information - aborting'
            self.log.error(msg)
            raise RuntimeError('Runstart: '+msg)
        
    def runend(self,data):
        """ Method to call at the end of pipe the pipe step call
            - Sends final log messages
        """
        # clear input arguments
        self.arglist = {}
        # update header (status and history)
        self.updateheader(data)
        self.log.info('Finished Reduction: Pipe Step %s' % self.name)

    def updateheader(self,data):
        """ Update the header for a single PipeData object
            - Sets the PROCSTAT and PROCLEVL keywords in the data header
            - Adds a history entry to the data header
        """
        
        data.setheadval('PROCSTAT','Intermediate','Processing Status')
        data.setheadval('PROCLEVL',self.name,'Processing Level')
        histmsg = 'Reduced: ' + self.name + ' v' + self.stepver + ' '
        histmsg += time.strftime('%Y-%m-%d_%H:%M:%S')
        data.header.add_history(histmsg)
        # update file name with .PipeStepName.fits
        data.filename = data.filenamebase + '.' + self.procname + '.fits'
        # Send log messages

    def getarg(self, parname):
        """ Returns the argument value of the parameter parname. The parameter
            is first searched for in self.arglist['parname'], then in
            config['stepname']['parname']. If the parameter is not found,
            the default value from parameter list is returned.
            Should the parameter name not have an entry in the parameter list,
            a error is returned and a KeyError is raised.
        """
        # list of strings that should parse to boolean true
        # we need to handle booleans separately, because bool("False")
        # evaluates to True
        booltrue = ['yes','true','1','t']

        parname = parname.lower() # so we don't have to worry about case
        
        # Get paramlist index and check if parameter is valid
        try:
            ind = [par[0].lower() for par in self.paramlist].index(parname)
        except ValueError:
            msg = 'GetParam: There is no parameter named %s' % parname
            self.log.error(msg)
            raise KeyError(msg)
        # get from arguments if possible
        if self.arglist.has_key(parname):
            if self.arglist[parname] != None:
                self.log.debug('GetParam: Done (%s)' % parname )
                return self.arglist[parname]
        # get from config if possible
        tmp = {}
        for k in self.config[self.name].keys():
            tmp[k.lower()] = self.config[self.name][k]
        if tmp.has_key(parname):
            value = tmp[parname]
            out = self.paramlist[ind][1]
            if isinstance(out,(tuple,list)): # if it's a sequence
                ret = []
                for i in xrange(len(self.paramlist[ind][1])):
                    par = self.paramlist[ind][1][i]
                    if isinstance(par,bool):
                        if value[i].lower() in booltrue:
                            ret.append(True)
                        else: # default to False
                            ret.append(False)
                    else:
                        ret.append(type(par)(value[i]))
                # convert to tuple
                self.log.debug('GetParam: Done (%s)' % parname)
                return type(out)(ret)
            else:
                self.log.debug('GetParam: Done (%s)' % parname)
                if isinstance(out,bool):
                    if value.lower() in booltrue:
                        return True
                    else:
                        return False
                else:
                    return type(out)(value)
        # get default from parameter list
        ret = self.paramlist[ind][1]
        # return parameter
        self.log.debug('GetParam: done (%s)' % parname)
        return ret

    def getparam(self, parname):
        """ DEPRECATED - use getarg instead
            Returns the value of the parameter parname. The parameter is
            first searched for in self.arglist['parname'], then in
            config['stepname']['parname']. If the parameter is not found,
            a warning is returned and a KeyError is raised.
        """
        default = 0
        # get from params if possible
        if self.arglist.has_key(parname):
            self.log.debug('GetParam: Done (%s)' % parname )
            return self.arglist[parname]
        # get from config
        try:
            ret = self.config[self.name][parname]
        except KeyError:
            msg = ('GetParam: Missing %s entry in config: using default = %s'
                   % (parname,str(default)) )
            self.log.error('GetParam: Missing %s entry in config: using default = %s'
                          % (parname,str(default)) )
            raise KeyError(msg)
        # return parameter
        self.log.debug('GetParam: done (%s)' % parname)
        return ret
    
    def undo(self):
        """ Reverts the pipeline step to its status before the last call
        """
        self.log.debug('Undo: done')
    
    def reset(self):
        """ Resets the step to the same condition as it was when it was
            created. Internal variables are reset, any stored data is
            erased.
        """
        self.log.debug('Reset: done')
        
    def execute(self):
        """ Runs the pipe step as called from the command line:
            The first argument is used as input file name. Other
            special arguments are:
            - config = name of the configuration file object
            - test = runs the test function using the input file
            - loglevel = name of logging level (INFO is default)
            Other arguments are used as parameters to the pipe step.
        """
        ### Read Arguments
        # Set up argument parser - Generic parameters
        parser = argparse.ArgumentParser(description="Pipeline Step")
        parser.add_argument('inputfiles', type=str, default='', nargs='*',
                            help='input files pathname',)
        parser.add_argument('-t','--test', action='store_true',
                            help='runs the selftest of this pipe step')
        parser.add_argument('--loglevel', default='INFO', type=str,
                            choices=['DEBUG','INFO','WARN',
                                     'ERROR','CRITICAL'],
                            help='log level (default = INFO)')
        parser.add_argument('--logfile', default='', type=str,
                            help='logging file (default = none)')
        parser.add_argument('--config', default='', type=str,
                            help='pipeline configuration file (default = none)')
        # Add step-specific parameters from parlist
        for param in self.paramlist:
            # Comment: default = None because getarg gets default value from
            #          paramlist
            parser.add_argument('--'+param[0], type=type(param[1]),
                                default=None, help=param[2])
        # Get arguments - store dict in arglist
        args = parser.parse_args()
        self.arglist=vars(args)
        print(self.arglist)
        ### Process generic arguments
        # Set logging (add file handler if logfile != '')
        level = getattr(logging, args.loglevel.upper(), None)
        logging.basicConfig(level=level)
        if len(args.logfile) > 0:
            hand = logging.FileHandler(args.logfile)
            fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            hand.setFormatter(logging.Formatter(fmt))
            logging.getLogger().addHandler(hand)
        # Set configuration (load if specified)
        if len(args.config) > 0:
            datain = PipeData(config = args.config)
            self.config = datain.config
        elif not args.test: # Set config unless test is requested
            self.config = ConfigObj()
            self.config[self.name]={}
        # Check for test
        if args.test:
            self.test()
            return
        ### Reduce data
        if len(args.inputfiles) > 0:
            self.execfiles(args.inputfiles)
        else:
            # Warning - no input file
            self.log.warn('Execute: Missing input File')
        self.log.info('Execute: done')
        
    def execfiles(self, inputfiles):
        """ Runs several files from execute.
            This function is overwritten in MISO and MIMO steps
        """
        for filename in inputfiles:
            # Read input file
            self.datain = PipeData(config = self.config)
            self.datain.load(filename)
            # Call start - run and call end
            self.runstart(self.datain,self.arglist)
            self.run()
            self.runend(self.dataout)
            # Write output file
            self.dataout.save()

    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step parent')
        # test get and set parameters
        print("sampar=%.3f" % self.getarg('sampar'))
        # test function call
        #testout=self(1) # should raise TypeError
        if self.config != None:
            testin = PipeData(config=self.config)
        else:
            testin = PipeData(config=self.testconf)
        testin.filename = 'this.file.type.fts'
        testout=self(testin, sampar=5.0)
        print(testout.header)
        print(testout.filename)
        # log message
        self.log.info('Testing pipe step parent - Done')
    
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
    StepParent().execute()

""" === History ===
    2014-3-14  Marc Berthoud: Renamed callstart/callend -> runstart/runend
               - Took out code for MI and MO steps (not separate objects)
    2013-12-18 Marc Berthoud: 
               - added default test configuration file specification
               - removed config parameter for test() functions
    2012-12-23 Marc Berthoud: Added command line capability for pipe steps with
                              execute(), added paramlist.
    2009-10-27 Marc Berthoud, Ver0.1.1: Split call into callstart - run - callend
    2008-11-5 Marc Berthoud, Ver0.1: Wrote and tested (with HAWC text)
    2010-10-21 Marc Berthoud, Now configuration is loaded from datain
"""
