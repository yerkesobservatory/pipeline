""" PIPE LINE - Version 1.0.0
    This module defines the HAWC data reduction pipeline object. The
    object identifies, loads, reduces and stores raw HAWC data.
"""

import os # os library
import logging # logging object library
import traceback # error traceback library
import argparse # Argument parsing library
import configobj # Config Object library
from drp.pipedata import PipeData # Pipeline Data object
from drp.stepparent import StepParent # Step Parent
from drp.stepmiparent import StepMIParent # Step Multiple Input Parent

class PipeLine(object):
    """ HAWC Pipeline Object
        The object receives filesnames or loaded pipedata objects,
        reduces them and adds them to the current data. The pipeline object
        identifies the instrument mode and loads the correct configuration
        file. The configuration file also contains information which pipe
        steps to use.
    """

    def __init__(self, filenames=None, config=None, pipemode=None,
                 runmode=None, savemem=False):
        """ Constructor: Initialize the pipeline
        """
        # Admin Variables
        self.pipemode = '' # pipe mode: determines which mode the pipeline
                            # is running, i.e. which steps are run
        self.runmode = 'auto' # pipe running modes: AUTOmatic or INTERactive
        self.errstop = 0 # flag for pipe to stop if pipe step returns error
        self.config = None # configuration object
        self.memFlag = savemem # set to True to conserve memory.  This will
                                # not save outputs of each step, and alter
                                # the behaviour of self.getresults()
        self.log = logging.getLogger('pipe.line')
        # Pipesteps Variables
        self.stepnames = [] # names of the pipeline steps (from step.name)
        self.steps = [] # list with the pipeline step objects
                        # ( steps[i] is of type stepname[i] )
                        # ( can't use dict b/c dicts have no order )
                        ### len(self.steps) > 0 indicates pipe is set up
        # Data Variables
        self.results = [] # list with results from each file for the most recent
                            #     step (Can be pipedata objects or lists thereof)
        self.outputs = [] # list of most recent outputs for each step
        self.inputs = [] # list of input data, one element for each step
                        # these are only used for Multiple Input steps
        self.finals = [] # list of the results from last few final steps
        # Files Variables
        self.openfiles = [] # file name list of non-reduced files (fifo)
        self.reducedfiles = [] # file name list of reduced files
        # set up the pipeline
        self.setup(filenames=filenames, config=config, pipemode=pipemode,
                   runmode=runmode)
        # message
        self.log.debug('Init - Done')

    def reset(self):
        """ Resets the pipeline and prepares it to reduce a new set of
            images. Any reduced results will be deleted (make sure you save
            before you clear). All pipestep objects will be destroyed.
        """
        # Reset settings
        self.pipemode = ''
        # Reset pipeline
        for step in self.steps:
            del step
        self.stepnames=[]
        self.steps=[]
        self.openfiles=[]
        self.reducedfiles=[]
        # Reset Data
        self.results=[]
        self.finals=[]
        # message
        self.log.debug('Cleared - Done')

    def setup(self, filenames=None, config=None, pipemode=None,
              runmode=None, errstop=None):
        """ Using the information in the (first) file in 'filenames', this
            function sets up the pipeline. If pipemode is None, the header
            keywords in the file are used to determine pipemode.
            Also sets configuration for the pipeline: if a filename is
            specified, the configuration file is read. The configuration
            object is returned.
            The configuration can always be changed, but once the pipemode and
            the pipe steps have been set they can only be changed by using the
            reset() command, followed by a new setup().
        """
        ### set runmode
        if runmode ==  'auto' or runmode == 'inter':
            self.runmode=runmode
        elif runmode != None:
            self.log.warn('Setup: Invalid Runmode (%s)' % runmode)
        ### set configuration
        if config != None:
            self.config = PipeData().setconfig(config)
        ### set errstop from input or from config
        if errstop != None:
            self.errstop = errstop
        else:
            try: # necessary because various errors are possible
                self.errstop = int(self.config['general']['errorstop'] )
            except:
                pass
        ### get file header
        if filenames == None or len(filenames) == 0:
            # no file -> say it
            headdata = None
            retmsg='no FITS files given'
        else:
            # try to get the length of filenames[0]
            try:
                namelen=len(filenames[0])
            except TypeError, error:
                self.log.error('Setup: invalid filenames type')
                raise error
            except IndexError, error:
                self.log.error('Setup: no input files')
                raise error
            # if len(filenames[0]) > 1 it should be a filename (from a list)
            if namelen > 1:
                firstfile = filenames[0]
            else:
                firstfile = filenames
            # load header
            headdata = PipeData(config=self.config)
            headdata.loadhead(firstfile)
            retmsg='Loaded FITS file header'
        ### get pipeline mode
        # If available: from input parameter
        if pipemode != None:
            if len(self.pipemode):
                if self.pipemode != pipemode:
                    self.log.warn('Attempt to overwrite current PIPEMODE ('
                                  +self.pipemode+') with '+pipemode
                                  +' - ignoring')
            else:
                self.pipemode = pipemode
        # else: get pipeline mode from file header
        elif headdata != None and len(self.pipemode) == 0 :
            self.pipemode = self.getpipemode(headdata)
            if len(self.pipemode) == 0:
                msg = "Setup: No Pipe Mode found matching initial data"
                self.log.error(msg)
                raise KeyError,msg
        retmsg=retmsg+' / PipeMode=%s' % self.pipemode
        ### read steps (if config and pipemode available and no steps loaded)
        if ( len(self.pipemode) > 0 and self.config != None and 
             len(self.steps) == 0 ):
            self.log.info('Pipe Mode = %s' %self.pipemode)
            conf_section = 'mode_' + self.pipemode
            # Update configuration based on mode preferences
            self.config = configobj.ConfigObj(self.config) # make copy to preserve orig
            modeconf = self.config[conf_section]
            for item in modeconf.sections:
                for subitem in modeconf[item]:
                    self.config[item][subitem]=modeconf[item][subitem]
                    self.log.debug('Setup: change conf[%s][%s] to %s' %
                                   (item,subitem,modeconf[item][subitem]))
            # get list of steps - copy to steplist
            try:
                steplist = list(self.config[conf_section]['stepslist'])
            except KeyError, error:
                self.log.error('Setup: Missing steplist item in configuration '
                               +' for PipeMode = ' + self.pipemode )
                raise error
            self.stepnames = steplist
            # get list of step packages
            try:
                steppacks = self.config['general']['steppacks']
                if isinstance(steppacks,str): # ensure we have a list if only
                    steppacks = [steppacks]    # 1 item, steppacks is str
            except KeyError, error:
                self.log.error('Setup: Missing steppacks item in configuration')
                raise error
            # initialize inputs and outputs
            self.inputs = []
            self.outputs = []
            # get the steps - loop through stepnames
            for stepname in self.stepnames:
                # Add step
                if stepname in 'load save':
                    self.steps.append(stepname)
                else:
                    # Loop through step packages
                    stepmodule = None
                    for pack in steppacks:
                        if pack == '.': # look in local directory
                            mod = stepname.lower()
                        else:
                            mod = '%s.%s' % (pack,stepname.lower())
                        self.log.debug('looking for %s' % mod)
                        # import the module
                        try:
                            stepmodule = __import__(mod, globals(), locals(),
                                                    [stepname])
                            self.log.debug('Pipe step %s found in %s' %
                                          (stepname,pack))
                            break
                        except ImportError,msg:
                            tmp = 'No module named %s' % stepname.lower()
                            if str(msg).startswith(tmp): # module not present in directory
                                self.log.debug('Pipe step %s not found in %s' %
                                               (stepname, pack))
                            else: # module tries to import a missing package
                                raise
                        except:   # print out import errors not due to missing
                            raise # modules (e.g., typos in code)
                    # If not found -> Raise error
                    if stepmodule == None:
                        msg = 'Could not find step=%s' % stepname
                        msg += ' in pipemode=%s' % self.pipemode
                        self.log.error('Setup: ' + msg)
                        raise ImportError(msg)
                    # Make a step instance and add to step list
                    try:
                        self.steps.append(stepmodule.__dict__[stepname]())
                    except KeyError, error:
                        msg = 'Pipe step=%s' % stepname
                        msg+= ' not found in module=%s' % mod
                        self.log.error('Setup: %s' % msg)
                        raise error
                # Add results
                self.results.append(PipeData(config=self.config))
            # add real stepnames from step.name, fill self.inputs
            for stepi in range(len(self.steps)):
                if not (self.stepnames[stepi] in 'load save'):
                    self.stepnames[stepi]=self.steps[stepi].name
                self.inputs.append([])
                self.outputs.append(None)
            retmsg=retmsg+' / Set up %s pipe steps' % steplist
        self.log.debug('Setup: done - %s' % retmsg)

    def getpipemode(self,data):
        """ Searches for an appropriate pipe mode in the config file, given
            the header values in the passed data. Tries to mach all key=value
            pairs in the datakeys value of the mode_ entries in the config
            file.  Returns name of first drp mode that matches
            the data.  Returns ' ' if no matching drp mode found.  
        """
        for sec in self.config.sections:
            if sec.startswith('mode_'):
                # Get the datakeys and make list of lists
                try:
                    datakeys = self.config[sec]['datakeys'].split('|')
                except KeyError:
                    self.log.warn("In configuration, missing datakeys for mode=%s" % sec)
                    continue
                datakeys = [dk.split('=') for dk in datakeys]
                #print(sec,datakeys)
                # Check all keywords in the file
                check = True
                for dk in datakeys:
                    try:
                        value = data.getheadval(dk[0].strip())
                        if value.upper().strip() != dk[1].upper().strip(): check = False
                        #print(value.upper().strip(),dk[1].upper().strip())
                    except KeyError:
                        self.log.debug('GetPipeMode: Data is missing key=%s', dk[0])
                        check = False
                if check:
                    self.log.debug('GetPipeMode: Found mode=%s' % sec[5:])
                    return sec[5:] # return mode name w/o 'mode_'
        # should only get here in no matching drp mode found
        msg = "GetPipeMode: No DRP mode found matching data = %s" % data.filename
        self.log.warn(msg)
        return ''

    def getsteps(self):
        """ Returns the list of the pipeline steps in this pipeline.
            If the list of steps is not yet known, an empty list is
            returned.
        """
        self.log.debug('GetSteps: done')
        return self.steps

    def getresult(self, stepname='final'):
        """ Returns the last result from the 'stepname' pipeline reduction
            step. Returns null if stepname is invalid or if that reduction
            step has not yet been run. If stepname is 'final' the result from
            the last step is returned. If stepname is 'finals' the list of
            all the final results is returned.
        """
        # get index
        if stepname.lower() == 'final':
            resind = len(self.steps)-1
        else:
            try:
                resind = self.stepnames.index(stepname.lower())
            except ValueError, error:
                self.log.error("GetResult: invalid step <%s> requested"
                               % stepname)
                raise error
        # check if data is available
        if len(self.reducedfiles) < 1:
            self.log.warning('GetResult: no results yet')
            return 0
        # get result
        if self.memFlag: # memory conserved, so returns self.results
            self.log.warning('GetResult: Memory conserved, so result of last step returned')
            ret = self.results[-1]
        else:
            ret = self.outputs[resind]
        if isinstance(ret, (list,tuple)):
            ret = ret[0]
        # return result
        self.log.debug('GetResult: done')
        return ret

    def addfiles(self, filenames):
        """ Prepares the listed file(s) for reduction. If no files
            have been loaded, the instrument mode is determined from the
            first file and the corresponding configuration file is loaded.
            The files have to be of the correct mode.
        """
        # try to get the length of filenames[0]
        try:
            namelen=len(filenames[0])
        except TypeError, error:
            self.log.error('Add: invalid filenames type')
            raise error
        # add the files to the list (check if each file exists)
        if namelen < 2: filenames=[filenames]
        for filename in filenames:
            if os.path.isfile(filename):
                self.openfiles.append(filename)
            else:
                self.log.warn('Add: could not find file (%s) - ignored'
                              % filename)
        self.log.debug('Add - Done')

    def save(self, filename=None, stepname='final'):
        """ Saves the result from the specified pipeline step under filename.
            If an error is raised if the step 'stepname' has not yet been
            run or if stepname is invalid. If stepname is 'final', the result
            from the last step is saved.
        """
        self.getresult(stepname=stepname).save(filename)
        self.log.debug('Save - Done')

    def __call__(self, filenames = None, config = None, pipemode = None,
                 runmode = None, force=False):
        """ Runs the pipeline with the current files.

            Files with instmode/instcfg that do not match the pipeline's
            instmode/instcfg, are skipped.  This behavior can be overwritten
            by setting force = True (i.e. pipeline will attempt non-matching
            files anway).
        """
        ### Setting up
        # Add files to file list
        if filenames != None:
            self.addfiles(filenames)
        # Configure
        self.setup(filenames=self.openfiles, config=config, pipemode=pipemode,
                   runmode=runmode)
        # check if steps are available
        if len(self.steps) == 0 :
            self.log.error('Call: pipeline is not set up - aborting')
            raise RuntimeError('Call: pipeline is not set up - aborting')
        ### Initialize results:
        #   - make a pipedata object for each file with filename and header
        #   - check if it's a valid file (instmode)
        # Reset resuls
        self.results = []
        # Loop through files in openfiles
        for filename in self.openfiles:
            # Make pipedata object
            data=PipeData(config=self.config)
            # Load file header
            data.loadhead(filename)
            fileshort = os.path.split(filename)[-1]
            self.log.info('Preparing file %s' % fileshort)
            # Check if instmode and instcfg match.  If not, skip file
            filepipemode = self.getpipemode(data)
            if self.pipemode != filepipemode:
                if force is True:
                    self.log.warn('Call: File %s does not fit PipeMode datakeys  - Will attempt to process' %filename)
                else:
                    self.log.warn('Call: File %s does not fit PipeMode datakeys - skipping' %filename)
                    self.openfiles.remove(filename)
                    continue
            # Add file to results
            self.results.append(data)
            # Add to reducedfiles
            self.reducedfiles.append(filename)
        # Remove the files from openfiles
        self.openfiles = []
        ### Loop over steps:
        if len(self.results) == 0:
            self.log.error('Call: Zero files to process - quitting')
            stepi = len(self.steps)
        else:
            stepi = 0
        while stepi < len(self.steps):
            # Get stepstart and stepend i.e. range of next SISO steps
            # Point stepend to first MI step after stepstart or last step
            stepstart = stepi
            nextMI = False # flag indicating next step is MI
            stepend = stepi
            while not nextMI and stepend < len(self.steps) :
                if isinstance(self.steps[stepend], StepMIParent):
                    nextMI = True
                else:
                    stepend += 1
            ### run all SISO steps for each file in results
            #   files with errors are removed from results
            filei = 0
            self.log.debug('Call: Running steps %d to %d for all %d files'
                           % (stepstart, stepend-1, len(self.results)))
            while filei < len(self.results):
                errorflag = False
                stepi = stepstart
                while stepi < stepend:
                    data = self.results[filei]
                    # catch special steps (load/save)
                    if self.stepnames[stepi] == 'load':
                        data.load()
                        if not self.memFlag:
                            self.outputs[stepi] = data
                    elif self.stepnames[stepi] == 'save':
                        data.save()
                        if not self.memFlag:
                            self.outputs[stepi] = data
                    # run step[stepi]
                    else:
                        self.log.debug('Call: Starting step %s on file %s' %
                                       ( self.stepnames[stepi], data.filename) )
                        if self.errstop > 0:
                            # if errors are not caught -> run the step
                            data = self.steps[stepi](data)
                        else:
                            # try the step catch and handle error
                            try:
                                data = self.steps[stepi](data)
                            except (ValueError, RuntimeError, TypeError,
                                    KeyError, IndexError, IOError), error:
                                # list warning
                                self.log.warn('Call: Step %s'  % self.stepnames[stepi]
                                    + ' for file %s' % data.filename + ' returned an error,'
                                    + ' message=%s' % str(error) + ' - skipping file')
                                # print traceback
                                traceback.print_exc()
                                # abort this file
                                errorflag = True
                                break # break out of loop over steps
                        self.steps[stepi].datain = None
                        self.steps[stepi].dataout = None
                    # store step output in outputs[stepi]
                    if not self.memFlag:
                        self.outputs[stepi] = data
                    self.results[filei] = data
                    stepi += 1
                # -- end of loop over steps between stepstart and stepend
                # If there's an error in file:
                if errorflag:
                    # Remove file from results forward stepi
                    self.results.remove(self.results[filei])
                    stepi = stepend
                else:
                    # set index for next file
                    filei += 1
            # -- end of loop over files
            ### Here would be a good moment to stop reduction and inquire
            #   If more files arrived to reduce (in that case restart the process)
            ### run an MI step if there is one
            if stepi < len(self.steps) :
                self.log.debug('Call: Running multi-input step %d (%s)' %
                               ( stepi, self.stepnames[stepi]) )
                # add current results to self.inputs[stepi]
                #     overwrite existing entries with identical filenames
                inpfilenames = [inp.filename for inp in self.inputs[stepi]]
                for res in self.results:
                    if res.filename in inpfilenames:
                        n = inpfilenames.index(res.filename)
                        self.inputs[stepi][n] = res
                    else:
                        self.inputs[stepi].append(res)
                # Run the next MI step with the inputs
                try:
                    data = self.steps[stepi](self.inputs[stepi])
                # if there's an error. stop and return an error (no matter what)
                #     make sure results = []
                except (ValueError, RuntimeError, TypeError,
                        KeyError, IndexError, IOError), error:
                    # list error
                    self.log.error('Call: Step %s'  % self.stepnames[stepi]
                                   + ' for file %s' % self.inputs[stepi][0].filename + ' returned an error,'
                                   + ' message=%s' %str(error) + ' - skipping file')
                    # print traceback
                    traceback.print_exc()
                    # raise error
                    raise(error)
                #self.steps[stepi].reset()
                self.steps[stepi].datain = None
                self.steps[stepi].dataout = None
                # store results in outputs[stepi] and results
                if not self.memFlag:
                    self.outputs[stepi] = data
                if isinstance(data, (list,tuple)):
                    self.results = data
                else:
                    self.results = [data]
                # increase step index
                stepi += 1
        # -- end loop over all steps
        # add output from last step to self.finals
        if len(self.results) > 0:
            self.finals.append(self.results[-1])
        ### Return Result
        self.log.debug('Call: Finished Reducing Queued Data')
        return self.getresult()

    def test(self):
        """ Test PipeLine Object
        """
        # log message
        self.log.info('Testing Pipe Line')
        # set setup
        config = StepParent().testconf
        self.setup(pipemode='chop',config=config,
                   runmode='auto', errstop = 1)
        #self.setup(filenames='080211_000_00Ha010.dog.raw.dogscan.fits',
        #self.setup(filenames='short.fft.raw.fits',
        #           config='pipe_test_conf.txt')
        # get results and test
        #print self.getsteps()
        #print self.getresult('demod')
        # add the files
        print('Reduction Steps = %s' % (self.stepnames,))
        #infiles=['mode_fft/short.fft.raw.fits']
        infiles = ['120306_000_00HA005.chop.raw.fits',
                   '120306_000_00HA006.chop.raw.fits']
        infiles = [os.path.join(self.config['testing']['testpath'],
                                'hawc1G/mode_chop', ifile) for ifile in infiles]
        self.addfiles(infiles)
        # reduce the files
        result = self()
        print('Shape of final Image = %s' % (result.image.shape,))
        # log message
        self.log.info('Testing Pipe Line - Done')


    def execute(self):
        """ Runs the pipeline as called from the command line:
            The first argument is used as configuration file name.
            Following arguments are used as input files. Other
            special arguments are:
            - test: runs the test function
            - loglevel: name of logging level (INFO is default)
            - logfile: log file name
        """
                ### Read Arguments
        # Set up argument parser - Generic parameters
        parser = argparse.ArgumentParser(description="Pipe Line")
        parser.add_argument('config', default = 'pipeconf.txt', type=str, nargs='?',
                            help='pipeline configuration file (default = pipeconf.txt)')
        parser.add_argument('inputfiles', type=str, default='', nargs='*',
                            help='input files pathname',)
        parser.add_argument('-t','--test', action='store_true',
                            help='runs the selftest of the pipeline')
        parser.add_argument('--loglevel', default='INFO', type=str,
                            choices=['DEBUG','INFO','WARN',
                                     'ERROR','CRITICAL'],
                            help='log level (default = INFO)')
        parser.add_argument('--logfile', default='', type=str,
                            help='logging file (default = none)')
        parser.add_argument('--pipemode', default='', type=str,
                            help='pipeline mode (default = none)')
        # Get arguments - store dict in arglist
        args = parser.parse_args()
        self.arglist=vars(args)
        ### Process generic arguments
        # Set logging (add file handler if logfile != '')
        level = getattr(logging, args.loglevel.upper(), None)
        logging.basicConfig(level=level)
        if len(args.logfile) > 0:
            fhand = logging.FileHandler(args.logfile)
            fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            fhand.setFormatter(logging.Formatter(fmt))
            logging.getLogger().addHandler(fhand)
        # Check for test
        if args.test:
            self.test()
            return
        # Get pipemode
        pipemode = None
        force = False
        if len(args.pipemode) > 0:
            pipemode = args.pipemode
            force = True
        # Set configuration
        datain = PipeData(config = args.config)
        self.config = datain.config
        ### Reduce data
        if len(args.inputfiles) > 0:
            self(args.inputfiles,pipemode=pipemode,force=force)
            self.save()
        else:
            # Warning - no input file
            self.log.warn('Execute: Missing input File')
        self.log.info('Execute: done')


if __name__ == '__main__':
    """ Main function to run the pipe line from command line on a file.
        Command:
          python stepparent.py pipeconf.txt input1.fits input2.fits -arg1 -arg2 . . .
        where pipeconf.txt stands for the (path)name of the pipeline
        configuration file and inputN.fits for the input file(s).

        Standard arguments:
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of commands.
    """
    PipeLine().execute()

""" === History ===
    2016-4-14 Marc Berthoud: Added ability to specify configuration changes in
               mode prefference: The mode_pipemode section of the pipeline can
               now contain subsections with changes to be applied to any other
               section in the pipeline configuration.
    2015-10-12 Marc Berthoud: PIPEMODE UPDATE: DRP now uses pipemode keyword with
               datakeys config file key to determine which pipeline mode should
               be run.
    2015-10-7 Marc Berthoud: Removed unused undo and resold 
    2014-7-31 Marc Berthoud, Added execute() to run pipeline from the
        command line.
    2012-9-27 Marc Berthoud, Added dynamic importing of modules for pipe steps
    2010-10-27 Marc Berthoud, Ver 0.1.1: Shift to PipeData 0.1.1 and upgrades:
        - Adapted to new requirement that pipe steps don't have config
        - Changed fileslist and reducedN to openfiles and reducedfiles
        - Changed to allow step name in configuration to be object name
        - Clean up call, error and undo mechanisms
    2008-11-21 Marc Berthoud, Ver 0.1: Wrote and Tested
"""
