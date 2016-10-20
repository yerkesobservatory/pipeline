""" PIPE LINE - Version 1.0.0
    This module defines the HAWC data reduction pipeline object. The
    object identifies, loads, reduces and stores raw HAWC data.
"""

import os # os library
import logging # logging object library
import traceback # error traceback library
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
    
    def __init__(self, filenames=None, config=None, instmode=None,
                 runmode=None):
        """ Constructor: Initialize the pipeline
        """
        # Admin Variables
        self.instmode = '' # instrument mode
        self.runmode = 'auto' # pipe running modes: AUTOmatic or INTERactive
        self.errstop = 0 # flag for pipe to stop if pipe step returns error
        self.config = None
        self.log = logging.getLogger('hawc.pipe.line')
        self.datalog = logging.getLogger('hawc.data.line')
        # Pipesteps Variables
        self.stepnames = [] # names of the pipeline steps (from step.name)
        self.steps = [] # list with the pipeline step objects
                        # ( steps[i] is of type stepname[i] )
                        # ( can't use dict b/c dicts have no order )
                        ### len(self.steps) > 0 indicates pipe is set up
        # Data Variables
        self.results = [] # list with results from each file for the most recent
                          #     step (Can be pipedata objects or lists thereof)
        self.resold = [] # list of previous complete result for each step NOT USED
        self.outputs = [] # list of most recent outputs for each step
        self.inputs = [] # list of input data, one element for each step
                         # these are only used for Multiple Input steps
        self.undone = 1 # flag indicating if undo function has been called NOT USED
        self.finals = [] # list of the results from last few final steps
        # Files Variables
        self.openfiles = [] # file name list of non-reduced files (fifo)
        self.reducedfiles = [] # file name list of reduced files
        # set up the pipeline
        self.setup(filenames=filenames, config=config, instmode=instmode,
                   runmode=runmode)
        # message
        self.log.debug('Init - Done')
    
    def reset(self):
        """ Resets the pipeline and prepares it to reduce a new set of
            images. Any reduced results will be deleted (make sure you save
            before you clear). All pipestep objects will be destroyed.
        """
        # Reset settings
        self.instmode=''
        # Reset pipeline
        for step in self.steps:
            del step
        self.stepnames=[]
        self.steps=[]
        self.openfiles=[]
        self.reducedfiles=[]
        # Reset Data
        self.results=[]
        self.resold=[]
        self.undone=1
        self.finals=[]
        # message
        self.log.debug('Cleared - Done')
    
    def setup(self, filenames=None, config=None, instmode=None, runmode=None,
              errstop=None):
        """ Using the information in the (first) file in 'filenames', this
            function sets up the pipeline. If instmode is None, the instrument
            mode in the file is used.
            Also sets configuration for the pipeline: if a filename is
            specified, the configuration file is read. The configuration
            object is returned.
            The configuration can always be changed, but once the instmode and
            the pipe steps have been set they can only be changed by using the
            reset() command, followed by a new setup().
            
            Setup returns the pipeline such that it can be used for initializing
        """
        ### set runmode
        if runmode ==  'auto' or runmode == 'inter':
            self.runmode=runmode
        elif runmode != None:
            self.log.warn('Setup: Invalid Runmode (%s)' % runmode)
        ### set configuration
        if config != None: self.config = PipeData().setconfig(config)
        ### set errstop from input or from config
        if errstop != None:
            self.errstop = errstop
        else:
            try: # necessary because various errors are possible
                self.errstop = int(self.config['general']['errorstop'] )
            except:
                pass
        ### get file header 
        if filenames is None:
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
            if namelen > 1: firstfile = filenames[0]
            else: firstfile = filenames
            # load header
            headdata = PipeData(config=self.config)
            headdata.loadhead(firstfile)
            retmsg='Loaded FITS file header'
        ### get instrument mode
        if instmode != None:
            if len(self.instmode):
                if self.instmode != instmode:
                    self.log.warn('Attempt to overwrite current INSTMODE ('
                                  +self.instmode+') with '+instmode
                                  +' - ignoring')
            else: self.instmode = instmode
        elif headdata != None and len(self.instmode) == 0 :
            try:
                self.instmode = headdata.getheadval('INSTMODE')
            except KeyError, error:
                self.log.error('Setup: Keyword INSTMODE not found in file %s'
                               % firstfile)
                raise error
        retmsg=retmsg+' / InstMode=%s' % self.instmode
        ### read steps (if config and instmode available and no steps loaded)
        if len(self.instmode)>0 and self.config!=None and len(self.steps)==0 :
            # get list of steps - copy to steplist
            conf_section='mode_'+self.instmode.lower()
            try:
                steplist = list(self.config[conf_section]['stepslist'])
            except KeyError, error:
                self.log.error('Setup: Missing steplist item in configuration '
                               +' for InstMode = ' + self.instmode )
                raise error
            self.stepnames = steplist
            # get list of step packages
            try:
                steppacks = self.config['general']['steppacks']
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
                            tmp = 'No module named %s' %stepname.lower()
                            if str(msg) == tmp: # module not present in directory
                                self.log.debug('Pipe step %s not found in %s' %
                                               (stepname, pack))
                            else: # module tries to import a missing package
                                raise
                        except:   # print out import errors not due to missing
                            raise # modules (e.g., typos in code)
                    # If not found -> Raise error
                    if stepmodule == None:
                        msg = 'Could not find step=%s' % stepname
                        msg += ' in instmode=%s' % self.instmode
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
                # Add results and resold
                self.results.append(PipeData(config=self.config))
                self.resold.append(PipeData(config=self.config))
            # add real stepnames from step.name, fill self.inputs
            for stepi in range(len(self.steps)):
                if not (self.stepnames[stepi] in 'load save'):
                    self.stepnames[stepi]=self.steps[stepi].name
                self.inputs.append([])
                self.outputs.append(None)
            retmsg=retmsg+' / Set up %s pipe steps' % steplist
        self.log.debug('Setup: done - %s' % retmsg)
        return self
    
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
    
    def save(self, filename, stepname='final'):
        """ Saves the result from the specified pipeline step under filename.
            If an error is raised if the step 'stepname' has not yet been
            run or if stepname is invalid. If stepname is 'final', the result
            from the last step is saved.
        """
        self.getresult(stepname=stepname).save(filename)
        self.log.debug('Save - Done')

    def __call__(self, filenames = None, config = None, instmode = None,
                 runmode = None):
        """ Runs the pipeline with the current files.
        """
        ### Setting up
        # Add files to file list
        if filenames != None:
            self.addfiles(filenames)
        # Configure
        self.setup(filenames=self.openfiles, config=config, instmode=instmode,
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
            # Check instmode - if it's wrong -> skip file
            try:
                fileinstmode = data.getheadval('INSTMODE')
                if self.instmode != fileinstmode:
                    self.log.warn('Call: File %s has wrong INSTMODE - skipping'
                                  % filename)
                    self.openfiles.remove(filename)
                    break
            except KeyError:
                self.log.warn('File %s does not have an INSTMODE keyword' %filename)

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
                        self.outputs[stepi] = data
                    elif self.stepnames[stepi] == 'save':
                        data.save()
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
                                # undo previous steps - NOT happening
                                # abort this file
                                errorflag = True
                                break # break out of loop over steps
                    # store step output in outputs[stepi]
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
            ### Here would be a good moment to stop reduction and enquire
            #   If more files arrived to reduce (in that case restart the process)
            ### run an MI step if there is one
            if stepi < len(self.steps) :
                self.log.debug('Call: Running multi-input step %d (%s)' % 
                               ( stepi, self.stepnames[stepi]) )
                # add current results to self.inputs[stepi]
                #     overwrite exitsting entries with identical filenames
                inpfilenames = [inp.filename for inp in self.inputs[stepi]]
                for res in self.results:
                    if res.filename in inpfilenames:
                        n = inpfilenames.index(res.filename)
                        self.infputs[stepi][n] = res
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
                                   + ' message=%s' % str(error) + ' - skipping file')
                    # print traceback
                    traceback.print_exc()
                    # raise error
                    raise(error)
                # store results in outputs[stepi] and results
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

    def call_Old(self, filenames=None, config=None, instmode=None, 
                 runmode=None):
        """ Reduces all files that have been added but not reduced.
            If filenames is specified it will add the files to the
            list of files to reduce. Returns the final reduced pipedata
            object.
        """
        ### Setting up
        # Add files to file list
        if filenames != None:
            self.addfiles(filenames)
        # Configure
        self.setup(filenames=self.openfiles, config=config, instmode=instmode,
                   runmode=runmode)
        # check if steps are available
        if len(self.steps) == 0 :
            self.log.error('Call: pipeline is not set up - aborting')
            raise RuntimeError('Call: pipeline is not set up - aborting')
        ### Loop through all files that have not yet been reduced
        while len(self.openfiles):
            # Load file header
            filename=self.openfiles[0]
            fileshort = os.path.split(filename)[-1]
            self.log.info('Start Reducing file %s' % fileshort)
            data=PipeData(config=self.config)
            data.loadhead(filename)
            # Check instmode - if it's wrong -> skip file
            if self.instmode != data.getheadval('INSTMODE') :
                self.log.warn('Call: File %s has wrong INSTMODE - skipping'
                              % filename)
                self.openfiles.remove(filename)
                break
            # Prepare reduction
            complete=1 # flag indicating last step was complete
            errorflag=0 # flag indicating an error in reducing the file
            # Store previous results
            for i in range(len(self.steps)):
                self.resold[i]=self.results[i]
            # Reduce file - go through steps
            for stepi in range(len(self.steps)) :
                # Catch special steps i.e. load / save
                if self.stepnames[stepi] == 'load':
                    data.load()
                    # save the result for this step
                    self.results[stepi]=data
                elif self.stepnames[stepi] == 'save':
                    data.save()
                    # save the result for this step
                    self.results[stepi]=data
                # if it is necessary to reduce the step
                # (i.e. previous steps all complete or running in interactive mode)
                elif complete or self.runmode == 'inter' :
                    ### reduce the data:
                    # if errors are not caught -> just run it and stop on errors
                    if self.errstop > 0:
                        data=self.steps[stepi](data)
                    # if errors are caught
                    else:
                        # run the step, see if it works
                        try:
                            data=self.steps[stepi](data)
                        except (ValueError, RuntimeError, TypeError,
                                KeyError, IndexError, IOError) , error: # in case of error
                            # list warning
                            self.log.warn('Call: Step %s'  % self.stepnames[stepi]
                                + ' for file %s' % filename + ' returned an error,'
                                + ' message=%s' % str(error) + ' - skipping file')
                            # print traceback
                            traceback.print_exc()
                            # undo previous steps
                            for j in range(stepi):
                                if not isinstance(self.steps[j],str):
                                    self.steps[j].undo()
                            # abort this file
                            errorflag = 1
                            break
                    # save the result for this step
                    self.results[stepi]=data
                    # check if the result is complete
                    #if not data.header['COMPLETE'] : complete = 0
            # -- end loop through steps
            # Remove the file from openfiles
            self.openfiles.remove(filename)
            # If reduction successful
            if not errorflag:
                self.undone = 0
                # add file to reduced list
                self.reducedfiles.append(filename)
                # Save Final result (if complete)
                if complete or self.runmode == 'inter':
                    self.finals.append(data)
                    # Remove oldest result if necessary
                    finalsaveN = int(self.config['pipeline']['finalsaveN'])
                    if len(self.finals) > finalsaveN:
                        del self.finals[0]
                # comment
                self.log.info('Reduced file %s' % fileshort)
            # If there was an error: restore pipeline to previous status
            else:
                self.undone = 1
                for i in range(len(self.steps)):
                    self.results[i]=self.resold[i]
        # -- end loop through files
        ### Return Result
        self.log.debug('Call: Finished Reducing Queued Data')
        return self.getresult()
      
    def undo(self):
        """ Overwrites current result for the specified data reduction
            step with the previous complete result.
        """
        # check if already undone
        if self.undone :
            self.log.warn('Undo: can not undo: last file already cleared')
            return
        # undo all steps
        for i in range(len(self.steps)):
            self.steps[i].undo()
            self.results[i]=self.resold[i]
        # pop from finals
        self.finals.pop()
        # remove last file from reducedfiles
        self.reducedfiles=self.reducedfiles[0:-1]
        # set undone=1
        self.undone=1
        self.log.debug('Undo: Done')
    
    def test(self):
        """ Test PipeLine Object
        """
        # log message
        self.log.info('Testing Pipe Line')
        # set setup
        config = StepParent().testconf
        self.setup(instmode='test',config=config,
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
                                'mode_chop', file) for file in infiles]
        self.addfiles(infiles)
        # reduce the files
        result = self()
        #self.undo()
        print('Shape of final Image = %s' % (result.image.shape,))
        # log message
        self.log.info('Testing Pipe Line - Done')
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    pipe=PipeLine()
    pipe.test()

""" === History ===
    2012-9-27 Marc Berthoud, Added dynamic importing of modules for pipe steps
    2010-10-27 Marc Berthoud, Ver 0.1.1: Shift to PipeData 0.1.1 and upgrades:
        - Adapted to new requirement that pipe steps don't have config
        - Changed fileslist and reducedN to openfiles and reducedfiles
        - Changed to allow step name in configuration to be object name
        - Clean up call, error and undo mechanisms
    2008-11-21 Marc Berthoud, Ver 0.1: Wrote and Tested
"""
