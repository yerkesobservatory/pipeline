#!/usr/bin/env python
""" PIPE STEP SUB PIPE - Version 1.0.0

    This module defines the HAWC pipeline step sub pipe object. This
    step runs a pipeline from command line as a separate process. It
    only works for pipes with SISO (single input single output) steps.
    This step is used to run pipeline steps which are resource hungry
    or otherwise unstable.
    
    @author: berthoud
"""

import os # os library
import subprocess # sys library
import logging # logging object library
import configobj # config object library
import drp.pipeline # for getting folder
from drp.pipedata import PipeData
from drp.stepparent import StepParent

class StepSubPipe(StepParent):
    """ HAWC Pipeline Step Sub Pipe Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version
    
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
        self.name='subpipe'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk' # taken from last pipestep of subpipe
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['pipeconf','',
                               "Pipeconf file to run subpipe. Default is '' i.e. use self.config"])
        self.paramlist.append(['pipemode','subtest','pipeline mode to run subpipe'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Get parameters
        pipemode = self.getarg('pipemode')
        pipeconf = self.getarg('pipeconf')
        # Save pipeline configuration file if no file given
        savedconf = False
        if len(pipeconf) < 1:
            pipeconf = os.path.join(os.path.split(self.datain.filename)[0],'subconf.txt')
            # Make a copy to not overwrite config
            conftemp = configobj.ConfigObj(self.config)
            conftemp.filename = pipeconf
            conftemp.write()
            savedconf = True
        # Load pipeconf into separate pipeobject
        confdata = PipeData(config=pipeconf)
        ### Save file if not saved
        if not os.path.exists(self.datain.filename):
            self.datain.save()
        ### Setup and run pipeline
        # Get pipeline command
        cmd = 'python ' + drp.pipeline.__file__
        # Add arguments
        logfile = os.path.join(os.path.split(self.datain.filename)[0],'sublog.txt')
        cmd += ' --loglevel DEBUG --logfile ' + logfile
        cmd += ' --pipemode ' + pipemode
        # Add config and filename
        cmd += ' ' + pipeconf
        cmd += ' ' + self.datain.filename
        # Run the pipeline
        self.log.debug('running command = %s' % cmd)
        subprocess.call(cmd, shell=True)
        ### Load log messages and submit them (optional)
        ### Get step name of last step of subpipe
        try:
            stepslist = list(confdata.config['mode_'+pipemode]['stepslist'])
            laststepind = len(stepslist)-1
            while(stepslist[laststepind] != 'save'): laststepind -= 1 # search for last save
            stepname = stepslist[laststepind-1]
        except:
            msg = 'Could not find last step for mode_%s in config=%s' % (pipemode, pipeconf)
            self.log.error(msg)
            raise ValueError(msg)        
        # Load the step (use existing confdata object)
        laststep = confdata.getobject(stepname)
        ### Load file to memory
        # Make filename
        self.procname = laststep.procname
        lastfile = self.datain.filenamebegin + laststep.procname.upper() + self.datain.filenameend
        # Load file
        self.dataout = PipeData(config=self.config)
        self.dataout.load(lastfile)
        ### Remove pipeconf if written
        if savedconf:
            os.remove(pipeconf)
    
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
    StepSubPipe().execute()

""" === History ===
"""
