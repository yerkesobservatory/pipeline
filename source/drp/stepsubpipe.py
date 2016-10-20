#!/usr/bin/env python
""" PIPE STEP SUB PIPE - Version 1.0.0

    This module defines the HAWC pipeline step sub pipe object. This
    step runs a pipeline from command line as a separate process. It
    only works for pipes with SISO (single input single output) steps.
    This step is used to run pipeline steps which are resource hungry
    or otherwise unstable.
    
    @author: berthoud
"""

import logging # logging object library
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
        self.name='parent'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk' # taken from pipemode of subpipe
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['pipemode','unk','ProcName for returned files'])
        self.paramlist.append()
        self.paramlist.append(['sampar', 1.0,
            'Sample Parameter - parent only - no practical use'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        self.log.debug('Running step %s' % self.name)
    
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
    blah().execute()

""" === History ===
"""
