#!/usr/bin/env python
""" PIPE STEP GET RAW - Version 1.0.0

    This module defines a utility pipe step that reverts the file name to
    the name of the original raw file by replacing the process name
    with "RAW". It can be used to make pipelines which reduce the same data
    multiple times by changing the steplist for a pipeline mode according
    to the following example:
    
    steplist = load, step1A, step1B, save, stepgetraw, load, step2A, step2B, save
    
    In this example the raw file is reduced with step1A and step1B, then
    reloaded and reduced with step2A and step2B.
    
    @author: berthoud
"""

import logging # logging object library
from drp.stepparent import StepParent

class StepGetRaw(StepParent):
    """ HAWC Pipeline Step Get Raw Object
        This pipe step reverts the file name to the name of the raw file by
        replacing the process name in the filename with RAW. Followed by a
        load in the steplist it will make the pipeline reload the raw file.
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
        self.name='getraw'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'raw'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # This step has no parameters

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        # Copy datain to dataout
        self.dataout = self.datain
        # No filenamechange is needed in this step, as the filename will be
        # changed in stepparent:updateheader()
        self.log.debug('Running step %s' % self.name)
        
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
    StepGetRaw().execute()

""" === History ===
    2016-3-19 Marc Berthoud: First version
"""
