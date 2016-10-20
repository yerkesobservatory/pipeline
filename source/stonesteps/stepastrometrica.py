#!/usr/bin/env python
""" PIPE STEP ASTROMETRICA - Version 1.0.0

    This pipe step calls the external program astrometrica to add
    WCS information to the data.
    @author: Prechelt / Berthoud
"""

import logging # logging object library
import tempfile # temporary file library
import os # library for operating system calls
import subprocess # library to run subprocesses
from drp.pipedata import PipeData
from drp.stepparent import StepParent

class StepAstrometrica(StepParent):
    """ HAWC Pipeline Step Parent Object
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
        self.name='astrometrica'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'WCS'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['astrocmd', 'cp %s %s',
                               'Command to call astrometrica, should contain 2' +
                               'string placeholders for intput and output file'])
        self.paramlist.append(['verbose',False,
            'log full astrometrica output at DEBUG level'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Preparation
        # construct a temp file name that astrometrica will output
        fp = tempfile.NamedTemporaryFile(suffix=".fits",dir=os.getcwd())
        # split off path name, because a path that is too long causes remap to
        # crash sometimes
        outname = os.path.split(fp.name)[1]
        fp.close()
        # Make sure input data exists as file
        if not os.path.exists(self.datain.filename) :
            self.datain.save(fpathname)
        # Make command string
        command = self.getarg('astrocmd') % (self.datain.filename, outname)

        ### Run Astrometrica
        self.log.debug('running command = %s' % command)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        output, error = process.communicate()
        if self.getarg('verbose'):
            self.log.debug(output)

        ### Post processing
        # Read output file
        self.dataout = PipeData(config=self.config)
        try:
            self.dataout.load(outname)
            os.remove(outname)
            self.dataout.filename = self.datain.filename
        except Exception, error:
            self.log.error("Unable to open Astrometrica output file = %s"
                           % outname)
            raise error
        self.log.debug('Run: Done')
    
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
    StepAstrometrica().execute()

""" === History ===
2016-10-15 First version
"""
