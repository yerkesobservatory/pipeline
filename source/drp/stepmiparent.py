''' PIPE STEP Multi-Input PARENT

    This module defines a pipeline step with multiple input data objects
    and one output data object. All such pipe steps are descendants of
    this one.

    @author: berthoud
'''

import logging # logging object library
from drp.pipedata import PipeData # pipeline data object
from drp.stepparent import StepParent # pipe step parent object

class StepMIParent(StepParent):
    """ Pipeline multiple data input parent object
    """
    
    stepver = '1.0' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
            calls the setup function.
        """
        # call superclass constructor (calls setup)
        super(StepMIParent,self).__init__()
        # Change datain
        self.datain = [PipeData()]
        # set iomode
        self.iomode = 'MISO'
    
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
        self.name='parentmi'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
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
        # Return the first datain element
        self.dataout = self.datain[0]
        
    def runstart(self, data, arglist):
        """ Method to call at the beginning of the pipe step call.
            Mostly refers to runstart of stepparent
        """
        # Check input data - should be a list/tuple with PipeData objects
        if isinstance(data, (list,tuple)):
            for d in data:
                if not isinstance(d, PipeData):
                    msg = 'Invalid input data type: PipeData object is required'
                    self.log.error(msg)
                    raise TypeError('Runstart: '+msg)
        else:
            msg = 'Invalid input data type: List object is required'
            self.log.error(msg)
            raise TypeError('Runstart: '+msg)
        # Call parent runstart
        super(StepMIParent,self).runstart(data[0],arglist)
        
    def execfiles(self, inputfiles):
        """ Runs several files from execute.
        """
        if len(inputfiles) > 0:
            # Read input files to datain
            self.datain = []
            for filename in inputfiles:
                # Read input file
                data = PipeData(config = self.config)
                data.load(filename)
                self.datain.append(data)
            # Call start - run and call end
            self.runstart(self.datain,self.arglist)
            self.run()
            self.runend(self.dataout)
            # Write output file
            self.dataout.save()
            self.log.info('Execute: Saved result %s' % self.dataout.filename)
        else:
            # Warning - no input file(s)
            self.log.warn('Execute: Missing input File(s)')

    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)
        # read configuration
        if self.config != None:
            datain = PipeData(config=self.config)
        else:
            datain = PipeData(config=self.testconf)
        # generate 2 files
        datain.filename = 'this.file.type.fts'
        datain = [datain,datain]
        # run function call
        dataout = self(datain)      
        # test output
        print(type(dataout))
        print(dataout.header)
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
    StepMIParent().execute()

""" === History ===
    2014-3-14 Marc Berthoud: Written and tested version 0.1
"""
