#!/usr/bin/env python
""" PIPE STEP MASTER BIAS - Version 1.1.0

    Code for StepMasterBias in pipeline: does the following
    
    !!!!!!!1 Add what the step needs as inputs, what it does and how, what the outputs are !!!!!!!!

    @author: Matt Merz, et al
"""

from drp.pipedata import PipeData # pipeline data object
from drp.stepmiparent import StepMIParent # pipe step parent object
from drp.datafits import DataFits # Data Fits object

class StepMasterBias(StepMIParent):
    """ Stone Edge Pipeline Step Master Bias Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version

    def __init__(self):
        """ Constructor: Initialize data objects and variables
        """
        # call superclass constructor (calls setup)
        super(StepRGB,self).__init__()
        # list of data
        self.datalist = [] # used in run() for every new input data file
        # set configuration
        self.log.debug('Init: done')
    
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
        self.name='masterbias'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'MBIAS'
        # Set Logger for this pipe step
        self.log = logging.getLogger('stoneedge.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters !!!! WHAT PARAMETERS ARE NEEDED ????? !!!!!
        self.paramlist.append(['minpercent', 0.05, 
                               'Specifies the percentile for the minimum scaling'])
        self.paramlist.append(['maxpercent', 0.999,
                               'Specifies the percentile for the maximum scaling'])
        self.paramlist.append(['reductionmethod','normal',
                               'Specifies how the data should be reduced options are bla, bla, bla'])

    def run(self):
        """ Runs the combining algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        # List all filesnames of input files
        for fin in self.datain:
            self.log.debug("Input filename = %s" % fin.filename)
        # Make a dummy dataout
        self.dataout = DataFits(config = self.config)
        # OR Get dataout as copy of first datain
        self.dataout = self.datain[0].copy()
        # Do a bunch of stuff to commbine self.datain datasets into self.dataout
        if self.getarg('reductionmethod') == 'normal':
            # use normal way
            continue
        else:
            # use abnormal way
            continue
        # rename filename
        self.dataout.filename = os.path.join(self.getarg('outputfolder'),
                                             os.path.split(self.dataout.filename)[1])
        
        
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
    StepMasterBias().execute()
        
        
""" === History ===
    2018-??-?? New step created based on StepRGB - Emily, Matt
"""
