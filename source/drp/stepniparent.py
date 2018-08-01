""" PIPE STEP No-Input PARENT

    This module defines a pipeline step with no input file (i.e. self.datain is
    ignored).  It creates a single self.dataout pipedata object.  This module
    is a child of StepMOParent.
    
    (We are the Knights who say NI!)
    
    @author: chapman
"""

import logging # logging object library
from configobj import ConfigObj
from drp.dataparent import DataParent #pipeline data object
from drp.datafits import DataFits
from drp.stepmoparent import StepMOParent # pipe step parent object

class StepNIParent(StepMOParent):
    """ Pipeline no data input and multiple data output parent object
    """
    
    stepver = '1.0' # pipe step version
    
    def __init__(self):
        """ Constructor: Initialize data objects and variables
            calls the setup function.
        """
        # call superclass constructor (calls setup)
        super(StepNIParent,self).__init__()
        # Change dataout
        self.dataout = [DataParent()]
        # set iomode
        self.iomode = 'NIMO'

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
        self.name='parentni'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['filloutput', False,
            'If TRUE the step will return one pipedata object - else an empty list is returned (default)'])
        
    def __call__(self,datain=None, **arglist):
        """ Object Call: returns reduced input data
            No datain is allowed, in this case the step expects to receive
            the configuration in another way Ex: step.config=config
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
        if issubclass(data.__class__, DataParent):
            self.config = data.config
        elif self.config == None:
            msg='Invalid configuration information - aborting'
            self.log.error(msg)
            raise RuntimeError('Runstart: '+msg)
        # Check configuration
        if not isinstance(self.config,ConfigObj):
            msg='Invalid configuration information - aborting'
            self.log.error(msg)
            raise RuntimeError('Runstart: '+msg)
    
    def run(self):
        """ Runs the data reduction algorithm. The result is in self.dataout.
            NI steps a NIMO so multiple outputs; they have to return a list
            of data objects.
        """
        # Return dataout
        if self.getarg('filloutput'):
            # Make a datafits object and put it into output
            dout = DataFits(config = self.config)
            dout.filename = 'StepNIOuput.fits'
            self.dataout = [dout]
        else:
            # Return empty list
            self.dataout = []

    def execfiles(self, inputfiles):
        """ Runs the step without an input file
        """
        #self.datain = DataParent(config = self.config)
        # Call start - run and call end
        #self.runstart(self.datain,self.arglist)
        self.runstart(None, self.arglist)
        self.run()
        self.runend(self.dataout)
        # Write output file
        if len(self.dataout)>0:
            for data in self.dataout:
                data.save()
                self.log.info('Execute: Saved result %s' % data.filename)
        
    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a set of basic tests on the object
        """
        # log message
        self.log.info('Testing pipe step parent')
        # test function call
        #testout=self(1) # should raise TypeError
        if self.config != None:
            testin = DataParent(config=self.config)
        else:
            testin = DataParent(config=self.testconf)
        testin.filename = 'this.file.type.fts'
        testout=self(testin, sampar=5.0)
        print(testout.header)
        print(testout.filename)
        # test get and set parameters
        print("sampar=%.3f" % self.getarg('sampar'))
        # log message
        self.log.info('Testing pipe step parent - Done')

if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepniparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of 
    """
    StepNIParent().execute()

""" === History ===
    2018-07-12 Matt Merz: Updated and tested version 1.0
    2014-12-30 Nicholas Chapman: Written and tested version 0.1
"""
