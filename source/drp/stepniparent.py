""" PIPE STEP No-Input PARENT

    This module defines a pipeline step with no input file (i.e. self.datain is
    ignored).  It creates a single self.dataout pipedata object.  This module
    is a child of StepParent.
    
    (We are the Knights who say NI!)
    
    @author: chapman
"""

import logging # logging object library
from drp.pipedata import PipeData # pipeline data object
from drp.stepparent import StepParent # pipe step parent object

class StepNIParent(StepParent):
    def execfiles(self, inputfiles):
        """ Runs the step without an input file
        """
        
        self.datain = PipeData(config = self.config)
        # Call start - run and call end
        self.runstart(self.datain,self.arglist)
        self.run()
        self.runend(self.dataout)
        # Write output file
        self.dataout.save()
        self.log.info('Execute: Saved result %s' % self.dataout.filename)

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
    2014-12-30 Nicholas Chapman: Written and tested version 0.1
"""
