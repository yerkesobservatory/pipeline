#!/usr/bin/env python
""" PIPE STEP SORTOBS- Version 1.0.0

    This pipe step sorts images taken by itzamna into 
    discrete folders based on observation sessions.
    
    @author: Josh

    NOTE: This step changes the filepath so that the files produced
    during data reduction will be output in the correct folder. It
    does not move/copy any RAW files and should be run at the 
    beginning of the pipeline.

"""

import re
import os
import logging
from darepype.drp import DataFits
from darepype.drp import StepParent

class StepSortObs(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.2' # pipe step version
    
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
        self.name='sortobs'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'RAW'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append Parameters
        self.paramlist.append(['pattern', '(^.+_([gri]-band|oiii|sii|clear|h-alpha))',
                               'Regex pattern used to get name by matching name_filter'])
        # Confirm end of setup
        self.log.debug('Setup: done')

    def run(self):
        # dataout will be identical to datain except for filepath
        self.dataout = self.datain.copy()
        
        # Regex matches object_filter_......fits  to get the name of observed object
        pattern = self.getarg('pattern')
        try: 
            orig_path = os.path.split(self.datain.filename)
            mtch = re.findall(pattern, orig_path[1])
            # re.findall() returns a list of length 1 containing a tuple of length n = # of matches,
            # so mtch[0][0] will be object_filter and mtch[0][1] will be filter based on the
            # grouping specified in pattern
            obj_name = mtch[0][0].replace(('_' + mtch[0][1]), '').strip()
        except:
            pass
        else:
            # If a folder matching the object name exists, place the file there. If not, 
            # make a folder first and then place the file there.
            new_path = os.path.join(os.path.split(orig_path[0])[0], obj_name)
            if os.path.exists(new_path):
                self.dataout.filename = os.path.join(new_path, orig_path[1])
            else:
                os.mkdir(new_path)
                self.dataout.filename = os.path.join(new_path, orig_path[1])

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
    StepSortObs().execute()

""" === History ===
2021-01-02  -Initial version
"""