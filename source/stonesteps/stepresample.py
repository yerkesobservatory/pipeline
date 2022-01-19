#!/usr/bin/env python
""" PipeStep ReSample

    This pipe step downsamples the image by a certain factor. It also
    updates scale keywords specified in the configuration. The step updates
    wcs coordinates in the image. If the size of the image can't be divided 
    by samplefac, an exception is raised.
    
    @author: Marc Berthoud
    
"""

import logging # logging object library
import numpy as np # Numeric library
from darepype.drp.stepparent import StepParent

class StepReSample(StepParent):
    """ DarePype Step ReSample Object
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
        # Name of the pipeline reduction step (same case as config file)
        self.name='resample'
        # Identifier for saved file names.
        self.procname = 'RESMP'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['samplefac', 2, 'Downsample factor - integer'])
        #self.paramlist.append(['subframe', ':,:', 'subframe to keep in numpy ( i:j, k:l ) notation. Default is :,: i.e. full image'])
        #self.paramlist.append(['keysmulti','', '| separated list of fits keywords to multiply with downsample'])
        self.paramlist.append(['divkeys', [], 'list of header keywords to divide by samplefac'])
        self.paramlist.append(['multkeys', [], 'list of header keywords to multiply by samplefac'])

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        self.log.debug('Running step %s' % self.name)
        # Minimal code set dataout to datain
        self.dataout = self.datain
        # Downsample the image
        sfac = self.getarg('samplefac')
        imgsiz = self.datain.image.shape
        imgout = np.empty([imgsiz[0]//sfac,imgsiz[1]//sfac,4],dtype=self.datain.image.dtype)
        imgout[...,0] = self.datain.image[range(0,imgsiz[0],sfac),:][:,range(0,imgsiz[1],sfac)]
        imgout[...,1] = self.datain.image[range(1,imgsiz[0],sfac),:][:,range(0,imgsiz[1],sfac)]
        imgout[...,2] = self.datain.image[range(0,imgsiz[0],sfac),:][:,range(1,imgsiz[1],sfac)]
        imgout[...,3] = self.datain.image[range(1,imgsiz[0],sfac),:][:,range(1,imgsiz[1],sfac)]
        self.dataout.image = np.median(imgout,axis=2)
        # Divide keywords by sfac
        for key in self.getarg('divkeys'):
            if not key in self.datain.header: continue
            self.dataout.setheadval(key,self.datain.getheadval(key)/sfac)
        # Multiply keywords by sfac
        for key in self.getarg('multkeys'):
            if not key in self.datain.header: continue
            self.dataout.setheadval(key,self.datain.getheadval(key)*sfac)
            
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of parameters
    """
    StepReSample().execute()

""" === History ===
    2021-4: Wrote first version - Marc Berthoud
"""
