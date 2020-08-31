#!/usr/bin/env python
""" PIPE STEP MASK IMAGE- Version 1.1.0

    This pipe step masks the image before being sent into either 
    local or web astrometry as specified in config file. This code
    was adapted from code written by students in Spring 2020.
    
    @author: Josh
"""

import logging
from darepype.drp import DataFits 
from darepype.drp import StepParent
import numpy as np
import pandas as pd
import sep

class StepMaskImage(StepParent):
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
        self.name='maskimage'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'MSK'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['edgefrac', 0.4,
                               'edgefrac specifies the fraction of surrounding edge sources to preserve, so for example' +
                               'if edgefrac = 0.4 only sources falling within the center 2/10 of each axis are removed' +
                               'since sources falling within the outer 4/10 of each axis length is preserved along each edge'])
        self.paramlist.append(['rel_thresh', 1.5,
                               'Relative threshold pixel value for source extraction'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def mask(self):
        ''' Masks the input image file by extracting sources using SEP
            and filling in the fit ellipse of sources near the center
            with a null value
        '''
        # Background extraction and subtraction
        bkg = sep.Background(self.datain.image.astype('int32'))
        primary = self.datain.image.astype('int32') - bkg
        self.log.debug('Mask: subtracted background')

        # Source extraction 
        rel_thresh = self.getarg('rel_thresh')
        objects = sep.extract(primary, rel_thresh, err=bkg.globalrms)
        df = pd.DataFrame()
        df['x'] = objects['x']; df['y'] = objects['y']; df['a'] = objects['a']; df['b'] = objects['b']; df['theta'] = objects['theta']; df['npix'] = objects['npix']; df['FLUX'] = objects['cflux']
        self.log.debug('Mask: successfully extracted sources')

        # edgefrac specifies the fraction of sources near the edge to preserve
        edgefrac = self.getarg('edgefrac')
        df = df[ (df['x'] > float(primary.shape[0]) * edgefrac) & (df['x'] < float(primary.shape[0]) * (1 - edgefrac)) & (df['y'] > float(primary.shape[1]) * edgefrac) & (df['y'] < float(primary.shape[1]) * (1 - edgefrac))]
        
        x = np.linspace(0, primary.shape[0]-1, primary.shape[0])
        y = np.linspace(0, primary.shape[1]-1, primary.shape[1])
        X, Y = np.meshgrid(x, y)

        # For each source to be removed, see which pixels fall within the fit ellipse and fill
        # in a null value for the ones that do
        for index, row in df.iterrows():
            x0, y0, a, b, theta = row['x'], row['y'], row['a'], row['b'], row['theta']
            a *= 6
            b *= 6
            terma = (X - x0) * np.cos(theta) + (Y - y0) * np.sin(theta) 
            termb = (X - x0) * np.sin(theta) - (Y - y0) * np.cos(theta) 
            mask = (terma**2 / a**2) + (termb**2 / b**2)
            primary [mask < 1.] = np.ma.masked

        self.log.debug('Mask: done')
        
        return primary

    def run(self):
        self.dataout = DataFits(config=self.config)
        self.dataout.header = self.datain.header
        self.dataout.image = self.mask()
        self.dataout.filename = self.datain.filename.replace('.fits', '_MSK.fits')
        
if __name__ == "__main__":
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of 
    """
    StepMaskImage().execute()

""" === History ===
    2020-08-27  -Removed dependency on VOApy.py, cleaned up, added documentation
    2020-07-31  -Initial version of StepMaskImage
"""