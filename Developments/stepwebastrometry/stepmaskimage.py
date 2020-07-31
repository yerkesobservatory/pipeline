#!/usr/bin/env python
""" PIPE STEP MASK IMAGE- Version 1.0.0

    This pipe step masks the image before being sent into either 
    local or web astrometry as specified in config file.
    
    @author: Josh
"""

import logging # logging object library
from astropy.io import fits # for manipulating fits objects
from darepype.drp import DataFits 
from darepype.drp import StepParent
import VOApy as vp # library for performing mask using SEP

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
        # confirm end of setup
        self.log.debug('Setup: done')

    def mask(self):
        '''
        Masks the input image file
        '''
        #Mask the image
        image_data = self.datain.image.astype('int32').byteswap(inplace=True).newbyteorder()
        mask = vp.unit(data = image_data, header = self.datain.header)
        mask.extract_bkg()
        mask.subtract_bkg()
        mask.set_primary('bkg_sub')
        mask.extract_sources()
        mask.build_sources_table()
        mask.filter_sources(edgefrac = 0.4)
        mask.mask_sources()
        masked_image = fits.PrimaryHDU(mask.primary, header = self.datain.header)
        mask_fp = self.datain.filename.replace('.fits', '_MSK.fits')
        masked_image.writeto(mask_fp)
        return mask_fp

    def run(self):
        self.dataout = DataFits(config=self.config)
        self.dataout.load(self.mask())
        
if __name__ == "__main__":
    StepMaskImage().execute()