#!/usr/bin/env python
""" PIPE STEP ASTROMETRY- Version 1.0.0

    This pipe step calls the external program astrometry.net to add
    WCS information to the data.
    
    @author: Josh / Prechelt / Berthoud
"""
# TODO more descriptive for scale_lower and scale_upper

import logging # logging object library
import tempfile # temporary file library
import os # library for operating system calls
from astropy import wcs # to get WCS coordinates
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u
from darepype.drp import DataFits
from darepype.drp import StepParent
from astroquery.astrometry_net import AstrometryNet

class StepWebAstrometry(StepParent):
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
        self.name='webastrometry'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'WCS'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['timeout', 300,
                               'Timeout for running astrometry (seconds)'])
        self.paramlist.append(['radius', 5.,
                               'Search within this many degrees of the center RA and Dec'])
        self.paramlist.append(['scale_lower', 0.5,
                               'lower scale'])
        self.paramlist.append(['scale_upper', 2.,
                               'upper scale'])
        self.paramlist.append(['scale_units', 'arcsecperpix',
                               'Units of degrees/pixel'])
        self.paramlist.append(['api_key', 'XXXXXXXX',
                               'API key used for interfacing with Astrometry.net'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def astrometrymaster(self):
        '''
        Runs Astrometry.net on the inputted image
        '''
        
        ast = AstrometryNet()
        ast.api_key = self.getarg('api_key')

        self.log.debug("Now running ast.solve, get comfy this'll take a while")

        try:
            ra = Angle(self.datain.getheadval('RA'), unit=u.hour).degree
            dec = Angle(self.datain.getheadval('DEC'), unit=u.deg).degree
        except:
            self.wcs_out = ast.solve_from_image(self.datain.filename, force_image_upload = True, solve_timeout = self.getarg('timeout'), 
                                                scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), 
                                                scale_units = self.getarg('scale_units'))
        else:
            self.wcs_out = ast.solve_from_image(self.datain.filename, force_image_upload = True, ra_key = 'RA', dec_key = 'DEC', 
                                                solve_timeout = self.getarg('timeout'), radius = self.getarg('radius'), 
                                                scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), 
                                                scale_units = self.getarg('scale_units'))

        try:
            #Check if WCS exists
            self.wcs_out
        except:
            self.log.error("Unable to perform web astrometry.")
        else:
            self.log.debug('Web astrometry successful')


    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        
        self.dataout = DataFits(config=self.config)
        self.astrometrymaster()

        ## Post-processing
        self.dataout = self.datain.copy()
        self.dataout.header.update(self.wcs_out)
        # Add RA from astrometry
        w = wcs.WCS(self.dataout.header)
        n1 = float( self.dataout.header['NAXIS1']/2 )
        n2 = float( self.dataout.header['NAXIS2']/2 )
        ra, dec = w.all_pix2world(n1, n2, 1)
        # No update because update may affect accuracy of WCS solution
        # self.dataout.header['CRPIX1']=n1
        # self.dataout.header['CRPIX2']=n2
        # self.dataout.header['CRVAL1']=float(ra)
        # self.dataout.header['CRVAL2']=float(dec)
        self.dataout.header['RA'] = Angle(ra,  u.deg).to_string(unit=u.hour, sep=':')
        self.dataout.header['Dec'] = Angle(dec, u.deg).to_string(sep=':')
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
    StepWebAstrometry().execute()

""" === History ===
2020-08-10  -Removed overwriting of RAW files
2020-08-07  -Removed unused code from previous version, combined webastrometry
             and astrometrymaster, API key can now be specified in config
2020-07-27  -No longer redundantly updates certain header keywords
            -Got rid of masking, going to move to separate pipestep
2020-07-24  -Initial version of StepWebAstrometry
"""
