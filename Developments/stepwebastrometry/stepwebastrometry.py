#!/usr/bin/env python
""" PIPE STEP WEBASTROMETRY- Version 1.2.0

    This pipe step calls the external program astrometry.net to add
    WCS information to the data.
    
    @author: Josh / Prechelt / Berthoud
"""
# TODO more descriptive for scale_lower and scale_upper

import logging # logging object library
from astropy import wcs # to get WCS coordinates
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.table import Table
import astropy.units as u
from darepype.drp import DataFits
from darepype.drp import StepParent
from astroquery.astrometry_net import AstrometryNet
import pandas as pd
import numpy as np

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
                               'Image plate scale units'])
        self.paramlist.append(['api_key', 'XXXXXXXX',
                               'API key used for interfacing with Astrometry.net'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def astrometrymaster(self):
        ''' Runs Astrometry.net on the inputted image/sources table
        '''
        
        ast = AstrometryNet()
        ast.api_key = self.getarg('api_key')


        # Check whether the file contains image or table data, prefer table data
        upload = True
        for hdu in fits.open(self.datain.filename):
            try:
                ### This try statement works with sextractor catalog files (...sex_cat.fits)
                # FITS files use big endian, so it must be converted to little endian before it can be used by numpy/pandas
                tbl = pd.DataFrame(np.array(hdu.data).byteswap(inplace=True).newbyteorder())
                # tbl.columns = ['X_IMAGE', 'Y_IMAGE', 'a', 'b', 'c', 'd', 'FLUX']
                tbl.columns = ['NUMBER', 'FLUX_APER', 'FLUX_AUTO', 'FLUXERR_AUTO', 'X_IMAGE', 'Y_IMAGE']
                # Astrometry.net requires the source table be sorted in descending order of flux
                tbl = tbl.sort_values(by='FLUX_APER', axis=0, ascending=False)
            except:
                try:
                    ### This try statement works with source extracted files containing image data (...SEXT.fits)
                    tbl = pd.DataFrame(np.array(hdu.data).byteswap(inplace=True).newbyteorder())
                    tbl.columns = ['ID', 'X_IMAGE', 'Y_IMAGE', 'FLUX', 'FLUX_ERR']
                    tbl = tbl.sort_values(by='FLUX', axis=0, ascending=False)
                except:
                    pass
                else:
                    upload = False
            else:
                upload = False


        # Determine the width/height of the image in pixels from the binning
        try:
            # The image width and height are stored in NAXIS1 and NAXIS 2 respectively
            imagew = float(self.datain.getheadval('NAXIS1'))
            imageh = float(self.datain.getheadval('NAXIS2'))
        except:
            # Binning 2 is pretty typical in a lot of cases
            imagew = 1024.
            imageh = 1024.
            self.log.debug('NAXIS1/2 keyword missing from header, assuming image width and height of %d and %d' % (imagew, imageh))
        else:
            self.log.debug('Image width is %d, image height is %d' % (imagew, imageh))


        if upload:
            self.log.debug('File to reduce has not been source extracted, uploading image to Astrometry.net')
        else:
            self.log.debug('File to reduce has been source extracted, uploading source table to Astrometry.net')

        self.log.debug("Now running ast.solve, get comfy this'll take a while")


        # Check if the header contains RA and Dec
        try:
            ra = Angle(self.datain.getheadval('RA'), unit=u.hour).degree
            dec = Angle(self.datain.getheadval('DEC'), unit=u.deg).degree
        except:
            # If the header is missing RA and Dec keywords, attempt to solve without them
            if upload:
                self.log.debug('Solving from image without RA/Dec')
                self.wcs_out = ast.solve_from_image(self.datain.filename, force_image_upload = True, solve_timeout = self.getarg('timeout'), 
                                                    scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), 
                                                    scale_units = self.getarg('scale_units'))
            else:
                self.log.debug('Solving from source list without RA/Dec')
                self.wcs_out = ast.solve_from_source_list(x=tbl['X_IMAGE'], y=tbl['Y_IMAGE'], image_width=imagew, image_height=imageh, 
                                                          solve_timeout=self.getarg('timeout'))
        else:
            # If the header contains RA and Dec, use them to solve
            if upload:
                self.log.debug('Solving from image with RA/Dec')
                self.wcs_out = ast.solve_from_image(self.datain.filename, force_image_upload = True, ra_key = 'RA', dec_key = 'DEC', 
                                                    solve_timeout = self.getarg('timeout'), radius = self.getarg('radius'), 
                                                    scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), 
                                                    scale_units = self.getarg('scale_units'))
            else:
                self.log.debug('Solving from source list with RA/Dec')
                self.wcs_out = ast.solve_from_source_list(x=tbl['X_IMAGE'], y=tbl['Y_IMAGE'], image_width=imagew, image_height=imageh, 
                                                          solve_timeout=self.getarg('timeout'), radius=self.getarg('radius'), center_ra=ra, 
                                                          center_dec=dec)


        # Check if the solution completed, though it should error if it times out anyway
        try:
            self.wcs_out
        except:
            self.log.error('Failed to find WCS solution')
        else:
            self.log.debug('Web astrometry completed successfully')


    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        
        self.dataout = DataFits(config=self.config)
        self.dataout = self.datain.copy()
        self.astrometrymaster()

        # Update RA/Dec from astrometry
        self.dataout.header.update(self.wcs_out)
        try:
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
        except:
            self.log.error('Run: Could not update RA/Dec from Astrometry')
        else:
            self.log.debug('Run: Updated RA/Dec from Astrometry')

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
2020-08-17  -Compatibility for files source extracted using sextractor 
2020-08-13  -Added ability to upload source tables to Astrometry
2020-08-10  -Removed overwriting of RAW files
2020-08-07  -Removed unused code from previous version, combined webastrometry
             and astrometrymaster, API key can now be specified in config
2020-07-27  -No longer redundantly updates certain header keywords
            -Got rid of masking, going to move to separate pipestep
2020-07-24  -Initial version of StepWebAstrometry
"""
