#!/usr/bin/env python
""" PIPE STEP ASTROMETRY- Version 1.0.0

    This pipe step calls the external program astrometry.net to add
    WCS information to the data.
    
    @author: Prechelt / Berthoud
"""
# TODO add additional documentation, clean up

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

import VOApy as vp

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
                               'radius'])
        self.paramlist.append(['scale_lower', 0.5,
                               'lower scale'])
        self.paramlist.append(['scale_upper', 2.,
                               'upper scale'])
        self.paramlist.append(['scale_units', 'arcsecperpix',
                               'scale units'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def astrometrymaster(self, mask = False, outname = None, outpath = None):
        '''
        Master function for calling upon different attempt styles for WCS matching.
        local decides whether local astrometry.net code is used or web API
        mask decides whether the image is masked first
        The other parameters are necessary on a case-by-case basis depending on which attempt style is used.
        '''
        if mask:
            #Mask the image
            masked_data = self.datain.copy()
            mask = vp.unit(data = masked_data.image, header = self.datain.header)
            mask.extract_bkg()
            mask.subtract_bkg()
            mask.set_primary('bkg_sub')
            mask.extract_sources()
            mask.build_sources_table()
            mask.filter_sources(edgefrac = 0.4)
            mask.mask_sources()
            masked_image = fits.PrimaryHDU(mask.primary, header = self.datain.header)
            #Create a new file to point the matching procedures to instead
            maskfp = tempfile.NamedTemporaryFile(suffix=".fits",dir=outpath)
            masked_image.writeto(maskfp.name, overwrite = True)
            name = maskfp.name
        else:
            name = self.datain.filename
        
        self.webastrometry(name)

        try:
            #Check if WCS exists
            self.wcs_out
        except:
            self.log.error("Unable to perform web astrometry.")
            return False
        else:
            if bool(self.wcs_out):
                #Check if it was successful
                self.log.debug('Web astrometry successful')
                return True
            else:
                self.log.error("Unable to perform web astrometry.")
                return False

    def webastrometry(self, inputfile):
        '''
        Run astrometry alignment through the web. Uploads the image and returns the solve result.
        '''
        ast = AstrometryNet()
        ast.api_key = 'xpolczmnfaxzkihm'
        try:
            ra = Angle(self.datain.getheadval('RA'), unit=u.hour).degree
            dec = Angle(self.datain.getheadval('DEC'), unit=u.deg).degree
        except:
            self.wcs_out = ast.solve_from_image(inputfile, force_image_upload = True, solve_timeout = self.getarg('timeout'), scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), scale_units = self.getarg('scale_units'))
        else:
            self.wcs_out = ast.solve_from_image(inputfile, force_image_upload = True, ra_key = 'RA', dec_key = 'DEC', solve_timeout = self.getarg('timeout'), radius = self.getarg('radius'), scale_lower = self.getarg('scale_lower'), scale_upper = self.getarg('scale_upper'), scale_units = self.getarg('scale_units'))

    def run(self):
        """ Runs the data reduction algorithm. The self.datain is run
            through the code, the result is in self.dataout.
        """
        ### Preparation
        # construct a temp file name that astrometry will output
        fp = tempfile.NamedTemporaryFile(suffix=".fits",dir=os.getcwd())
        # split off path name, because a path that is too long causes remap to
        # crash sometimes
        outname = os.path.split(fp.name)[1]
        fp.close()
        # Add input file path to ouput file and make new name
        outpath = os.path.split(self.datain.filename)[0]
        # Make sure input data exists as file
        if not os.path.exists(self.datain.filename) :
            self.datain.save()
        origimg = self.datain.imageget()
        self.dataout = DataFits(config=self.config)
        #Create attempts. First bool is whether it's local. Second is whether the image is masked
        for attempt in [False, True]:
            #Loop through each created attempt until one is successful
            if attempt:
                success = self.astrometrymaster(mask = attempt, outpath = outpath)
            else:
                success = self.astrometrymaster(mask = attempt)
            if success:
                break
        assert success, "Unable to successfully match astrometry with any method."

        ## Post-processing
        self.dataout = self.datain.copy()
        self.dataout.header.update(self.wcs_out)
        self.dataout.filename = self.datain.filename
        self.dataout.save()
        # Add RA from astrometry
        w = wcs.WCS(self.dataout.header)
        n1 = float( self.dataout.header['NAXIS1']/2 )
        n2 = float( self.dataout.header['NAXIS2']/2 )
        ra, dec = w.all_pix2world(n1, n2, 1)
        self.dataout.header['CRPIX1']=n1
        self.dataout.header['CRPIX2']=n2
        self.dataout.header['CRVAL1']=float(ra)
        self.dataout.header['CRVAL2']=float(dec)
        self.dataout.header['RA'] = Angle(ra,  u.deg).to_string(unit=u.hour, sep=':')
        self.dataout.header['Dec']= Angle(dec, u.deg).to_string(sep=':')
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
2018-10-12 MGB: - Add code to try different --downsample factors
                - Add timeout for running astrometry.net
                - Renamed StepAstrometry from StepAstrometrica
2016-10-15 First version
"""
