import sep
import pandas as pd
import VOApy as vp
from astropy import wcs
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table
from darepype.drp import DataFits 
from darepype.drp import StepParent
import logging

class StepSEP(StepParent):
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
        self.name='sep'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'SEP'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        # confirm end of setup
        self.log.debug('Setup: done')

    def source_extract(self):
        bkg = sep.Background(self.datain.image.astype('int32'))
        primary = self.datain.image.astype('int32') - bkg
        objects = sep.extract(primary, 1.5, err=bkg.globalrms)
        df = pd.DataFrame()
        df['x'] = objects['x']; df['y'] = objects['y']; df['a'] = objects['a']; df['b'] = objects['b']; df['theta'] = objects['theta']; df['npix'] = objects['npix']; df['FLUX'] = objects['cflux']
        table_image = Table(df.values)
        print(repr(df))
        table_fp = self.datain.filename.replace('.fits', '_TABLE.fits')
        table_image.write(table_fp, format='fits')
        return table_fp

    def run(self):
        self.dataout = DataFits(config=self.config)
        self.dataout.load(self.source_extract())
        self.dataout.header['RA'] = self.datain.header['RA']
        self.dataout.header['Dec'] = self.datain.header['Dec']
        self.dataout.save()
        
if __name__ == "__main__":
    StepSEP().execute()

