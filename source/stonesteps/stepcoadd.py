""" PIPE STEP COADD - Version 1.0.0

    This module combines images and creates mosaics from input images.
    Should be run with like data inputs (same binning, etc.) with high quality WCS headers (like from astrometry.net)

    @author: Matt Merz
"""

import os # os library
import sys # sys library
import numpy # numpy library
import logging # logging object library
from astropy.wcs import WCS
from drizzle import drizzle as drz
from grizli.utils import make_wcsheader #function for making north-up WCS header and object
from drp.stepmiparent import StepMIParent
from drp.datafits import DataFits

class StepCoadd(StepMIParent):
    """ Stone Edge Pipeline Step Master Bias Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '1.0' # pipe step version
    
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
        self.name='coadd'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'coadd'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['kernel','square',
                               'Specifies the kernel used to determine spreading of input pixels onto output pixels \
                               - options are square, point, gaussian, smoothing, tophat - default is square'])
        self.paramlist.append(['pixfrac', 1.,
                               'The fraction of an output pixel(s) that an input pixels flux is confined to - default is 1.'])
        self.paramlist.append(['resolution', 1.,
                               'Plate-scale divisor for output image (higher gives more resolution, lower gives less) \
                               - default is 1.'])
        self.paramlist.append(['pad', 0,
                               'Extra padding outside maximum extent of inputs - default is 0'])
        self.paramlist.append(['fillval', '0.',
                               'Value for filling in the area(s) in the output where there is no input data - default is 0.'])
        self.paramlist.append(['drizzleweights','uniform',
                               'How each input image should be weighted when added to the output \
                               - options are exptime, expsq and uniform - default is uniform'])

    def run(self):
        """ Runs the mosaicing algorithm. The self.datain is run
        through the code, the result is in self.dataout.
        """
        #get all parameters
        if self.getarg('kernel') == 'smoothing':
            kernel = 'lanczos3'
        else:
            kernel = self.getarg('kernel')
        pixfrac = self.getarg('pixfrac')
        resolution = self.getarge('resolution')
        pad = self.getarg('pad')
        fillval = self.getarg('fillval')
        if self.getarg('drizzleweights') == 'uniform':
            driz_wt = ''
        else:
            driz_wt = self.getarg('drizzleweights')
        
        #set up lists of data and wcs from self.datain
        im_list = []
        wcs_list = []
        for file in self.datain:
            im_list.append(file.data)
            wcs_list.append(WCS(file.header))
            
        #calculate platescale of first input image
        det = np.linalg.det(wcs_list[0].wcs.cd)
        pscale = np.sqrt(np.abs(det))*3600.
        
        #creates north-up wcs which spans the needed area to contain all input images
        px = []
        py = []
        for wcs in wcs_list:
            px.extend(wcs.calc_footprint()[:,0])
            py.extend(wcs.calc_footprint()[:,1])
        x0 = (max(px)+min(px))/2.
        y0 = (max(py)+min(py))/2.
        sx = (max(px)-min(px))*np.cos(y0/180*np.pi)*3600 # arcsec
        sy = (max(py)-min(py))*3600 # arcsec
        fullwcs = grizli.utils.make_wcsheader(ra=x0, dec=y0, size=(sx+pad*2, sy+pad*2), pixscale=pixel_scale/resolution, get_hdu=False, theta=0)[1]
        
        #create drizzle object and add input images
        
        driz = drz.Drizzle(outwcs = fullwcs, pixfrac=pixfrac, kernel=kernel, fillval = fillval, wt_scl = driz_wt)
        for i in range(len(im_list)):
            driz.add_image(im_list[i], wcs_list[i])
        
        #creates output fits file from drizzle output
        self.dataout = DataFits(config = self.config)
        self.dataout.header=self.datain[0].header
        self.dataout.imageset(self.bias)
        self.dataout.filename = self.datain[0].filename
        
        # Add history
        self.dataout.setheadval('HISTORY','Coadd: %d files combined with %s kernel, pixfrac %f at %f times resolution' % \
            (len(filelist), kernel, pixfrac, resolution))

        
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
    StepCoadd().execute()  
    
""" === History ===
    2019-05-2 New step created for mosaicing and combining images - Matt Merz
"""