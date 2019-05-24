#!/usr/bin/env python
""" PIPE STEP COADD - Version 1.0.0

    This module combines images and creates mosaics from input images.
    Should be run with like data inputs (same binning, etc.) with high quality WCS (e.x. from astrometry.net)

    @author: Matt Merz
"""
import os # os library
import sys # sys library
import numpy as np# numpy library
import logging # logging object library
from astropy.io import fits
from astropy import wcs as wcs
from drizzle import drizzle as drz
from drp.stepmiparent import StepMIParent
from drp.datafits import DataFits

class StepCoadd(StepMIParent):
    """ Stone Edge Pipeline Step Master Bias Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '1.1' # pipe step version
    
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
                               - options are square, point, gaussian, smoothing, tophat'])
        self.paramlist.append(['pixfrac', 1.,
                               'The fraction of an output pixel(s) that an input pixel\'s flux is confined to'])
        self.paramlist.append(['resolution', 1.,
                               'Pixel scale divisor for output image (higher gives more resolution, lower gives less)'])
        self.paramlist.append(['pad', 0,
                               'Extra padding outside maximum extent of inputs'])
        self.paramlist.append(['fillval', '0',
                               'Value for filling in the area(s) in the output where there is no input data'])
        self.paramlist.append(['drizzleweights','exptime',
                               'How each input image should be weighted when added to the output \
                               - options are exptime, expsq and uniform'])
        self.paramlist.append(['outputangle', 0,
                               'angular deviation from north-up in degrees (currently broken - do not vary)'])

    def run(self):
        """ Runs the mosaicing algorithm. The self.datain is run
        through the code, the result is in self.dataout.
        """
        #calculate platescale of first input image
        det = np.linalg.det(wcs.WCS(self.datain[0].header).wcs.cd)
        pscale = np.sqrt(np.abs(det))*3600./self.getarg('resolution')
        
        #calculations necessary for updating wcs information
        px = []
        py = []
        
        for f in self.datain:
            px.extend(wcs.WCS(f.header).calc_footprint()[:,0])
            py.extend(wcs.WCS(f.header).calc_footprint()[:,1])
        x0 = (max(px)+min(px))/2.
        y0 = (max(py)+min(py))/2.
        sx = (max(px)-min(px))*np.cos(y0/180*np.pi) # arcsec
        sy = (max(py)-min(py)) # arcsec
        size = (sx*3600+self.getarg('pad')*2, sy*3600+self.getarg('pad')*2)
        xpix = size[0]//pscale
        ypix = size[1]//pscale
        cdelt = [pscale/3600.]*2
        
        #create self.dataout and give it a copy of an input's header
        self.dataout = DataFits(config = self.config)
        self.dataout.header = self.datain[0].header.copy()
        
        #update header wcs information
        self.log.info('Creating new WCS header')
        
        self.dataout.header['CRPIX1'] = xpix/2
        self.dataout.header['CRPIX2'] = ypix/2
        self.dataout.header['CRVAL1'] = x0
        self.dataout.header['CRVAL2'] = y0
        self.dataout.header['CD1_1'] = -cdelt[0]
        self.dataout.header['CD1_2'] = self.dataout.header['CD2_1'] = 0.
        self.dataout.header['CD2_2'] = cdelt[1]
        self.dataout.header['NAXIS1'] = int(xpix)
        self.dataout.header['NAXIS2'] = int(ypix)
        self.dataout.header['CTYPE1'] = 'RA---TAN-SIP'
        self.dataout.header['CTYPE2'] = 'DEC--TAN-SIP'
        self.dataout.header['RADESYS'] = 'ICRS'
        self.dataout.header['EQUINOX'] = 2000
        self.dataout.header['LATPOLE'] = self.datain[0].header['CRVAL2']
        self.dataout.header['LONPOLE'] = 180
        self.dataout.header['PIXASEC'] = pscale
        
        theta_rad = np.deg2rad(self.getarg('outputangle'))
        rot_matrix = np.array([[np.cos(theta_rad), -np.sin(theta_rad)], 
                        [np.sin(theta_rad),  np.cos(theta_rad)]])
        rot_cd = np.dot(rot_matrix, np.array([[self.dataout.header['CD1_1'], 0.],[0., self.dataout.header['CD2_2']]]))
        for i in [0,1]:
            for j in [0,1]:
                self.dataout.header['CD{0:d}_{1:d}'.format(i+1, j+1)] = rot_cd[i,j]
        
        #check drizzle arguments
        if self.getarg('kernel') == 'smoothing':
            kernel = 'lanczos3'
        elif self.getarg('kernel') in ['square', 'point', 'gaussian', 'tophat']:
            kernel = self.getarg('kernel')
        else:
            self.log.error('Kernel name not recognized, using default')
            kernel = 'square'
        if self.getarg('drizzleweights') == 'uniform':
            driz_wt = ''
        elif self.getarg('drizzleweights') in ['exptime', 'expsq']:
            driz_wt = self.getarg('drizzleweights')
        else:
            self.log.error('Drizzle weighting not recognized, using default')
            driz_wt = ''
                        
        #create drizzle object and add input images
        fullwcs = wcs.WCS(self.dataout.header)
        self.log.info('Starting drizzle')
        driz = drz.Drizzle(outwcs = fullwcs, pixfrac=self.getarg('pixfrac'), \
                           kernel=kernel, fillval=self.getarg('fillval'), wt_scl=driz_wt)
        for f in self.datain:
            self.log.info('Adding %s to drizzle stack' % f.filename)
            driz.add_image(f.imgdata[0], wcs.WCS(f.header))
        
        #creates output fits file from drizzle output
        self.dataout.imageset(driz.outsci)
        self.dataout.imageset(driz.outwht,'OutWeight', self.dataout.header)
        self.dataout.filename = self.datain[0].filename

        #add history
        self.dataout.setheadval('HISTORY','Coadd: %d files combined with %s kernel, pixfrac %f at %f times resolution' \
                                % (len(self.datain), kernel, self.getarg('pixfrac'), self.getarg('resolution')))

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
    2019-05-22 New step created for mosaicing and combining images - Matt Merz
    2019-05-24 Minor changes, notably no longer utilizing Grizli for making expanded WCS - Matt Merz
"""
