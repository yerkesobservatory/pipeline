#!/usr/bin/env python
"""
    Pipestep FluxCal

    This module defines the pipeline step to flux calibrate data files.
    The pipe step looks at previously extracted sources from the LTS table 
    and compares itentified sources with values from the StSci guide star catalog.

    Author: Amanda Pagul / Marc Berthoud / Will Rehmus

    export PYTHONPATH=/Users/berthoud/edu/outreach/Telescopes/pipeline/source

"""
import os # os library
# import sys # sys library
import numpy as np # numpy library
from scipy.optimize import minimize # scipy library
# import string # string library
import logging # logging object library
# import subprocess # running a subprocess library
import requests # http request library
# import astropy.table # Read astropy tables
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord # To make RA/Dec as float
from astropy import units as u # To help with SkyCoord
from astropy import wcs
from astropy.time import Time
import matplotlib # to make plots
matplotlib.use('Agg') # Set pixel image
import matplotlib.pyplot as plt # matplotlib library for plotting
from darepype.drp import StepParent # pipestep stepparent object

class StepFluxCal(StepParent):
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
    """

    stepver = '0.2' # pipe step version


    def setup(self):
        """ ### Names and Parameters need to be Set Here ###
            Sets the internal names for the function and for saved files.
            Defines the input parameters for the current pipe step.
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
        ### Set Misc.
        # Set name of the pipeline reduction step
        self.name='fluxcal'
        
        # Set shortcut for pipeline reduction step and identifier for saved file names.
        self.procname = 'FCAL'
        
        # Set Logger
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        
        
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        
        # Append parameters
        self.paramlist.append(['filtermap', 'g-band=g|r-band=r|i-band=i|z-band=z',
                               'Mapping from telescope filter names to SDSS filter names. ' +
                               'Data from multiple filters can be calibrated using the same band. ' +
                               'Example: "telg=g|telr=r|telclear=r"'])
        self.paramlist.append(['zeropercent', 30.0,
                               'Percentile for BZERO value - Currently unused'])
        self.paramlist.append(['fitplot',False,
                               'Flag for making png plot of the fit'])
        self.log.debug('Setup: done')


    def run(self):
        """ Runs the calibrating algorithm. The calibrated data is
            returned in self.dataout
        """
        ### Import sources from table created during source extraction
        sep_catalog = self.datain.tableget('LTS')
        X = sep_catalog['X']
        Y = sep_catalog['Y']
        seo_Mag = -2.5*np.log10(sep_catalog['Uncalibrated Flux'])
        seo_MagErr = (2.5/np.log(10)*(sep_catalog['Uncalibrated Flux Error']/sep_catalog['Uncalibrated Flux']))
        w = wcs.WCS(self.datain.header)
        ra, dec = w.all_pix2world(X,Y, 0)
        seo_radec = SkyCoord(ra = ra*u.deg, dec = dec*u.deg)


        ### Query and extract reference sources from Guide Star Catalog
        # Get RA & Dec of image center
        ra_center =  self.datain.getheadval('RA').split(':')
        dec_center = self.datain.getheadval('DEC').split(':')
        ra_cent =  ' '.join([str(s) for s in ra_center])
        dec_cent = ' '.join([str(s) for s in dec_center])
        center_coordinates = SkyCoord(ra_cent + ' ' + dec_cent, unit=(u.hourangle, u.deg))
        self.log.debug('Using RA/Dec = %s / %s' % (center_coordinates.ra, center_coordinates.dec))
        
        # Query guide star catalog 2 with center coordinates
        gsc2_query = 'https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?'
        gsc2_query += 'RA='+str(center_coordinates.ra.value)
        gsc2_query += '&DEC='+str(center_coordinates.dec.value)
        gsc2_query += '&DSN=+&FORMAT=CSV&CAT=GSC241&SR=0.5&'
        self.log.debug('Running URL = %s' % gsc2_query)
        gsc2_result = requests.get(gsc2_query)
        gsc_table = ascii.read(gsc2_result.text)
        
        # Convert filter names from SEO to SDSS convention (i.e. 'i-band' -> 'i')
        filter_map = self.getarg('filtermap').split('|')
        filter_name = filter_tel = self.datain.getheadval('FILTER')
        for fil in filter_map: 
            entry = fil.split('=')
            if entry[0] == filter_tel:
                try:
                    filter_name = entry[1]
                except:
                    self.log.error("Badly formatted filter mapping. No '=' after %s"
                                   % filter_tel)
        table_filter = 'SDSS'+filter_name+'Mag'
        table_filter_err = 'SDSS'+filter_name+'MagErr'
        
        # Filter to only include guide stars with 0 < magnitude < 22 in the imaged band
        query_mask = (gsc_table[table_filter]<22) & (gsc_table[table_filter]>0)
        GSC_ID = gsc_table['hstID'][query_mask]
        GSC_Mag = gsc_table[table_filter][query_mask]
        GSC_MagErr = gsc_table[table_filter_err][query_mask]
        GSC_RA = gsc_table['ra'][query_mask]
        GSC_DEC = gsc_table['dec'][query_mask]
        GSC_radec = SkyCoord(ra=GSC_RA*u.deg, dec=GSC_DEC*u.deg)
        self.log.debug('Received %d entries from Guide Star Catalog' % len(GSC_RA))
        

        ### Match Guide Star Catalog data with Source Extractor data
        idx, d2d, d3d = GSC_radec.match_to_catalog_sky(seo_radec)
        
        # Filter to only include pairs within one pixel of each other
        binning = self.datain.getheadval('XBIN')
        pixel_scale = self.datain.getheadval('TELSCALE') * self.datain.getheadval('PIXSIZE1') / 1000
            # Factor of 1000 to convert between microns and mm
        dist_value = 1*pixel_scale*binning/3600. #Maximum distance is 1 pixel
            # Factor of 3600 to convert between arcsec and deg
        mask = d2d.value<dist_value
        if(np.sum(mask) < 2):
            self.log.warning('Only %d sources match between image and guide star catalog, fit may not work' %
                          np.sum(mask))
        self.log.debug('Distance_Value = %f, Min(distances) = %f, Mask length = %d' %
                       (dist_value, np.min(d2d.value), np.sum(mask)))
        

        ### Calculate the fit correction between the guide star and the extracted values
            # The fit finds m_ml and b_ml where seo_Mag = b_ml + m_ml * GSC_Mag
            # then takes the value of this fit at the median magnitude star as the offset 
        # Make lambda function to be minimized
        nll = lambda *args: residual(*args)
        
        # Get combined errors
        eps_data = np.sqrt(GSC_MagErr**2+seo_MagErr[idx]**2)
        
        # Make estimate for intercept to give as initial guess
        b_ml0 = np.ma.median(seo_Mag[idx][mask]-GSC_Mag[mask])
        self.log.debug('Offset guess is %f mag' % b_ml0)
        
        # Calculate distance from that guess and get median
        guessdistances = np.abs(b_ml0 - (seo_Mag[idx] - GSC_Mag))
        guessdistmed = np.ma.median(guessdistances[mask])
        
        # Update mask to ignore sources with SEO magnitudes far from the expected value 
            # based on the initial offset guess 
            # Uses 5 * median difference as the cutoff
        mask = np.logical_and(mask, guessdistances < 5 * guessdistmed)
        arguments = (GSC_Mag[mask], seo_Mag[idx][mask], eps_data[mask])
        self.log.debug(f'Initial guess residuals = {residual((1, b_ml0), *arguments)}')
        self.log.debug('Median of distance to guess = %f, Mask length = %d' %
                       (guessdistmed, np.sum(mask)))
        
        # Solve linear equation
        result = minimize(nll, [1, b_ml0], args=arguments)
        m_ml, b_ml = result["x"]
        self.log.info('Fitted offset is %f mag, fitted slope is %f' % (b_ml, m_ml))
        self.log.debug(f'Best fit residuals = {residual((m_ml, b_ml), *arguments)}')
        
        # Calculate b_ml_corr 
            # find the value of the best fit at the median magnitude source, takes that as true calibration
        b_ml_corr = b_ml + (m_ml-1) * np.ma.median(GSC_Mag[mask])
        self.log.info('Corrected offset is %f mag' % b_ml_corr)
        self.log.debug(f'Corrected offset residuals = {residual((1, b_ml_corr), *arguments)}')


        ### Make table with data which was fit
        # Get order, sorted by GSC magnitude
        sort = np.argsort(GSC_Mag[mask])
        
        # Collect data columns
        cols = []
        cols.append(fits.Column(name='GSCII ID', format='10A', array=GSC_ID[mask][sort]))
            # hstID from GSC database
            # Can be searched at https://gsss.stsci.edu/webservices/GSC2/WebForm.aspx
            # If the above fails, replace the id number in the following query:
                # https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?FORMAT=HTML&id=NBQI004317
        cols.append(fits.Column(name='RA (GSC)', format='D', array=GSC_RA[mask][sort],
                                unit='deg'))
        cols.append(fits.Column(name='Dec (GSC)', format='D', array=GSC_DEC[mask][sort],
                                unit='deg'))
        cols.append(fits.Column(name='Diff_Deg', format='D', array=d2d[mask][sort],
                                unit='deg'))
        cols.append(fits.Column(name='GSC_Mag', format='D', array=GSC_Mag[mask][sort], 
                                unit='magnitude'))
        cols.append(fits.Column(name='Img_Mag', format='D',
                                array=seo_Mag[idx][mask][sort],
                                unit='magnitude'))
        cols.append(fits.Column(name='Error', format='D', array=eps_data[mask][sort],
                                unit='magnitude'))
        
        # Make table
        c = fits.ColDefs(cols)
        fitdata_table = fits.BinTableHDU.from_columns(c)
        
        
        ### Make output data
        # Copy data from datain
        self.dataout = self.datain
            # This copies the pointer, so modifying dataout will also modify datain
        
        # Add Photometric Zero point magnitude
        self.dataout.setheadval('PHTZPRAW', -b_ml_corr, 'Photometric zeropoint for RAW data')
        self.dataout.setheadval('PTZRAWER', 0.0, 'Uncertainty of the RAW photometric zeropoint')
        self.dataout.setheadval('PHOTZP', 8.9,  'Photometric zeropoint MAG=-2.5*log(data)+PHOTZP')
        self.dataout.setheadval('BUNIT', 'Jy/pixel', 'Units for the data')
        # bzero = np.nanpercentile(self.dataout.image,self.getarg('zeropercent'))
            # This is an incomplete attempt at computing b_zero. It is unused
        
        # Add flux calibrated image (using b_ml_corr as calibration)
        bscale = 3631. * 10 ** (b_ml_corr/2.5) 
            # 3631 is the zero-point for AB magnitudes in Janskys (Jy)
            # b_ml_corr accounts for scaling difference between true magnitude 
                # and the magnitude derived from the BDF image counts
            # bscale is thus a conversion factor between pixel counts and Jy
        self.dataout.image = self.datain.image * bscale
        
        # Add flux calibrated background-subtracted image (same calibration)
        dataname = "CALIMSUB"
        self.dataout.imageset(self.datain.imageget("IMSUB") * bscale, imagename=dataname)
        self.dataout.setheadval('HISTORY', 'CALIMSUB',
                                dataname=dataname)
        
        # Add fitdata table
        self.dataout.tableset(fitdata_table.data,'Fit Data',fitdata_table.header)


        ### Add calibrated data to LTS table
        self.dataout.tableaddcol(colname="RA", array=ra, tablename="LTS")
        self.dataout.tableaddcol(colname="Dec", array=dec, tablename="LTS")
        lts_mag = sep_catalog['Uncalibrated Flux'] * bscale
        lts_mag_err = sep_catalog['Uncalibrated Flux Error'] * bscale
        self.dataout.tableaddcol(colname="Calibrated Flux", array=lts_mag, tablename="LTS")
        self.dataout.tableaddcol(colname="Calibrated Flux Error", array=lts_mag_err, tablename="LTS")

        ### If requested make a plot of the fit and save as png
        if self.getarg('fitplot'):
            # Set up plot
            plt.figure(figsize=(10,7))
            ax = plt.axes()

            # Plot the datapoints
            reject_mask = np.logical_and(d2d.value<dist_value, ~mask)
            ax.errorbar(GSC_Mag[reject_mask],
                        seo_Mag[idx][reject_mask],
                        yerr=np.sqrt(eps_data[reject_mask]**2),fmt='x',linestyle='none',
                        c="xkcd:red", label='Rejected Data')
            ax.errorbar(GSC_Mag[mask],seo_Mag[idx][mask],
                        yerr=np.sqrt(eps_data[mask]**2),fmt='.',linestyle='none',
                        c="k", label='Accepted Data')
            
            # Fix plot scale to only consider datapoints
            ax.autoscale_view()
            ax.autoscale(False)
            
            # Plot 5sigma error range
            gmin, gmax = ax.get_xlim()
            ax.fill([gmin,gmin,gmax,gmax],[gmin+b_ml0-5*guessdistmed, gmin+b_ml0+5*guessdistmed,
                                            gmax+b_ml0+5*guessdistmed, gmax+b_ml0-5*guessdistmed], 
                     c="xkcd:green", alpha=0.2, label='Acceptance Range')
            
            # Plot fits
            ax.plot([gmin, gmax],[gmin, gmax]+b_ml0, label='Initial Guess', c="xkcd:green")
            ax.plot([gmin, gmax],np.array([gmin, gmax])*m_ml+b_ml, label='Linear fit', c="xkcd:orange")
            ax.plot([gmin, gmax],[gmin, gmax]+b_ml_corr, label='Corrected offset (Used)', c="xkcd:blue")
            
            ax.legend()
            ax.set_ylabel('Source extrator magnitude (SEO observation data)')
            ax.set_xlabel('Guide Star Catalog magnitude')
            ax.set_title('Calibration Fit for file\n' + os.path.split(self.dataout.filename)[1])
            # Plot the fit
            # Axis and labels
            # Save the image
            pngname = self.dataout.filenamebegin + 'FCALplot.png'
            plt.savefig(pngname)
            plt.clf()
            self.log.debug('Saved fit plot under %s' % pngname)


def residual(params, x, data, errors):
    """ Fitting function for scipy.optimize.minimize
    """
    m, c = params
    model = m*x+c
    inv_sigma2 = 1.0/(errors**2)
    return 0.5*(np.sum(((data-model)**2)*inv_sigma2))


if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
        python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
        --config=ConfigFilePathName.txt : name of the configuration file
        --test : runs the functionality test i.e. pipestep.test()
        --loglevel=LEVEL : configures the logging output for a particular level
    """
    StepFluxCal().execute()

'''HISTORY:
2018-09-019 - Started based on Amanda's code. - Marc Berthoud
2026-06 - Updated for increased clarity and ease of information access - Will Rehmus
'''
