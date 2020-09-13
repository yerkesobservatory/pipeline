#!/usr/bin/env python
"""
    Pipestep SrcExtPy

    This module defines the pipeline step to extract sources from data files.
    The pipe step runs the python library SEP on the data in order to extract
    the sources and calculate a first order approximation of flux and magnitude

    For more info check out the read the docs for SEP: https://sep.readthedocs.io/

   This uses the SEP python library for source extraction.
    Author: Amanda Pagul / Marc Berthoud/ Daniel Sharkey

    Update 8/6/20 by Daniel Sharkey
    Here I have attempted to convert the original Fluxcalsex step into a 
    Sextract step that uses the python library SEP.



"""
import os # os library
import sys # sys library
import numpy as np # numpy library
import scipy # scipy library
import string # string library
import logging # logging object library
import subprocess # running a subprocess library
import requests # http request library
import astropy.table # Read astropy tables
import sep # Extracts Sources and Calculates Flux
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord # To make RA/Dec as float
from astropy import units as u # To help with SkyCoord
from astropy.stats import mad_std
import matplotlib # to make plots
matplotlib.use('Agg') # Set pixel image
import pylab as plt # pylab library for plotting
from lmfit import minimize, Parameters # For brightness correction fit
from darepype.drp import StepParent # pipestep stepparent object


class StepSrcExtPy(StepParent):
    """ Pipeline Step Object to extract sources from image files
    """

    stepver = '0.1' # pipe step version


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
        ### Set Names
        # Name of the pipeline reduction step
        self.name='SrcExtPy'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'sepext'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['sourcetable',False,
                               'Flag for making txt table of all sources'])
        self.paramlist.append(['sourcetableformat','csv',
                               'txt table format (see astropy.io.ascii for options)'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def run(self):
        """ Runs the calibrating algorithm. The calibrated data is
            returned in self.dataout
        """
        ### Preparation
        binning = self.datain.getheadval('XBIN')
        ### Perform Source Extraction
        # Make sure input data exists as file
        if not os.path.exists(self.datain.filename) :
            self.datain.save()
        
        #Open data out of fits file for use in SEP
        psimage = self.datain.image
        image = psimage.byteswap().newbyteorder()

   
        #These variables are used for the background analysis.
        maskthresh = 0.0
        bw, bh = 16, 16
        fw, fh = 3, 3
        fthresh = 0.0

        #Create the background image and it's error
        bkg = sep.Background(image, maskthresh=maskthresh,bw=bw, bh=bh, fw=fw,
        fh=fh, fthresh=fthresh) #have sep determine the background of the image

        bkg_image=bkg.back()

        bkg_rms =bkg.rms()
        #Subtract the background from the image
        image_sub = image - bkg_image
        #Calculate the Median and STD of the Subtracted Image
        imsubmed = np.nanmedian(image_sub)
        imsubmad = mad_std(image_sub)

		#Create variables that are used during source Extraction
        extract_thresh = 2.0
        bright_factor= 10.0
        deblend_nthresh =256
        extract_err = bkg_rms

        #Extract sources from the subtracted image. It extracts a low threshold list and a high threshold list
        sources = sep.extract(image_sub, extract_thresh, err=extract_err, deblend_nthresh= deblend_nthresh)
        sourcesbri= sep.extract(image_sub, extract_thresh*bright_factor, err=extract_err)

        ### Sort sources by descending isophotal flux. (Taken from Dr. Harper's Explore SEP Notebook)
        ind = np.argsort(sources['flux'])
        reverser = np.arange(len(ind) - 1,-1,-1)
        rev_ind = np.take_along_axis(ind, reverser, axis = 0)
        objects = np.take_along_axis(sources, rev_ind, axis = 0)

        indbri = np.argsort(sourcesbri['flux'])
        reverserbri = np.arange(len(indbri) - 1,-1,-1)
        rev_indbri = np.take_along_axis(indbri, reverserbri, axis = 0)
        briobjects = np.take_along_axis(sourcesbri, rev_indbri, axis = 0)

        ###Do basic uncalibrated measurments of flux for use in step astrometry. 
        '''
        First we calculate flux using Ellipses. In order to do this we need to calculate
        the Kron Radius for the ellipses the Extract process identified using the ellipse 
        parameters it gives. 
        R is equal to 6 as that is the default used in Source Extractor
        '''
        kronrad, krflag = sep.kron_radius(image_sub, objects['x'], objects['y'], 
        	objects['a'], objects['b'], objects['theta'], r=6.0)

        brikron, brikflag= sep.kron_radius(image_sub, briobjects['x'], briobjects['y'], 
            briobjects['a'],briobjects['b'], briobjects['theta'], r=6.0)

        '''
        This is the equivalent of the flux_auto rmin param for Source Extractor. 
        It is 3.5 in the param file from the original version of the step which use Sextractor
        '''
        r_min=3.5

        #Using this Kron radius we calculate the flux
        #This is equivalent to FLUX_AUTO in SExtractor
        flux_elip, fluxerr_elip, flag = sep.sum_ellipse(image_sub, objects['x'], objects['y'], objects['a'], 
                                      objects['b'], objects['theta'], r= 2.5*kronrad, err=bkg_rms,
                                      subpix=1)
        flux_ebri, fluxerr_ebri, flag = sep.sum_ellipse(image_sub, briobjects['x'], briobjects['y'], briobjects['a'], 
                                      briobjects['b'], briobjects['theta'], r= 2.5*brikron, err=bkg_rms,
                                      subpix=1)

        #Now we want to calculate the Half-flux Radius. This will be reported later
        dx = (objects['xmax'] - objects['xmin']) / 2
        dy = (objects['ymax'] - objects['ymin']) / 2

        dxb = (briobjects['xmax'] - briobjects['xmin']) / 2
        dyb = (briobjects['ymax'] - briobjects['ymin']) / 2


        rmax = np.sqrt(dx*dx + dy*dy)
        rmaxbri= np.sqrt(dxb*dxb + dyb*dyb)
        '''Frac is the percentage of flux we want contained within the radius,
        since we want half flux radius, frac is .5 '''
        
        frac=0.5
        rh, rh_flag = sep.flux_radius(image_sub, objects['x'], objects['y'], rmax, frac)
        rhb, rhb_flag = sep.flux_radius(image_sub, briobjects['x'], briobjects['y'], rmaxbri, frac)

   
        #Sort the individual arrays so that the final table is sorted by flux
        #create sorting index by using flux. This is for the boosted threshold
        indb = np.argsort(flux_ebri)
        reverserb = np.arange(len(indb) - 1,-1,-1)
        rev_indb = np.take_along_axis(indb, reverserb, axis = 0)
        flux_ebri = np.take_along_axis(flux_ebri, rev_indb, axis = 0)
        #now apply it to all the axis
        fluxerr_ebri = np.take_along_axis(fluxerr_ebri, rev_indb, axis = 0)
        briobjects = np.take_along_axis(briobjects, rev_indb, axis = 0)
        rhb = np.take_along_axis(rhb, rev_indb, axis = 0)

        #now for normal threshold

        ind = np.argsort(flux_elip)
        reverser = np.arange(len(ind) - 1,-1,-1)
        rev_ind = np.take_along_axis(ind, reverser, axis = 0)
        flux_elip = np.take_along_axis(flux_elip, rev_ind, axis = 0)
        #now apply it to all the axis
        fluxerr_elip = np.take_along_axis(fluxerr_elip, rev_ind, axis = 0)
        objects = np.take_along_axis(objects, rev_ind, axis = 0)
        rh = np.take_along_axis(rh, rev_ind, axis = 0)

        
        # Select only the stars in the image: circular image and S/N > 10
        #Establish an elongation limit
        elim=1.5
        #Create cuts
        elongation = (objects['a']/objects['b'])<elim
        seo_SN = (elongation) & ((flux_elip/fluxerr_elip)<1000) & (fluxerr_elip != 0)

        #Now do this for the low threshold sources
        elongbri = (briobjects['a']/briobjects['b'])<elim
        seo_SNB = (elongbri) & ((flux_ebri/fluxerr_ebri)<1000) & (fluxerr_ebri != 0)


        self.log.debug('Selected %d high thershold stars from Source Extrator catalog' % np.count_nonzero(seo_SN))
        self.log.debug('Selected %d low thershold stars from Source Extrator catalog' % np.count_nonzero(seo_SNB))

        #Calculate mean RH, its STD, and mean Elongation to report in header
        rhmean, rhstd = np.nanmean(rh[seo_SN]), mad_std(rh[seo_SN], ignore_nan = True)
        elmean= np.nanmean(objects['a'][seo_SN]/objects['b'][seo_SN])

        
        ### Make table with all data from source extractor
        # Collect data columns
        cols = []
        num = np.arange(1, len(objects['x'][seo_SN]) + 1 )
        cols.append(fits.Column(name='ID', format='D',
                                array=num))
        cols.append(fits.Column(name='X', format='D',
                                array=objects['x'][seo_SN],
                                unit='pixel'))
        cols.append(fits.Column(name='Y', format='D',
                                array=objects['y'][seo_SN],
                                unit='pixel'))
        cols.append(fits.Column(name='Uncalibrated Flux', format='D',
                                array=flux_elip[seo_SN],
                                unit='flux'))
        cols.append(fits.Column(name='Uncalibrated Fluxerr', format='D',
                                array=fluxerr_elip[seo_SN], unit='flux'))
        cols.append(fits.Column(name='Half-light Radius', format='D',
                                array=rh[seo_SN], unit='pixel'))



        # Now lets make a table using the Brighter threshold
        grid = []
        numbri = np.arange(1, len(briobjects['x'][seo_SNB]) + 1 )
        grid.append(fits.Column(name='ID', format='D',
                                array=numbri))
        grid.append(fits.Column(name='X', format='D',
                                array=briobjects['x'][seo_SNB],
                                unit='pixel'))
        grid.append(fits.Column(name='Y', format='D',
                                array=briobjects['y'][seo_SNB],
                                unit='pixel'))
        grid.append(fits.Column(name='Uncalibrated Flux', format='D',
                                array=flux_ebri[seo_SNB],
                                unit='flux'))
        grid.append(fits.Column(name='Uncalibrated Fluxerr', format='D',
                                array=fluxerr_ebri[seo_SNB], unit='flux'))
        grid.append(fits.Column(name='Half-light Radius', format='D',
                                array=rhb[seo_SNB], unit='pixel'))


        # Make table
        c = fits.ColDefs(cols)

        cbri=fits.ColDefs(grid)

        sources_table = fits.BinTableHDU.from_columns(c)

        brisource_table= fits.BinTableHDU.from_columns(cbri)

        
        ### Make output data
        # Copy data from datain
        self.dataout = self.datain
        self.dataout.setheadval ('RHALF',rhmean, 'Mean half-power radius of stars (in pixels)') 
        self.dataout.setheadval ('RHALFSTD', rhstd, 'STD of masked mean of half-power radius')
        self.dataout.setheadval ('ELONG',elmean, 'Mean elongation of accepted sources')
        self.dataout.tableset(sources_table.data,'Low Threshold Sources',sources_table.header)
        self.dataout.tableset(brisource_table.data, 'High Threshold Sources', brisource_table.header)
        self.dataout.setheadval ('EXTRACT_THRESH', extract_thresh, 'Extraction Thershold for Low Thershold Table')
        self.dataout.setheadval ('BRIGHT_FACTOR', bright_factor, 'Multiplier to create High Threshold Table')
       
        ### If requested make a text file with the sources list
        if self.getarg('sourcetable'):

            # Save region file

            filename = self.dataout.filenamebegin + 'FCALsources.reg'
            with open(filename, 'w+') as f:
                f.write("# Region file format: DS9 version 4.1\n")
                f.write("""global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 image\n""")
                for i in range(len(objects['x'][seo_SN])):
                    f.write("circle(%.7f,%.7f,0.005) # text={%i}\n"%(objects['x'][seo_SN][i],objects['y'][seo_SN][i],num[i]))

            # Save the table
            txtname = self.dataout.filenamebegin + 'FCALsources.txt'
            ascii.write(self.dataout.tableget('Low Threshold Sources'),txtname,
                        format = self.getarg('sourcetableformat'))
            self.log.debug('Saved sources table under %s' % txtname)


if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
        python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
        --config=ConfigFilePathName.txt : name of the configuration file
        --test : runs the functionality test i.e. pipestep.test()
        --loglevel=LEVEL : configures the logging output for a particular level
    """
    StepSrcExtPy().execute()

'''HISTORY:
2018-09-019 - Started based on Amanda's code. - Marc Berthoud
'''
