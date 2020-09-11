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
    """ Pipeline Step Object to calibrate Bias/Dark/Flat files
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
        ### Run Source Extractor
        # Make sure input data exists as file
        if not os.path.exists(self.datain.filename) :
            self.datain.save()
        
        #Open data out of fits file for use in SEP
        psimage = self.datain.image
        image = psimage.byteswap().newbyteorder()

        #Set values for variables used later
        #These variables are used for the background analysis. bw and bh I found just testing various numbers
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

        imsubmed = np.nanmedian(image_sub)
        imsubmad = mad_std(image_sub)



		#Create variables that are used during source Extraction
        extract_thresh1 = 20
        extract_thresh2 = 2
        extract_err = bkg_rms

        #Extract sources from the subtracted image
        sources = sep.extract(image_sub, extract_thresh1, err=extract_err)
        lowthresh= sep.extract(image_sub, extract_thresh2, err=extract_err)


        ## Sort by descending isophotal flux. (Taken from Dr. Harper's SEP Notebook)
        ind = np.argsort(sources['flux'])
        reverser = np.arange(len(ind) - 1,-1,-1)
        rev_ind = np.take_along_axis(ind, reverser, axis = 0)
        objects = np.take_along_axis(sources, rev_ind, axis = 0)

        indlow = np.argsort(lowthresh['flux'])
        reverserlow = np.arange(len(indlow) - 1,-1,-1)
        rev_indlow = np.take_along_axis(indlow, reverserlow, axis = 0)
        lowobjects = np.take_along_axis(lowthresh, rev_indlow, axis = 0)

        ###Do basic uncalibrated measurments of flux for use in step astrometry. 
        '''
        First we calculate flux using Ellipses. In order to do this we need to calculate
        the Kron Radius for the ellipses Extract identified using the ellipse 
        parameters it gives. 
        R is equal to 6 as that is the default used in Source Extractor
        '''
        kronrad, krflag = sep.kron_radius(image_sub, objects['x'], objects['y'], 
        	objects['a'], objects['b'], objects['theta'], r=6.0)

        lowkron, lowkrflag= sep.kron_radius(image_sub, lowobjects['x'], lowobjects['y'], 
            lowobjects['a'],lowobjects['b'], lowobjects['theta'], r=6.0)

        '''
        This is the equivalent of the flux_auto rmin param for Source Extractor. 
        It is 3.5 in the param file from the original version of the step
        '''
        r_min=3.5
        
        

        #Using this Kron radius we calculate the flux
        #This is equivallent to FLUX_AUTO in SExtractor
        flux_elip, fluxerr_elip, flag = sep.sum_ellipse(image_sub, objects['x'], objects['y'], objects['a'], 
                                      objects['b'], objects['theta'], r= 2.5*kronrad, err=bkg_rms,
                                      subpix=1)
        flux_elow, fluxerr_elow, flag = sep.sum_ellipse(image_sub, lowobjects['x'], lowobjects['y'], lowobjects['a'], 
                                      lowobjects['b'], lowobjects['theta'], r= 2.5*lowkron, err=bkg_rms,
                                      subpix=1)

		#Then we calculate it using Circular Apetures
        #This will be used to remove sources that are too elipitical
        #It is equivalent to FLUX_APER in Sextractor
        flux_circ, fluxerr_circ, lflag = sep.sum_circle(image_sub,
			objects['x'], objects['y'], r=2.5, err = bkg_rms, subpix=1)

        flux_lowc, fluxerr_lowc, cflag = sep.sum_circle(image_sub,
            lowobjects['x'], lowobjects['y'], r=2.5, err = bkg_rms, subpix=1)


        #Now we want to calculat the Half-flux Radius.

        dx = (objects['xmax'] - objects['xmin']) / 2
        dy = (objects['ymax'] - objects['ymin']) / 2

        dxl = (lowobjects['xmax'] - lowobjects['xmin']) / 2
        dyl = (lowobjects['ymax'] - lowobjects['ymin']) / 2


        rmax = np.sqrt(dx*dx + dy*dy)
        rmaxlow= np.sqrt(dxl*dxl + dyl*dyl)
        #Frac is the amount of flux we want for the radius, since we want half flux it is .5
        frac=0.5
        rh, rh_flag = sep.flux_radius(image_sub, objects['x'], objects['y'], rmax, frac)
        rhl, rhl_flag = sep.flux_radius(image_sub, lowobjects['x'], lowobjects['y'], rmaxlow, frac)

   
        #Sort the individual arrays so that the final table is sorted by flux
        #create sorting index by using flux. This is for low threshold
        indl = np.argsort(flux_elow)
        reverserl = np.arange(len(indl) - 1,-1,-1)
        rev_indl = np.take_along_axis(indl, reverserl, axis = 0)
        flux_elow = np.take_along_axis(flux_elow, rev_indl, axis = 0)
        #now apply it to all the axis
        fluxerr_elow = np.take_along_axis(fluxerr_elow, rev_indl, axis = 0)
        lowobjects = np.take_along_axis(lowobjects, rev_indl, axis = 0)
        rhl = np.take_along_axis(rhl, rev_indl, axis = 0)

        #now for high threshold

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
        elonglow = (lowobjects['a']/lowobjects['b'])<elim
        seo_SNL = (elonglow) & ((flux_elow/fluxerr_elow)<1000) & (fluxerr_elow != 0)


        self.log.debug('Selected %d high thershold stars from Source Extrator catalog' % np.count_nonzero(seo_SN))
        self.log.debug('Selected %d low thershold stars from Source Extrator catalog' % np.count_nonzero(seo_SNL))

        #Calculate mean RH, its STD, and mean Elongation to report in header
        rhmean, rhstd = np.nanmean(rhl[seo_SNL]), mad_std(rhl[seo_SNL], ignore_nan = True)
        elmean= np.nanmean(lowobjects['a'][seo_SNL]/lowobjects['b'][seo_SNL])


       




        
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



        # Now lets make a table using the lower threshold
        grid = []
        numlow = np.arange(1, len(lowobjects['x'][seo_SNL]) + 1 )
        grid.append(fits.Column(name='ID', format='D',
                                array=numlow))
        grid.append(fits.Column(name='X', format='D',
                                array=lowobjects['x'][seo_SNL],
                                unit='pixel'))
        grid.append(fits.Column(name='Y', format='D',
                                array=lowobjects['y'][seo_SNL],
                                unit='pixel'))
        grid.append(fits.Column(name='Uncalibrated Flux', format='D',
                                array=flux_elow[seo_SNL],
                                unit='flux'))
        grid.append(fits.Column(name='Uncalibrated Fluxerr', format='D',
                                array=fluxerr_elow[seo_SNL], unit='flux'))
        grid.append(fits.Column(name='Half-light Radius', format='D',
                                array=rhl[seo_SNL], unit='pixel'))


        # Make table
        c = fits.ColDefs(cols)

        clow=fits.ColDefs(grid)

        sources_table = fits.BinTableHDU.from_columns(c)

        lowsource_table= fits.BinTableHDU.from_columns(clow)

        
        ### Make output data
        # Copy data from datain
        self.dataout = self.datain
        self.dataout.setheadval ('RHALF',rhmean, 'Mean half-power radius of stars (in pixels)') 
        self.dataout.setheadval ('RHALFSTD', rhstd, 'STD of masked mean of half-power radius')
        self.dataout.setheadval ('ELONG',elmean, 'Mean elongation of accepted sources')
        self.dataout.tableset(sources_table.data,'High Threshold Sources',sources_table.header)
        self.dataout.tableset(lowsource_table.data, 'Low Threshold Sources', lowsource_table.header)


        
    
        

       
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
