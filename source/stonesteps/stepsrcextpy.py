#!/usr/bin/env python
"""
    Pipestep SrcExtPy

    This module defines the pipeline step to extract sources from FITS image data.
    Normally, it would be used after bias/dark/flat correction and bad-pixel masking.
    However, it can operate on any image file.
    It uses the python library SEP to extract a list of sources and associated data.
    In addition to the tables and image described below, keyword data are added to
    the header of the DataFits output object for the mean elongation of isolated stars,
    the half-light radii of isolated stars, and for their STDs of the mean.

    Inputs: 
        A FITS image file and input parameters defined in a config file.
    Table outputs (in DataFits output object):
        A table of extracted sources (from sep.extract() with additional columns with
            other data derived from the SEP output).
        A table of sources deemed to be stars, along with their fluxes (in ADU) and
            x and y positions (in physical pixel coordinates). The table is ordered
            by descending flux and is optimized to provide a list of stars to be
            analyzed by a subsequent astrometry step to determine a WCS.
        ASCII versions of the second and third tables above in CSV format,
            saved as separate files.
    Image outputs (HDUs in DataFits output object):
        The primary image HDU passed through from the input file.
        An image which is the difference of the original image and the background image.
    Optional image HDU outputs:
        An image of the background provided by the function sep.Background().
        An image of the rms background noise.

    For more info check out the read the docs for SEP: https://sep.readthedocs.io
      Also sep uses Sextractor, more info at https://astromatic.github.io/sextractor
    
    Bugfixes:
    - Exception: object deblending overflow: limit of 1024 sub-objects reached while
                 deblending. Decrease number of deblending thresholds or increase the
                 detection threshold.
      => Increase the ext_thresh parameter
    
    Authors: Amanda Pagul / Daniel Sharkey/ Al Harper/ Alexa Bukowski/ Alexandra Masegian
"""
import os # os library
import sys # sys library
import numpy as np # numpy library
#import scipy # scipy library
import string # string library
import logging # logging object library
#import subprocess # running a subprocess library
#import requests # http request library
import astropy.table # Read astropy tables
import sep # Extracts Sources and Calculates Flux
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import mad_std
from darepype.drp import StepParent # pipestep stepparent object
from darepype.drp.datafits import DataFits


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
        '''
        NAMES:
        Set the internal name for the function and the procedure identifier to be incorporated in output file names. 
        name: The name should be all lower-case and should be identical to 
            the module name (the module name is used to identify pipe steps
            in the config file, so if it doesn’t match case an error will occur). 
        procname: The procname should be all upper-case and, if possible, 
            restricted to three or four alphabetic characters(the code would 
            automatically capitalize, but it’s good form to also capitalize it here).
        '''
        # Name of the pipeline reduction step
        self.name='srcextpy'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'SEP'
        # Set Logger for this pipe step
        self.log = logging.getLogger('pipe.step.%s' % self.name)
        ### Set Parameter list
        '''
        PARAMETERS:
        Parameters are stored in a list containing the following information:
        name: The name for the parameter. This name is used when calling the 
            pipe step from command line or python shell. It is also used to 
            identify the parameter in the pipeline configuration file.
        default: A default value for the parameter. 
            If nothing, set ‘ ’ for strings 0 for integers, and 0.0 for floats.
        help: A short description of the parameter.
        '''
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['sourcetable',False,
                               'Flag for making txt table of all sources'])
        self.paramlist.append(['sourcetableformat','csv',
                               'txt table format (see astropy.io.ascii for options)'])
        self.paramlist.append(['bkg_maskthreshold', 0.0,
                                'mask threshold for background detection'])
        self.paramlist.append(['bkg_wh', [16,16],
                                'background box width and height for background detection'])
        self.paramlist.append(['bkg_filterwh', [3,3],
                                'filter width and height for background detection'])
        self.paramlist.append(['bkg_fthreshold', 0.0,
                                'filter threshold for background detection'])
        self.paramlist.append(['ext_thresh', 2.0,
                                'extraction threshold for source extration'])
        self.paramlist.append(['ext_bfactor', 10.0,
                                'brightness factor for creating highlevel threshold'])
        self.paramlist.append(['ext_deblend', 256,
                                'deblend threshold for source extraction'])
        self.paramlist.append(['phot_kronf', 2.5,
                                'factor multiplied into kronrad to get radius for integration'])
        self.paramlist.append(['save_background', True,
                                'option to save the background as a seprate hdu'])
        #self.paramlist.append(['byte_swap', False,
        # 	                    'says if the bytes should be swapped or not for the image'])
        # confirm end of setup
        self.log.debug('Setup: done')

    def run(self):
        """ Runs the source extraction algorithm. The extracted object data are
            returned in a table. The input image, and background-subtracted image
            are stored in HDUs in self.dataout. If desired, the background image
            is also stored in an HDU.
        """
        
        '''Prepare data. Byte swap and convert to float64, if necessary.'''
        
        binning = self.datain.getheadval('XBIN')
        psimage = self.datain.image.copy()   # Get image data from fits file.
        psimage = psimage + 1000.
        
        # Byteswap if required
        maxexp = np.max(np.abs(np.frexp(psimage)[1]))
        self.log.debug('Initial max exponent value: %d' % maxexp)
        self.log.debug('Isnative: %s' % repr(psimage.dtype.isnative))
        # if maxexp > 19:
        if not psimage.dtype.isnative:
            self.log.debug("Performing byte swap")
            image = psimage.byteswap(inplace=True)
            print('Byte swap performed')
        else:
            self.log.debug("No byte swap required")
            image = psimage
            print('No byte swap required')
        # Convert to float64, if necessary.
        self.log.debug("Initial data type: %s" % psimage.dtype)
        if psimage.dtype != np.float64:
            self.log.debug("Converting to float64")
            image = image.astype(np.float64)
            print('Conversion to float64 required')
                        
        '''Get variable values used to compute background (from paramlist).'''
        
        maskthresh = self.getarg('bkg_maskthreshold')
        backwh= self.getarg('bkg_wh')
        filwh= self.getarg('bkg_filterwh')
        bw, bh = backwh[0], backwh[1]
        fw, fh = filwh[0], filwh[1]
        fthresh = self.getarg('bkg_fthreshold')

        '''Create the background image and its error image.'''
        
        bkg = sep.Background(image, maskthresh=maskthresh,bw=bw, bh=bh, fw=fw,
            fh=fh, fthresh=fthresh) 
        bkg_image=bkg.back()
        bkg_rms =bkg.rms()
        self.log.debug('Background image global')
        if bkg.globalback < np.nanmin(image) or bkg.globalback > np.nanmax(image):
            self.log.warn('Background has out of bounds values - image may not reduce')
        print('globalback =', bkg.globalback)
            
        '''Subtract the background from the image. Calculate its median and mad.'''
        
        image_sub = image - bkg_image
        imsubmed = np.nanmedian(image_sub)
        imsubmad = mad_std(image_sub)
        
        ''' Create variable values needed for source extraction.'''
        
        extract_thresh = self.getarg('ext_thresh')
        deblend_nthresh = self.getarg('ext_deblend')
        kfactor = self.getarg('phot_kronf')
        extract_err = bkg_rms
        
        '''Extract sources from the background-subtracted image.'''
        
        sources = sep.extract(image_sub, extract_thresh, err=extract_err, deblend_nthresh= deblend_nthresh)
        print('Length of sources =', len(sources))

        ''' Eliminate very small sources.'''
                
        alow, blow = 1.25, 1.25
        smalls = []
        for i, s in enumerate(sources):
            if s['a']<-alow or s['b']<=blow:
                smalls.append(i)
        print('length of smalls =', len(smalls))
        sources = np.delete(sources, smalls, 0)
        print('length of sources after removing smalls =', len(sources))
        
        '''
        Delete any rows containing NaNs in a or b columns. This code was added
        to prevent "invalid aperture parameters" errors which crashed the pipeline.
        '''
        
        nan_rows = []
        for i, s in enumerate(sources):
            if np.isnan(s['a']) or np.isnan(s['b']):
                nan_rows.append(i)
        print('length of nan_rows =', len(nan_rows))
        sources = np.delete(sources, nan_rows, 0)
        print('length of sources after removing rows with NaNs =', len(sources))
        
        '''Sort sources by descending isophotal flux'''

        ind = np.argsort(sources['flux'])
        reverser = np.arange(len(ind) - 1,-1,-1)
        rev_ind = np.take_along_axis(ind, reverser, axis = 0)
        objects = np.take_along_axis(sources, rev_ind, axis = 0)
        
        '''
        The following code prevents floating point errors arising from theta
        values that are too close to pi/2.
        '''
        objects['theta'] = np.where(abs(objects['theta'] - np.pi/2) < 0.001, np.pi/2, objects['theta'])
        
        '''
        Calculate fluxes using sep.sum_ellipse. This is equivalent to FLUX_AUTO
        in SExtractor. Set hard-coded value of R=6, since that is the default
        value used in SExtractor.
        '''
        
        # First, calculate kron radii (needed as input arguments for sep.sum_ellipse)
        kronrad, krflag = sep.kron_radius(image_sub, objects['x'], objects['y'], 
        	objects['a'], objects['b'], objects['theta'], r=6.0)
        	
        ## Print some diagnostics.
        #print( kr[0], krflags[0], len(kr))
        kronrad_notfinite = len(kronrad) - len(np.isfinite(kronrad))
        #print('Number kronrad not finite =', kronrad_notfinite)
        numnans = np.sum(np.isnan(image_sub))
        #print('Number of nans in image_sub =', numnans)
        numnans = np.sum(np.isnan(objects['theta']))
        #print('Number of nans in theta =', numnans) 
        numnans = np.sum(np.isnan(objects['b']))
        #print('Number of nans in b =', numnans)
        #for i in range(len(objects)):
        #    print(i, objects['a'][i], objects['b'][i])
               
        # Now calculate the elliptical fluxes.
        flux_elip, fluxerr_elip, flag = sep.sum_ellipse(image_sub, objects['x'], objects['y'], objects['a'], 
        objects['b'], objects['theta'], r= kfactor*kronrad, err=bkg_rms, subpix=1)
        
        '''
        Calculate circular half-flux radii using sep.flux_radius.
        '''
        # First we need a maximum radius (Rmax) for the integral.
        dx = (objects['xmax'] - objects['xmin']) / 2
        dy = (objects['ymax'] - objects['ymin']) / 2
        rmax = np.sqrt(dx*dx + dy*dy)
        # Now we need the argument for the percentage of flux contained within the radius.
        frac=0.5
        # Now do the calculation.
        rh, rh_flag = sep.flux_radius(image_sub, objects['x'], objects['y'], rmax, frac)

        '''Check for nans in rh. Diagnostic for issues experienced with some files.'''
        
        numnans = np.sum(np.isnan(rh))
        print('Number of nans in rh =', numnans)

        '''
        Sort all the arrays that will go into the output table (including the
        list of extracted sources) by descending elliptical flux.
        '''
        ind = np.argsort(flux_elip)
        reverser = np.arange(len(ind) - 1,-1,-1)
        rev_ind = np.take_along_axis(ind, reverser, axis = 0)
        flux_elip = np.take_along_axis(flux_elip, rev_ind, axis = 0)
        fluxerr_elip = np.take_along_axis(fluxerr_elip, rev_ind, axis = 0)
        rh = np.take_along_axis(rh, rev_ind, axis = 0)
        objects = np.take_along_axis(objects, rev_ind, axis = 0)
        
        '''
        Create data cuts aimed at identifying isolated stars in the object list.
        Select for circular, unsaturated stellar images with good signal to noise.
        List will be passed to StepFluxcal to match with GSC calibration stars,
        '''
        elim=1.5    # Set limit on ellipticity of extracted, unsaturated sources
        # Create cuts
        a2b= (objects['a']/objects['b'])
        semimajor = objects['a'] < 1.0
        semiminor = objects['b'] < 1.0
        smallmoment = (semimajor) & (semiminor)
        elong = a2b<elim
        seo_SN = (elong)  & ((flux_elip/fluxerr_elip)<1000) & (fluxerr_elip > 0) & (flux_elip > 0) 
        self.log.debug('Selected %d low threshold stars from Source Extrator catalog' % np.count_nonzero(seo_SN))
                        
        '''
        Make table with the selected subset of data from SEP
        '''
        # Collect data columns
        cols = []
        num = np.arange(1, len(objects['x'][seo_SN]) + 1 )
        #This is done so that the list is indexed properly
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
        cols.append(fits.Column(name='Uncalibrated Flux Error', format='D',
                                array=fluxerr_elip[seo_SN],
                                unit='flux'))
        cols.append(fits.Column(name='Half-light Radius', format='D',
                                array=rh[seo_SN], unit='pixel'))

        # Make selected subset table
        c = fits.ColDefs(cols)
        sources_table = fits.BinTableHDU.from_columns(c)
        
        #Calculate mean RH, its STD, and mean a/b to report in header
        rhmean = np.nanmean(rh[seo_SN])
        rhstd = np.std(rh[seo_SN])
        elmean= np.nanmean(objects['a'][seo_SN]/objects['b'][seo_SN])

        '''Make output DataFits object (self.dataout)'''
        
        # Copy data from datain
        self.dataout = self.datain.copy()
         
        # Make a table which includes all extracted objects and add to output data
        self.dataout.tableset(objects, tablename='SEP_objects')
        self.dataout.tableaddcol('rh', rh, 'SEP_objects')
        self.dataout.tableaddcol('kflux', flux_elip, 'SEP_objects')
        self.dataout.tableaddcol('a2b', a2b, 'SEP_objects')
        
        # Add table with selected subset of data
        self.dataout.tableset(sources_table.data,'LTS',sources_table.header)

        # Save the background-subtracted image in a second image HDU
        dataname = "IMSUB"
        self.dataout.imageset(image_sub, imagename=dataname)
        self.dataout.setheadval('HISTORY', 'IMSUB',
                                dataname=dataname)
                                
        #If save_background is true, save background image in a third image HDU
        if self.getarg('save_background'):
            dataname = "BACKGROUND"
            self.dataout.imageset(bkg_image, imagename=dataname)
            self.dataout.setheadval('HISTORY', 'BACKGROUND', 
                                    dataname=dataname)
            dataname = "BACKGROUND RMS"
            self.dataout.imageset(bkg_rms, imagename=dataname)
            self.dataout.setheadval('HISTORY', 'BACKGROUND RMS', 
                                    dataname=dataname)                                    
        #Add header keywords
        self.dataout.setheadval ('RHALF',rhmean, 'Mean half-power radius of LTS stars (pixels)') 
        self.dataout.setheadval ('RHALFSTD', rhstd, 'std of half-power radii of LTS stars (pixels)')
        self.dataout.setheadval ('ELONG',elmean, 'Mean elongation ratio of accepted sources')
        self.dataout.setheadval ('ETHRESH', extract_thresh, 'Extraction threshold for SEP_objects table')
       
        '''
        If requested, make a text file with the selected sources list
        in DS9 region file format
        '''
        
        if self.getarg('sourcetable'):
            # Save sources from LTS table as a DS9 *.reg file
            filename = self.dataout.filenamebegin + 'FCALsources.reg'
            with open(filename, 'w+') as f:
                f.write("# Region file format: DS9 version 4.1\n")
                f.write("""global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
                image \n""")
                for i in range(len(objects['x'][seo_SN])):
                    f.write("circle(%.7f,%.7f,5) # text={%i}\n"%(objects['x']
                    [seo_SN][i]+1,objects['y'][seo_SN][i]+1,num[i]))
                    
            # Save the LTS table as a text file
            txtname = self.dataout.filenamebegin + 'FCALsources.txt'
            ascii.write(self.dataout.tableget('LTS'),txtname,
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
2022-12-06 - Updated documentation. Eliminated HTS table. Added work-around for
             "invalid aperture parameters" issue. - Al Harper
'''
