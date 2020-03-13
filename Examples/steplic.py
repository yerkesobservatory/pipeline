""" PIPE STEP LIC - Version 1.0.0
    Line-integral-convolution map plotting code. Takes a (PMP) fits file, 
    convoles it, and produces an intensity map overlayed with polarization
    streams.
    
    @author: Dylan Sukay
    LINE INTEGRAL CONVOLUTION (LIC) BASED ON CABRAL ET AL. (1993) AND IDL CODE 
    FROM FALCETA-GONCALVES
    Translated to Python by Joseph Michail
"""

import numpy as np
import scipy as sp
import logging # logging object library
from drp.stepparent import StepParent
from drp.dataparent import DataParent # Pipeline Data object
from astropy.io import fits
from astropy import wcs as astwcs
import matplotlib
matplotlib.use('Agg') # Set pixel image (necessary for remote/scripted use)
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
#from matplotlib import cm
import matplotlib.colors as colors
#import colorcet
import lic as lic
import aplpy
import os
import sys
from matplotlib.patches import Ellipse, Circle


class StepLIC(StepParent):
    """ HAWC Pipeline Step Region Object
    The object is callable. It requires a valid configuration input
    (file or object) when it runs.
    """
    stepver = '1.0' # pipe step version

    def setup(self):
        """ ### Names and Prameters need to be Set Here ###
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
        self.name = 'lic'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'lic'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        #$OBJECT BAND $passband
        self.paramlist.append(['filetype', 'png', 'Type of file to output. ' +\
                               'Supported formats: eps, pdf, pgf, png, ' +\
                               'ps, raw, rgba, svg, svgz.'])
        self.paramlist.append(['calc_lic', True,  'Calculate the LIC? The ' +\
                               'LIC must have been previous calculated and ' +\
                               'the fits file that was outputted must be ' +\
                               'used as datain. Note that the bgimage, ' +\
                               'rebin, polepol, croppolepol, and showbin ' +\
                               'will be the parameters that were used to ' +\
                               'create the original. (T/F)'])
        self.paramlist.append(['show_beam', True,  'Show beam size? (T/F)'])
        self.paramlist.append(['bgimage', '', 
                               'calc_lic must be true. Supply a background ' +\
                               'image. Use an empty string to make an ' +\
                               'intensity map with the HAWC data.'])
        self.paramlist.append(['rebin', 15.0, 'Rebinning factor for the ' +\
                               'final image'])
        self.paramlist.append(['colormap', 'binary_r', 'Colormap to be used.' +\
                               '(i.e. \'planck\', \'viridis\', \'inferno\')'])    
        self.paramlist.append(['rotated', True, 'Use 90-deg rotated pol '+\
                               'data? (T/F)'])
        self.paramlist.append(['croppolepol', True, 'Mask out map areas ' +\
                               'with P/sigP < polepol? (T/F)'])
        self.paramlist.append(['showmap', False, 'Show color map in masked' +\
                               'out areas with P/sigP < polepol? (T/F)'])
        self.paramlist.append(['polepol', 3.0, 'Value of polepol to be' +\
                               'used for the masking. pol > (polepol)epol. ' +\
                               'NOTE THAT THE POLEPOL IN THE SUBTITLE IS '+\
                               'THE ONE USED TO PLOT POLARIZATION VECTORS,'+\
                               '(if ploting pol vectors, otherwise, ignore '+\
                               'this) WHICH COMES FROM CALCULATIONS DONE IN '+\
                               'STEPREGIONS. IF PUBLISHING, ENSURE THESE '
                               'NUMBERS ARE THE SAME.'])
        self.paramlist.append(['centercrop', [], 
                               'Center and crop at specific coordinates, ' +\
                               'e.g. [275.097,-16.2105,0.0385,0.040] or ' +\
                               'an empty list [] to not center crop. ' +\
                               'values are RA center (deg), DEC center '+\
                               '(deg), Width (deg), Height (deg)'])
        self.paramlist.append(['lowhighscale', 'automatic', 'Contrast ' +\
                               'adjustment. Use [low,high], or /"automatic/"'])
        self.paramlist.append(['alphalic', 0.6, 'Transparency level for ' +\
                               'the lic overlay (0.5-0.7 generally).'])
        self.paramlist.append(['add_contours', False, 'Add contours to the ' +\
                               'map? (T/F)'])
        self.paramlist.append(['add_streamplots', False, 'Add streamplots? '+\
                               '(T/F)'])
        self.paramlist.append(['hide_axis', False, 'Hide graph axis? (T/F)'])
        ###POLMAP PARAMS###
        self.paramlist.append(['polvec', True, 'Plot polarization vectors? ' +\
                                '(T/F). If F, remaining ' +\
                                'parameters irrelevant.'])
        self.paramlist.append(['scalevec', 0.00011, 'Scale factor for ' +\
                               'vector sizes.'])
        self.paramlist.append(['scale',True,'Set to False to make ' +\
                               'all vectors the same length'])
        self.paramlist.append(['debias', True, 'Use debiased polarizations'])
        self.paramlist.append(['colorvec','red', 'Pol Vector color'])
        self.paramlist.append(['less_vec', 1, 'Removes 1/n vectors. ' +\
                               'To display all vectors, let n = 1.'])
        #END POLMAP PARAMS###
        self.paramlist.append(['filenameout', '', 'Filename'+\
                               ' for the licmap. Leave as an empty string' +\
                               'to automatically generate a filename.'])
        self.paramlist.append(['title', '', 'Title, pulled from the fits ' +\
                               'header if given an empty string.'])
        self.paramlist.append(['subtitle', '', 'Subtitle, pulled from the ' +\
                               'fits header if given an empty string.'])
        #self.paramlist.append(['streamplot_color', 'w', 'Color of the ' +\
                               #'plots. \'w\', \'b\', etc.'])
        
    def make_lic(self):
        """
        Create and save the LIC overlay.
        """         
        
        #Checking to make sure the input data doesn't already have an LIC saved
        #to it.
        try:
            I = self.dataout.imageget("LIC Background Image")
        except ValueError:
            I = ''
            
        if type(I) != str:
            message = "The input data already has an LIC saved to it. " +\
            "Supply the orginal file or set calc_lic to false."
            self.log.error(message)
            raise RuntimeError(message)

        #Import polarization info
        if self.getarg('rotated'):
            ang = self.datain.imageget('ROTATED POL ANGLE') 
        else:
            ang = self.datain.imageget('POL ANGLE')
        pol = self.datain.imageget('PERCENT POL')
        epol = self.datain.imageget('ERROR PERCENT POL')
 
        # Find NaNs inside the map and assign them to zeros 
        #(ok for the pupose of drawing the LIC) 
        mask = np.where(np.isnan(pol))
        ang[mask] = 0. 
        pol[mask] = 0. 
        epol[mask] = 9999999.
                       
        #use the last header as an example header, so we can reset some values.
        exheader = self.datain.imgheads[-1]
        
        #Recalculate WCS to account for the rebinning.
        naxis1 = self.datain.header['NAXIS1']
        naxis2 = self.datain.header['NAXIS2']

        exheader['NAXIS1'] = int(self.getarg('rebin') \
                           * naxis1)
        exheader['NAXIS2'] = int(self.getarg('rebin') \
                           * naxis2)
        exheader['CDELT1'] /= self.getarg('rebin')
        exheader['CDELT2'] /= self.getarg('rebin')
        exheader['CRPIX1'] = ((self.datain.header['CRPIX1'] / 
                           naxis1) * exheader['NAXIS1'])
        exheader['CRPIX2'] =  ((self.datain.header['CRPIX2'] / 
                           naxis2) * exheader['NAXIS2'])
        exheader['BMAJ'] =  self.datain.header['BMAJ']
        exheader['BMIN'] =  self.datain.header['BMIN']
        
        self.exheader = exheader

        #Get background image, smooth data, reprojecting image if necessary.
        if self.getarg('bgimage') == '':
            I = self.dataout.imageget('Stokes I')
            I = sp.ndimage.interpolation.zoom(I, self.getarg('rebin'), 
                                              prefilter=False)  
        else:
            ''' If you are running the pipeline on your Mac PC, reproject will
            likely give you an error asking for you to give a path to the 
            C complier. The suggested command won't work, use
            sudo ln -s /usr/bin/clang /usr/bin/gcc-4.2 '''
            from reproject import reproject_interp
            bgfile = fits.open(self.getarg('bgimage'))
            bgfile_primary = bgfile[0]

            reprojected_bg, reprojected_footprint = reproject_interp(
                    bgfile_primary, self.exheader)
            
            #this sets the values that are not behind the lic to 0
            #reprojected_bg[mask] = 0.
            I = reprojected_bg

        # Show only values above polepol.
        if self.getarg('croppolepol'):
            ''' First we make a liced without the masked out components,
            this way we can make a mask later that allows us to set only
            the parts outside of the map to be white and the inside masked 
            out parts to be grey. '''
            # Px/Py Components
            px = -1.0*pol*np.sin(np.radians(ang))
            py = pol*np.cos(np.radians(ang))

            self.px = sp.ndimage.interpolation.zoom(px, self.getarg('rebin'), 
                                               prefilter=False)
            self.py = sp.ndimage.interpolation.zoom(py, self.getarg('rebin'), 
                                               prefilter=False)

            #Create noise image as base for convolution.
            speckle = np.random.rand(self.px.shape[0], self.px.shape[1])

            #Use the LIC module to smooth
            liced0 = lic.lic(self.px, self.py, speckle, niter=2, length=int(
                             0.10*self.px.shape[1]), normalize=True, 
                             log=self.log)
            
            #Mask out the parts that are not signficant
            mask = np.where(pol/epol < self.getarg('polepol')) 
            pol[mask] = 0.
            #Now everything that uses pol needs to be redone, 
            #so we can get an LIC that is not masked.

        # Px/Py Components
        px = -1.0*pol*np.sin(np.radians(ang)) 
        py = pol*np.cos(np.radians(ang))

        #smooth data, reprojecting image if necessary
        self.px = sp.ndimage.interpolation.zoom(px, self.getarg('rebin'), 
                                           prefilter=False)
        self.py = sp.ndimage.interpolation.zoom(py, self.getarg('rebin'), 
                                           prefilter=False)

        #Create noise image as base for convolution
        speckle = np.random.rand(self.px.shape[0], self.px.shape[1])
        #Use the LIC module to smooth
        liced = lic.lic(self.px, self.py, speckle, niter=2, length=int(0.10*
                        self.px.shape[1]), normalize=True, log=self.log)
        
        if self.getarg('croppolepol'):
            #Mask_outside has the locations of 0 values outside of the LIC.
            mask_outside = np.where(liced0 == 0.0)
      
            #Mask has the locations of all zero values, 
            #inside and outside of the LIC.
            mask = np.where(liced == 0.0)
      
            #actual_data has all the locations with LIC data.
            actual_data = np.where(liced != 0.0)
            liced[mask_outside] = 1.0  
            liced[mask] += 1.0
            mask_inside_map = np.where(liced == 1.0)
   
        #Masking the data according to the croppolepol and showmap params.
        if self.getarg('croppolepol') and self.getarg('showmap'):
            '''
            Instead of setting all 0 values to white, and hiding the color 
            map below the LIC, this code takes the parts inside the LIC 
            and sets them to the median value of the data so that it 
            blends in, but is still clearly not data.
            '''
            medianval = np.median(liced[actual_data])
      
            '''
            (1) This fixes the bit where we set values to 1 and added 1 
            further up in the code and makes the outside white. Also, if 
            you use liced instead of liced0, you get bad edges around the 
            masked bits since the LIC had no data to use around there.
            '''
            liced0[mask_inside_map] = medianval
            liced0[mask_outside] = np.nan
      
            '''
            (2) We want to use the LIC that was done without cropped data, 
            since it is smooth. We will still be masking out the low 
            signal to noise spots.
            '''
            liced = liced0
      
        elif self.getarg('croppolepol'):
            #Instead of setting to low signal to noise spots to a 
            #median value, here we white them out.
            liced0[mask_inside_map] = np.nan
            
            #See (1) in if statement.
            I[mask_inside_map] = np.nan
            liced0[mask_outside] = np.nan
            I[mask_outside] = np.nan
            
            #See (2) in if statement.
            liced = liced0
      
        else:
            #Here, we don't mask anything. We just set the 0's to np.nan, 
            #so they can plot correctly.
            mask = np.where(liced == 0.0)
            I[mask] = np.nan
            liced[mask] = np.nan
            
        #By saving this data, we can make LIC images without redoing the
        #above calculations by providing the outputted file next time.        
        self.dataout.imageset(liced, "LIC Map", exheader)
        self.dataout.imageset(I, "LIC Background Image", exheader)

    def plot_lic(self):
        """
        Plotting the lic over a supplied image, if no image is supplied, 
        over the datain intensity.
        """         
        if not self.getarg("calc_lic") and self.getarg('bgimage') != '':
            message = "bgimage is added in the calc_lic section. " +\
                     "Leave bgimage as an empty string OR set calc_lic=True."
            self.log.error(message)  
            raise RuntimeError(message)

        try:
            I = self.dataout.imageget("LIC Background Image")
        except ValueError:
            message = "You need to supply a fits image with lic " +\
                  "data saved in the hdu ( _LIC.fits should be at the end " +\
                  "of the filename) OR set calc_lic to true."
            self.log.error(mesage)
            raise RuntimeError(message)
        
        ### Set text for title and subtitle
        #Read data cuts from stepregion
        minI = self.dataout.getheadval("CUTMINI")
        minP = self.dataout.getheadval("CUTMINP")
        sigma = self.dataout.getheadval("CUTPSIGP")
        minisigi = self.dataout.getheadval("CUTISIGI")
        maxp = self.dataout.getheadval("CUTMAXP")

        obj = self.dataout.getheadval('OBJECT')
        band = self.dataout.getheadval('SPECTEL1')
        nvec = self.dataout.getheadval("NVECCUT")

        #obj = obj.replace('_','\_')
        #band = band.replace('_','\_')
        
        if self.getarg('rotated'):
            eorb = "B"
        else:
            eorb = "E"
            
        #Add polarization vectors
        #The following code comes mostly from steppolmap
        if self.getarg('polvec'):
            nhwp = self.dataout.getheadval('nhwp')
    
            # Check if any pol. vect. was found after data cuts
            try:
                poldata = self.dataout.tableget('FINAL POL DATA')
            except ValueError:
                poldata = None
            #if nhwp == 1 or poldata is None:
            if nhwp == 1:
                self.log.info('No polarization data, so skipping plot ' +\
                              'polarization vectors.')
            else:
                # Check if any pol. vect. was found after data cuts
                poldataexist = True
                try:
                    poldata = self.dataout.tableget('FINAL POL DATA')
                except ValueError:
                    poldataexist = False
                    message = "No polarization data passed the elimation " +\
                             "critera , so no LIC can be made."
                    self.log.error(message)  
                    raise RuntimeError(message)
                
            # Title
            if self.getarg('title') == '':
                titlestr = "Object: %s, Band: %s, Polarization %s vectors "
                title = titlestr%(obj,band[-1],eorb)
            else:
                title = self.getarg('title')
         
            #Subtitle 
            if self.getarg('subtitle') == '':
                if self.getarg("calc_lic"):
                    fn = self.datain.filename.replace('.fits', 
                                                      '_LIC.fits')
                    fname = os.path.basename(fn)
                else:
                    fname = os.path.basename(self.datain.filename)
                subtitle = "Filename: %s"%fname
                if self.getarg('bgimage') != '':
                    fname2 = os.path.basename(self.getarg('bgimage'))
                    subtitle2 = "Background Image: %s"%fname2
                            
                #Subtitle 
                if self.getarg('subtitle') == '':
                    if nvec > 0.5:
                        sub2 ="Pol. data selection (vectors): $p/\sigma p >$"+\
                        " %s ; %s $< p (%s) <$ %s ; $I/peak(I) >$ %s ; " +\
                        "$I/\sigma I >$ %s ; N. vectors = %s"
                        subtitle2 = sub2%(str(sigma), str(minP), str("\%"),
                                          str(maxp), str(minI), str(minisigi),
                                          str(nvec))
                    else:
                        sub2 = "Pol. data selection: $p/\sigma p >$ %s ; " +\
                        "%s $< p (%s) <$ %s ; $I/peak(I) >$ %s ; " +\
                        "$I/\sigma I >$ %s ; No vectors found after cuts"
                        subtitle2 = sub2%(str(sigma), str(minP), str("\%"), 
                                          str(maxp), str(minI), str(minisigi))
                    subtitle = subtitle + ' \n ' + subtitle2
                else:
                    subtitle = self.getarg('subtitle')
    
                # Set parameters for Pol vectors
                cdelt = abs(float(self.dataout.getheadval('CDELT1')))
                if poldataexist:
                    ra = poldata['Right Ascension']
                    dec = poldata['Declination']
                    theta = poldata['Theta']
                    if self.getarg('debias'):
                        pol = poldata['debiased percent pol']
                    else:
                        pol = poldata['percent pol']
                    if self.getarg('rotated'):
                        theta = poldata['rotated theta']
                    else:
                        theta = poldata['theta']
                    if self.getarg('scale') is False: 
                        polplot = self.getarg('scalevec')*np.ones(
                                theta.shape[0])*5.0
                    else: 
                        polplot = self.getarg('scalevec')*pol
    
                #Sets up where pol vectors go.
                if self.getarg('bgimage') == '':
                    data = fits.PrimaryHDU(data=self.dataout.imageget(
                            'Stokes I'), header=self.exheader)
                    fig = aplpy.FITSFigure(data)
                else:
                    fig = aplpy.FITSFigure(self.getarg('bgimage'),
                                           hdu='Stokes I')
                    exdata = DataParent(config=self.config).load(
                            self.getarg('bgimage'))
                    data = fits.PrimaryHDU(data=exdata.imageget('Stokes I'),
                                           header=exdata.header)

                if self.getarg('lowhighscale') == 'automatic':
                    if self.getarg('bgimage') == '':
                        image = self.dataout.imageget('Stokes I')
                        med = np.nanmedian(image)
                        sd = np.nanstd(image)
                        lowscale = med - 1.*sd
                        highscale = med + 3.*sd
                    else:
                        image = exdata.imageget('Stokes I')
                        med = np.nanmedian(image)
                        sd = np.nanstd(image)
                        lowscale = med - 1.*sd
                        highscale = med + 3.*sd
                else:
                    lowscale,highscale = self.getarg('lowhighscale')
                
                # Plot vectors
                if poldataexist: 
                    line_list=[]
                    smallp_listra = [] 
                    smallp_listdec = [] 
                    smallp_radius = [] 
                    for i in range(0,len(ra)):
                        if pol[i] >= 0.0: #temporarily disabling small dot plot
                            ra1 = ra[i]-0.5*polplot[i]*np.sin(theta[i]*np.pi/
                                    180.)/np.cos(dec[i]*np.pi/180.)
                            dec1 = dec[i]-0.5*polplot[i]*np.cos(
                                    theta[i]*np.pi/180.)
                            ra2 = ra[i]+0.5*polplot[i]*np.sin(theta[i]*
                                    np.pi/180.)/np.cos(dec[i]*np.pi/180.)
                            dec2 = dec[i]+0.5*polplot[i]*np.cos(
                                    theta[i]*np.pi/180.)
                     
                            tempwcs = astwcs.WCS(self.exheader)
                            rap1, decp1 = tempwcs.wcs_world2pix(ra1,dec1,1)
                            rap2, decp2 = tempwcs.wcs_world2pix(ra2,dec2,1)
                            
                            line = np.array([[rap1, rap2], [decp1, decp2]])
                            line_list.append(line)
                        else:
                            smallp_listra.append(float(ra[i]))
                            smallp_listdec.append(float(dec[i]))
                            smallp_radius.append(0.2*cdelt)
        
                    if len(smallp_listra) > 0:
                        fig.show_circles(smallp_listra, smallp_listdec,
                                         smallp_radius, 
                                         color=self.getarg('colorvec'),
                                         facecolor = self.getarg('colorvec'))
        else:
            data = fits.PrimaryHDU(data=self.dataout.imageget(
                    'Stokes I'), header=self.exheader)
            fig = aplpy.FITSFigure(data)

            #Title
            if self.getarg('title') == '':
                titlestr = "Object: %s, Band: %s, Polarization %s"
                title = titlestr%(obj,band[-1],eorb)
            else:
                title = self.getarg('title')
     
            #Subtitle 
            if self.getarg('subtitle') != '':
                sub2 = "Pol. data selection: $p/\sigma p >$ %s ; " +\
                "%s $< p (%s) <$ %s ; $I/peak(I) >$ %s ; $I/\sigma I >$ %s"
                subtitle2 = sub2%(str(self.getarg('polepol')), str(minP), str("\%"), 
                                  str(maxp), str(minI), str(minisigi))
                subtitle = self.getarg('subtitle') + ' \n ' + subtitle2
            else:
                subtitle = self.getarg('subtitle')

        #Axes
        if self.getarg('hide_axis'):
            plt.tick_params(axis='x',which='both',bottom='off',
                            top='off',labelbottom='off')
            plt.tick_params(axis='y',which='both',left='off',
                            right='off',labelleft='off')
        else:   
            ax = plt.axes(projection=astwcs.WCS(self.exheader))
            ra = ax.coords[0]
            dec = ax.coords[1]
            ra.set_major_formatter('hh:mm:ss.s')
            dec.set_major_formatter('dd:mm:ss.s')
            plt.xlabel('RA (J2000)')
            plt.ylabel('Dec (J2000)')

        if self.getarg('lowhighscale') == 'automatic' and self.getarg(
                'bgimage') == '':
            std = np.arcsinh(np.nanmedian(self.dataout.imageget('ERROR I')))
            Imax = np.arcsinh(np.nanmax(self.dataout.imageget('STOKES I')))
            n_levels = 20
            levels = np.linspace(std,Imax,n_levels)
            lowscale = np.min(levels)
            highscale = np.max(levels)*(0.8)
        #For non-hawc images, this will probably need to be adjusted
        elif self.getarg('lowhighscale') == 'automatic' and self.getarg(
                'bgimage') != '':
            near = np.percentile(np.arcsinh(I),5)
            far = np.percentile(np.arcsinh(I),95)
            n_levels = 20
            levels = np.linspace(near, far,n_levels)
        else:
            lowscale, highscale = self.getarg('lowhighscale')

        #Optionally, include grid in image.
        plt.grid()
    
        #This loads in the planck RGB colormap, if that is the colormap used.
        if self.getarg('colormap') == 'planck':
            from matplotlib.colors import ListedColormap
            from io import StringIO   # StringIO behaves like a file object
            planckRGB = StringIO(u"  0   0 255\n  0   2 255\n  0   5 255\n  "
                                 + "0   8 255\n  0  10 255\n  0  13 255\n  "
                                 + "0  16 255\n  0  18 255\n  "
                                 + "0  21 255\n  0  24 255\n  0  26 255\n  "
                                 + "0  29 255\n  0  32 255\n  0  34 255\n  "
                                 + "0  37 255\n  0  40 255\n  0  42 255\n  "
                                 + "0  45 255\n  0  48 255\n  0  50 255\n  "
                                 + "0  53 255\n  0  56 255\n  0  58 255\n  "
                                 + "0  61 255\n  0  64 255\n  0  66 255\n  "
                                 + "0  69 255\n  0  72 255\n  0  74 255\n  "
                                 + "0  77 255\n  0  80 255\n  0  82 255\n  "
                                 + "0  85 255\n  0  88 255\n  0  90 255\n  "
                                 + "0  93 255\n  0  96 255\n  0  98 255\n  "
                                 + "0 101 255\n  0 104 255\n  0 106 255\n  "
                                 + "0 109 255\n  0 112 255\n  0 114 255\n  "
                                 + "0 117 255\n  0 119 255\n  0 122 255\n  "
                                 + "0 124 255\n  0 127 255\n  0 129 255\n  "
                                 + "0 132 255\n  0 134 255\n  0 137 255\n  "
                                 + "0 139 255\n  0 142 255\n  0 144 255\n  "
                                 + "0 147 255\n  0 150 255\n  0 152 255\n  "
                                 + "0 155 255\n  0 157 255\n  0 160 255\n  "
                                 + "0 162 255\n  0 165 255\n  0 167 255\n  "
                                 + "0 170 255\n  0 172 255\n  0 175 255\n  "
                                 + "0 177 255\n  0 180 255\n  0 182 255\n  "
                                 + "0 185 255\n  0 188 255\n  0 190 255\n  "
                                 + "0 193 255\n  0 195 255\n  0 198 255\n  "
                                 + "0 200 255\n  0 203 255\n  0 205 255\n  "
                                 + "0 208 255\n  0 210 255\n  0 213 255\n  "
                                 + "0 215 255\n  0 218 255\n  0 221 255\n"
                                 + "  6 221 254\n 12 221 253\n 18 222 252\n"
                                 + " 24 222 251\n 30 222 250\n 36 223 249\n"
                                 + " 42 223 248\n 48 224 247\n 54 224 246\n"
                                 + " 60 224 245\n 66 225 245\n 72 225 244\n"
                                 + " 78 225 243\n 85 226 242\n 91 226 241\n"
                                 + " 97 227 240\n103 227 239\n109 227 238\n"
                                 + "115 228 237\n121 228 236\n127 229 236\n"
                                 + "133 229 235\n139 229 234\n145 230 233\n"
                                 + "151 230 232\n157 230 231\n163 231 230\n"
                                 + "170 231 229\n176 232 228\n182 232 227\n"
                                 + "188 232 226\n194 233 226\n200 233 225\n"
                                 + "206 233 224\n212 234 223\n218 234 222\n"
                                 + "224 235 221\n230 235 220\n236 235 219\n"
                                 + "242 236 218\n248 236 217\n255 237 217\n"
                                 + "255 235 211\n255 234 206\n255 233 201\n"
                                 + "255 231 196\n255 230 191\n255 229 186\n"
                                 + "255 227 181\n255 226 176\n255 225 171\n"
                                 + "255 223 166\n255 222 161\n255 221 156\n"
                                 + "255 219 151\n255 218 146\n255 217 141\n"
                                 + "255 215 136\n255 214 131\n255 213 126\n"
                                 + "255 211 121\n255 210 116\n255 209 111\n"
                                 + "255 207 105\n255 206 100\n255 205  95\n"
                                 + "255 203  90\n255 202  85\n255 201  80\n"
                                 + "255 199  75\n255 198  70\n255 197  65\n"
                                 + "255 195  60\n255 194  55\n255 193  50\n"
                                 + "255 191  45\n255 190  40\n255 189  35\n"
                                 + "255 187  30\n255 186  25\n255 185  20\n"
                                 + "255 183  15\n255 182  10\n255 181   5\n"
                                 + "255 180   0\n255 177   0\n255 175   0\n"
                                 + "255 172   0\n255 170   0\n255 167   0\n"
                                 + "255 165   0\n255 162   0\n255 160   0\n"
                                 + "255 157   0\n255 155   0\n255 152   0\n"
                                 + "255 150   0\n255 147   0\n255 145   0\n"
                                 + "255 142   0\n255 140   0\n255 137   0\n"
                                 + "255 135   0\n255 132   0\n255 130   0\n"
                                 + "255 127   0\n255 125   0\n255 122   0\n"
                                 + "255 120   0\n255 117   0\n255 115   0\n"
                                 + "255 112   0\n255 110   0\n255 107   0\n"
                                 + "255 105   0\n255 102   0\n255 100   0\n"
                                 + "255  97   0\n255  95   0\n255  92   0\n"
                                 + "255  90   0\n255  87   0\n255  85   0\n"
                                 + "255  82   0\n255  80   0\n255  77   0\n"
                                 + "255  75   0\n251  73   0\n247  71   0\n"
                                 + "244  69   0\n240  68   0\n236  66   0\n"
                                 + "233  64   0\n229  62   0\n226  61   0\n"
                                 + "222  59   0\n218  57   0\n215  55   02\n"
                                 + "211  54   0\n208  52   0\n204  50   0\n"
                                 + "200  48   0\n197  47   0\n193  45   0\n"
                                 + "190  43   0\n186  41   0\n182  40   0\n"
                                 + "179  38   0\n175  36   0\n172  34   0\n"
                                 + "168  33   0\n164  31   0\n161  29   0\n"
                                 + "157  27   0\n154  26   0\n150  24   0\n"
                                 + "146  22   0\n143  20   0\n139  19   0\n"
                                 + "136  17   0\n132  15   0\n128  13   0\n"
                                 + "125  12   0\n121  10   0\n118   8   0\n"
                                 + "114   6   0\n110   5   0\n107   3   0\n"
                                 + "103   1   0\n100   0   0")
            planck_map = ListedColormap(np.loadtxt(planckRGB)/255.)
            planck_map.set_bad("white")
            self.paramlist[5][1] = planck_map

        plt.imshow(np.arcsinh(I), interpolation='bilinear', 
                   cmap=self.getarg('colormap'), origin='lower')
                   #vmin=8.25, vmax=9.6)

        if self.getarg('add_contours'):
            #This currently only works for HAWC images when lowhighscale is 
            #automatic. You need to make levels for other images.
            plt.contour(np.arcsinh(I), levels=levels, colors='white', 
                        origin='lower',linewidths = 0.2, fmt='--')

        #plot lic as alpha map
        LIC = self.dataout.imageget("LIC Map")/np.nanmax(
                self.dataout.imageget("LIC Map"))

        plt.imshow(LIC, interpolation='catrom', cmap='binary', origin='lower', 
                   alpha = self.getarg('alphalic'))

        plt.title(title + " \n" + subtitle)
        
        if self.getarg("polvec") and poldataexist:
            #Plotting the polarization vectors and 
            #removing one out of every less_vec (a supplied #) vectors
            z=0
            less_vec = self.getarg('less_vec')
            for line in line_list:
                z += 1
                if z % less_vec == 0:
                    #line[0] is [x cord1, xcord 2] in pixels
                    plt.plot(line[0], line[1], c=self.getarg('colorvec'))

        if self.getarg('centercrop') != []:
            tempwcs = astwcs.WCS(self.exheader)
            xpix,ypix = tempwcs.wcs_world2pix(self.getarg('centercrop')[0],
                                              self.getarg('centercrop')[1],1)
            sizex = (self.getarg('centercrop')[2]/self.exheader['CDELT1'])/2.  
            plt.xlim(xpix + sizex,xpix - sizex)
            sizey = (self.getarg('centercrop')[3]/self.exheader['CDELT1'])/2.  
            plt.ylim(ypix + sizey,ypix - sizey)
            ax.relim() 
            
            if self.getarg("polvec"):
                #plotting the scalebar for polaization vectors
                #scalevec is in degrees, we need it in pixels
                #so first we find pixels per degree in x
                ppdx = sizex/self.getarg('centercrop')[2]
                sc = ppdx * self.getarg('scalevec') * 5.0
                
                #x and y are in pixels
                x = xpix + sizex + abs(sizex*2)*.92
                y = ypix + sizey + abs(sizey*2)*0.98
                plt.plot([x, x+abs(sc)], [y, y], c='black')
            #calculating beam if centercroping
            if self.getarg("show_beam"):
                """
                we want degrees/pixels, and can't get it from the fits file,
                since we are cropping the image. The 3rd and 4th centercrop
                items give the degrees on the x and y axis respectively. sizex
                and sizey are half the pixels on the x and y axis respectively.
                """
                dppx = self.getarg('centercrop')[2]/(2*sizex)
                dppy = self.getarg('centercrop')[3]/(2*sizey)
                
                #BMAJ and BMIN are in degrees, dpp is in degrees/pixels so
                #x and y are in pixels. They give the x and y dims of the beam.
                x = (self.dataout.header['BMAJ']/2)/dppx
                y = (self.dataout.header['BMIN']/2)/dppy
                
                radius = (abs(x)+abs(y))/2
                
                xlen = abs(2*sizex)
                ylen= abs(2*sizey)
                #r is a number between 0 and 1, representing 
                #(beam radius in pixels)/(total pixels being plotted)
                r = radius/((xlen+ylen)/2)
                
                #plotting here in absolute coordinates, from 0 to 1.
                ax.add_artist(Ellipse((r+0.02,r+0.02), (x/xlen)*2, (y/ylen)*2, transform = ax.transAxes, 
                                      facecolor ='r', edgecolor = 'black'))
                plt.text((r/2)+0.01, 0.01, "Beam FWHM", transform = ax.transAxes, fontsize=5)
            
        #calculating beam if not centercropping                  
        elif self.getarg("show_beam"):
            #BMAJ and BMIN are in degrees, CDELT1 is in degrees/pixels
            #so x and y are in pixels
            x = (self.dataout.header['BMAJ']/2)/self.exheader['CDELT1']
            y = (self.dataout.header['BMIN']/2)/self.exheader['CDELT2']
            
            radius = (abs(x)+abs(y))/2
            #NAXIS1 is the number of pixels on the x-axis
            xlen = self.exheader['NAXIS1']
            ylen = self.exheader['NAXIS2']

            avgpix = (xlen+ylen)/2
            #r is a number between 0 and 1, representing 
            #(beam radius in pixels)/(total pixels being plotted)
            r = radius/avgpix

            beam = plt.Circle((r+0.02,r+0.02), r, transform = ax.transAxes,
                  facecolor ='r', edgecolor = 'black')
            ax.add_artist(beam)

            plt.text((r/2)+0.01, 0.01, "Beam FWHM", 
                     transform = ax.transAxes, fontsize=5)
        
        #plots text under pol vector scale
        plt.text(0.9, 0.96, "p = 5%", transform = ax.transAxes, fontsize=8)
        if self.getarg("polvec") and self.getarg('centercrop') == []:
            #plotting the scalebar for polaization vectors
            #scalevec is in degrees, we need it in pixels
            #cdelt is degrees/pixel, so scalevec/cdelt gives scalevec in pixels
            svp = (self.getarg('scalevec')/self.exheader['CDELT1']) * 5.0
            
            #x and y are in pixels
            x = self.exheader['NAXIS1']*.945
            y = self.exheader['NAXIS2']*.98
            plt.plot([x, abs(x+svp)], [y, y], c='black')
        
        if self.getarg('add_streamplots'):
            #calculate px and py (from calc_lic)
            if not self.getarg('calc_lic'):
                #Import polarization info
                if self.getarg('rotated'):
                    ang = self.datain.imageget('ROTATED POL ANGLE')
                else:
                    ang = self.datain.imageget('POL ANGLE')
                pol = self.datain.imageget('PERCENT POL')
                
                naxisold = self.datain.header['NAXIS1']
                naxisnew = self.exheader['NAXIS1']
                rebin = naxisnew/naxisold
                
                # Find NaNs inside the map and assign them to zeros 
                #(ok for the pupose of drawing the LIC) 
                mask = np.where(np.isnan(pol))
                ang[mask] = 0.
        
                # Px/Py Components
                px = -1.0*pol*np.sin(np.radians(ang)) 
                py = pol*np.cos(np.radians(ang))
           
                #smooth data, reprojecting image if necessary
                self.px = sp.ndimage.interpolation.zoom(px, rebin, prefilter=False)
                self.py = sp.ndimage.interpolation.zoom(py, rebin, prefilter=False)
                
            grid = np.indices(self.px.shape)
            #if you want to use arrows, set arrowsize to a larger number
            plt.streamplot(grid[1],grid[0],self.px,self.py,  
                           color = self.getarg('streamplot_color'), 
                           arrowsize=0.0001)
        
        if self.getarg("calc_lic"):
            if self.getarg("filenameout") == '':
                name = self.datain.filename.replace('.fits', '_licmap')
            else:
                num = len(self.datain.filename)-len(os.path.basename(
                        self.datain.filename))
                folder = self.datain.filename[0:num]
                name = folder+self.getarg("filenameout"
                                          )+"."+self.getarg("filetype")
        else:
            if self.getarg("filenameout") == '':
                length = len(self.datain.filename)
                name = self.datain.filename[0:(length-8)]+(
                        'licmap.'+self.getarg("filetype"))
            else:
                num = len(self.datain.filename)-len(os.path.basename(
                        self.datain.filename))
                folder = self.datain.filename[0:num]
                name = folder+self.getarg("filenameout"
                                          )+"."+self.getarg("filetype")
        
        plt.savefig(name,dpi=300)
        #fig.set_system_latex(False) 
        self.log.info("Saved "+self.getarg("filetype")+" LIC map under %s"%name)
        
    def run(self):
        """
        Runs the data reduction algorithm. The self.datain is run
        through the code, the result is in self.dataout.
        """
        # Initialize self.dataout
        self.dataout = self.datain.copy()
        # Run make_lic if needed
        if self.getarg("calc_lic"):
            self.make_lic()
        else:
            self.exheader = self.datain.getheader('LIC Map')
            
        # Run plot_lic
        self.plot_lic()
        
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
            python stepregion.py input.fits -arg1 -arg2 . . .
            Standard arguments:
        --config=ConfigFilePathName.txt : name of the configuration file
        -t, --test : runs the functionality test i.e. pipestep.test()
        --loglevel=LEVEL : configures the logging output for a particular level
        -h, --help : Returns a list of
    """
    StepLIC().execute()
