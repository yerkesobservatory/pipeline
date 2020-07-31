import os
import math
import sep
import utilities as ut
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
import scipy.optimize as opt
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from astropy.utils.data import download_file
from datetime import datetime
import astroalign as aa
from scipy.integrate import dblquad

import warnings
warnings.filterwarnings("ignore")

# Preliminary Integrated fleXible Astronomy Reduction Dynamics (PIXARD)
# 5/3/20 -Siegel, J

class stack:
    def __init__(self, notes = None):
        self.comps = []
        self.notes = notes
        self.summary = None

    def summarize(self, typ = 'sum'):
        data = [s.primary for s in self.comps]
        if typ == "sum":
            data = sum(data)
        elif typ == "mean":
            data = np.mean(data, axis = 0)
        elif typ == 'median':
            data = np.median(data, axis = 0)

        self.summary = unit(data = data, header = typ, notes = self.notes)

    def smooth_to_psf(self):
        psfs = [o.source_table[o.source_table['good']]['fwhm'].median() for o in self.comps]
        mpsf = max(psfs)
        print("Observed PSFs: {}".format(psfs))
        print("Fitting to: {}".format(mpsf))
        smoothed = []
        for c in self.comps:
            c.smooth(sigma = mpsf)
            smoothed.append(c)
        self.comps = smoothed

    def get_time_deltas(self, t0 = None):
        times = np.array([t.time for t in self.comps])
        if not t0:
            t0 = times.min()
        deltas = times - t0
        self.time_deltas = deltas

    def mean_time(self, tmean = None):
        if not tmean:
            times = np.array([t.time for t in self.comps])
            self.time = times.mean()
        else:
            self.time = tmean

    def subtract(self, basefile=None):
        if basefile:
            baseunit = [f for f in self.comps if f.notes == basefile]
            if len(baseunit) != 1:
                print ("FAULTY BASEFILE. Will use earliest observation instead.")
                basefile = None
        if not basefile:
            t0 = np.array([t.time for t in self.comps]).min()
            baseunit = [f for f in self.comps if f.time == t0][0]
        
        others = [f for f in self.comps if f.notes != baseunit.notes]
        self.subtracted = []
        for obs in others:
            s =  unit(data = obs.primary - baseunit.primary, header = 'SUBTRACTED', notes = baseunit.notes +"-"+obs.notes)
            s.time = obs.time - baseunit.time
            self.subtracted.append(s)

    def split_stack(self, d = 86400):
        times = np.array([t.time for t in self.comps])
        t0 = times.min()
        storage = {}
        for obs in self.comps:
            delta = obs.time - t0
            storage[obs.notes] = delta
        sortd = []
        stacks = []
        n = 1.
        while len(sortd) < len(times):
            st = stack(notes = "")
            for key, t in storage.items():
                obs = [f for f in self.comps if f.notes == key][0]
                if t < d * n and key not in sortd:
                    sortd.append(key)
                    st.comps.append(obs)
                    st.notes += obs.notes + '-'
            if len(st.comps) > 0:
                stacks.append(st)
            n += 1

        return stacks

    def order_by_time(self):
        storage = {}
        for obs in self.comps:
            if obs.time in storage.keys():
                storage[obs.time].append(obs)
            else:    
                storage[obs.time] = [obs]
        sortd = []
        for t in sorted(storage.keys()):
            sortd += storage[t]
        self.comps = sortd.copy()
    
    def align(self, base = None, manualbase = False):
        if not manualbase:
            base = self.comps[0].primary
        for i in range(len(self.comps)):
            new, footprint = aa.register(self.comps[i].primary, base)
            self.comps[i].primary = new.copy()

class unit:
    def __init__(self, path = None, data = None, header = None, notes = None):
        self.path = path
        self.notes = notes

        if path:
            self.raw, self.header = self.openfits(path)
            self.primary = self.raw.copy()
            self.center = (self.raw.shape[0] /2, self.raw.shape[1] /2)
            self.primary_name = "raw"
        else:
            try:
                self.raw = data.byteswap(inplace=True).newbyteorder()
                self.header = header
                self.center = (self.raw.shape[0] /2, self.raw.shape[1] /2)
                self.primary = self.raw.copy()
                self.primary_name = "raw"              
            except:
                self.raw, self.header = None, None
                self.center = None
                self.primary = None
                self.primary_name = None         

        self.bkg = None
        self.bkg_sub = None
        self.sources = None
        self.source_table = None

    def openfits (self, path, changebyte = False):
        if path[0:4] == "http":
            image_file = download_file(path)
            image_data = fits.open(image_file)[0]
            data, header = image_data.data, image_data.header
        else:
            image_data = fits.open(path)[0]
            data, header = image_data.data, image_data.header

            if changebyte:
                data = data.byteswap().newbyteorder()
            data = data.astype(float)
        return data, header

    def set_primary(self, typ):
        if typ == "raw":
            self.primary = self.raw.copy()
            self.primary_name = typ
        elif typ == "bkg_sub":
            self.primary = self.bkg_sub.copy()
            self.primary_name = typ
        elif typ == "bkg":
            self.primary = np.array(self.bkg)
            self.primary_name = typ
        else:
            print("{} is not a recognized data name.".format(typ))

    def check_primary(self):
        print ("Primary is set to: {}".format(self.primary_name))
    
    def summary_stats (self, typ = 'primary'):
        if typ == 'raw':
            return np.nanmean(self.raw), np.nanstd(self.raw)
        elif typ == 'bkg':
            return np.nanmean(np.array(self.bkg)), np.nanstd(np.array(self.bkg))
        elif typ == 'bkg_sub':
            return np.nanmean(self.bkg_sub), np.nanstd(self.bkg_sub)
        elif typ == 'primary':
            return np.nanmean(self.primary), np.nanstd(self.primary)

    def extract_bkg(self):
        self.bkg = sep.Background(self.raw)

    def subtract_bkg(self):
        self.bkg_sub = self.raw - self.bkg

    def extract_sources(self, sig = 1.5):        
        self.sources = sep.extract(self.primary, sig, err=self.bkg.globalrms)

    def build_sources_table(self):       
        df = pd.DataFrame()
        objects = self.sources
        df['x'] = objects['x']; df['y'] = objects['y']; df['a'] = objects['a']; df['b'] = objects['b']; df['theta'] = objects['theta']; df['npix'] = objects['npix']; df['cflux'] = objects['cflux']
        
        argu = (df['y'] - self.center[1]) / (df['x'] - self.center[0])
        df['phi'] = np.arctan(argu)
        
        self.source_table = df

    def filter_sources(self, maxpixs = None, minpixs = None, fluxsigmas = None, MajorMinorRatio = None, centerfrac = None, edgefrac = None):
        data = self.primary
        df = self.source_table
        
        m, s = np.mean(data), np.std(data)

        if MajorMinorRatio:
            df = df[ abs(1- df['a'] / df['b'])  < MajorMinorRatio] # Remove non-circular sources
        if minpixs:
            df = df[ df['npix'] > minpixs] # Remove smallest sources
        if maxpixs:
            df = df[ df['npix'] < maxpixs] # Remove largest sources
        if fluxsigmas:
            df = df[ df['cflux'] > fluxsigmas * s + m] # Remove darkest sources

        if centerfrac:
            df = df[ (df['x'] < data.shape[0] * centerfrac) | (df['x'] > data.shape[0] * (1 - centerfrac)) ]
            df = df[ (df['y'] < data.shape[1] * centerfrac) | (df['y'] > data.shape[1] * (1 - centerfrac)) ]

        if edgefrac:
            df = df[ (df['x'] > float(data.shape[0]) * edgefrac) & (df['x'] < float(data.shape[0]) * (1 - edgefrac )) & (df['y'] > float(data.shape[1]) * edgefrac) & (df['y'] < float(data.shape[1]) *(1 - edgefrac))] # Remove sources near the center
       
        self.source_table = df

    def display(self, sources = False, cmap = "Greys", vmin = 1, vmax = 1, figsize = (5,5), savename = None):

        data = self.primary

        m, s = np.mean(data), np.std(data)

        fig, ax = plt.subplots(figsize = figsize)

        
        vmin=m-s*vmin
        vmax = m+s*vmax
        im = ax.imshow(data, interpolation='nearest', cmap=cmap,
                    vmin=vmin, vmax=vmax, origin='lower')

        # plot an ellipse for each object
        theta = np.arange(0, 2*np.pi, 0.01)
        if sources:
            df = self.source_table
            for index, row in df.iterrows():
                x, y, a, b, theta = row['x'], row['y'], row['a'], row['b'], row['theta'] 

                e = Ellipse(xy=(x, y),
                            width=6*a,
                            height=6*b,
                            angle=theta * 180. / np.pi)
                e.set_facecolor('none')
                e.set_edgecolor('red')
                ax.add_artist(e)

        plt.colorbar(im, ax = ax, aspect=10)
        plt.show()
        if savename:
            fig.savefig(savename+".png", dpi=100, box_inches='tight')

        plt.close()

    def fit_gauss(self, DIM = 5, QUALITY = 0.3, BKG_SIG = 1, MAX_SIG = 10, NFAIL = 5, display = False, integrate = False):
        
        data = self.primary
        df = self.source_table

        m, s = self.summary_stats("bkg")

        xt = np.linspace(0, 2*DIM-1, 2*DIM)
        yt = np.linspace(0, 2*DIM-1, 2*DIM) 
        xa, ya = np.meshgrid(xt, yt)

        bkgs = []
        ints = []
        wxs = []
        wys = []
        rts = []
        good = []

        for index, row in df.iterrows():
            TEMP = data.copy()
            X, Y = np.meshgrid(np.arange(0,data.shape[1],1), np.arange(0,data.shape[1],1))
            x, y = row['x'], row['y']
            # R = np.sqrt((X - x)**2 + (Y - y)**2)
            # TEMP[ R > r] = m
            TEMP = TEMP[int(y-DIM+1):int(y+DIM+1),int(x-DIM+1):int(x+DIM+1)]

            wy, wx, rot, height, bg = ut.fitgaussian2d(TEMP,0,False)
            
            gauss = ut.gaussian2d(height, x - int(x-DIM+1), y - int(y-DIM+1), wy, wx, bg, rot)(xa, ya)
            
            if integrate:
                try:

                    fn = ut.gaussian2d(height, x - int(x-DIM+1), y - int(y-DIM+1), wy, wx, bg, rot)
                    INT = dblquad(fn, -20, 20, -20, 20)

                    fn =  ut.gaussian2d(height, x - int(x-DIM+1), y - int(y-DIM+1), wy, wx, 0.0, rot)
                    BKG_SUB_INT = dblquad(fn, -20, 20, -20, 20)
    
                    bkgs.append(BKG_SUB_INT[0])
                    ints.append(INT[0])

                except:
                    print("FAILED")
                    bkgs.append(None)
                    ints.append(None)

            quality = np.abs(gauss - TEMP) / np.abs(gauss)
            nfail = quality[ quality > QUALITY].size
            if wy < MAX_SIG and wx < MAX_SIG and nfail <= NFAIL:
                if display:
                    ut.display_fit(TEMP, gauss, TEMP-gauss, x - int(x-DIM+1), y - int(y-DIM+1), QUALITY, cmap = 'plasma')  
                good.append(True)
            else:
                good.append(False)
                # if display:
                #     ut.display_fit(TEMP, gauss, quality, x - int(x-DIM+1), y - int(y-DIM+1), QUALITY, cmap = "rainbow") 

            wxs.append(wx); wys.append(wy); rts.append(rot)

        if integrate:
            df['BKG_SUB_INT'] = np.array(bkgs)
            df['INT'] = np.array(ints)
            
        df['sigx'] = np.array(wxs)
        df['sigy'] = np.array(wys)
        df['rot'] = np.array(rts)
        df['sigma'] = ( df['sigx'] + df['sigy'] ) / 2
        df['fwhm'] = df['sigma'] * 2.355
        df['good'] = good

        self.source_table = df

    def smooth(self, sigma = 1.0):
        self.primary =  gaussian_filter(self.primary, sigma = sigma)

    def time_it(self):
        date = self.header['DATE-OBS']
        yr = date[:4]; month = date[5:7]; day = date[8:10]; hour = date[11:13]; minute = date[14:16]; second = date[17:19]
        t = datetime(int(yr), int(month), int(day), int(hour), int(minute)).timestamp()
        self.time = t

    def mask_sources(self):
        masked = self.primary.copy()
        df = self.source_table
        x = np.linspace(0, masked.shape[0]-1, masked.shape[0])
        y = np.linspace(0, masked.shape[1]-1, masked.shape[1])
        X, Y = np.meshgrid(x, y)
        for index, row in df.iterrows():
            x0, y0, a, b, theta = row['x'], row['y'], row['a'], row['b'], row['theta']
            a *= 6
            b *= 6
            terma = (X - x0) * np.cos(theta) + (Y - y0) * np.sin(theta) 
            termb = (X - x0) * np.sin(theta) - (Y - y0) * np.cos(theta) 
            mask = ( terma**2 / a **2 ) + ( termb**2 / b **2 )
            masked [mask < 1.] = np.ma.masked #np.nan
        self.primary = masked.copy()

    def integrate_sources(self):
        df = self.source_table
        x = np.linspace(0, self.primary.copy().shape[0]-1, self.primary.copy().shape[0])
        y = np.linspace(0, self.primary.copy().shape[1]-1, self.primary.copy().shape[1])
        X, Y = np.meshgrid(x, y)
        integrated = []
        area = []
        for index, row in df.iterrows():
            masked = self.primary.copy()
            x0, y0, a, b, theta = row['x'], row['y'], row['a'], row['b'], row['theta']
            a *= 6
            b *= 6
            terma = (X - x0) * np.cos(theta) + (Y - y0) * np.sin(theta) 
            termb = (X - x0) * np.sin(theta) - (Y - y0) * np.cos(theta) 
            mask = ( terma**2 / a **2 ) + ( termb**2 / b **2 )
            masked [mask > 1.] = 0.0
            integrated.append(np.sum(masked))
            
            mask[mask < 1] = 1
            mask[mask != 1] = 0
            area.append(np.sum(mask))

        df['simple_sum'] = integrated
        df['area'] = area
        self.source_table = df.copy()
