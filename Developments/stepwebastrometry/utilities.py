import os
import math
import sep
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
import scipy.optimize as opt
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
from astropy.utils.data import download_file

def gather_files(dir):  
    return list(os.listdir(dir))

def gaussian2d(height, center_x, center_y, width_x, width_y, bgoffset, rotation):
    """ Returns a gaussian function with the given parameters.
        The following parameters are used:
    """
    width_x = float(width_x)
    width_y = float(width_y)
    rot = np.deg2rad(rotation)
    #center_x = center_x * np.cos(rotation) - center_y * np.sin(rotation)
    #center_y = center_x * np.sin(rotation) + center_y * np.cos(rotation)

    def gauss2d(x,y):
        x_off = (x - center_x) * np.cos(rot) - (y - center_y) * np.sin(rot)
        y_off = (x - center_x) * np.sin(rot) + (y - center_y) * np.cos(rot)
        #x_rot = x * np.cos(rotation) - y * np.sin(rotation)
        #y_rot = y * np.sin(rotation) + y * np.cos(rotation)
        g = height*np.exp(-((x_off/width_x)**2 + (y_off/width_y)**2)/2.) + bgoffset
        return g
    
    return gauss2d

def moments2d(data):
    """ Returns (height, x, y, width_x, width_y, bgoffset) the gaussian parameters
        of a 2D distribution by calculating its moments. This function is used
        to get an initial estimate for the parameters.
    """
    total = np.nansum(data)
    bgoffset = np.nanmedian(data)
    Y, X = np.indices(data.shape)
    ysize, xsize = data.shape
    if abs(total) == 0: total = 0.1
    x = np.nansum(X*data)/total
    y = np.nansum(Y*data)/total
    # print(" Moments: x=%f y=%f tot=%f bgoff=%f" % (x, y, total, bgoffset))
    # If the initial guess for x and y is outside the array,
    # will assume a guess in the center of the image
    # (this situation might happen in case there is no source
    # in the image. The assumption is used to avoid an index
    # problem in the col and row definitions below
    if x >= xsize or x <= 0 or np.isnan(x):
        x = xsize/2.
    if y >= ysize or y <= 0 or np.isnan(y):
        y = ysize/2.
    col = data[int(y), :]
    width_x = np.sqrt(np.nansum(abs((np.arange(col.size)-y)**2*col)) / abs(np.nansum(col)))
    row = data[:, int(x)]
    width_y = np.sqrt(np.nansum(abs((np.arange(row.size)-x)**2*row)) / abs(np.nansum(row)))
    if abs(np.nanmax(data)) > abs(np.nanmin(data)):
        height = np.nanmax(data)
    if abs(np.nanmax(data)) < abs(np.nanmin(data)):
        height = np.nanmin(data)
    # print(" Moments: returning h=%f x=%f y=%f wx=%f wy=%f bg=%f" % (height, x, y, width_x, width_y, bgoffset))
    return height, x, y, width_x, width_y, bgoffset, 0.0

def fitgaussian2d(data, nanpix, fillbadpix):
    """ Returns (height, x, y, width_x, width_y) the gaussian parameters
        of a 2D distribution found by a fit
    """
    params = moments2d(data)
    data_img = data.copy()
    # If fillbadpix is False,
    # assume that for bad pixels (NaNs), the value on the array is equal to the model
    # (so that the difference -- errorfunction -- will be zero)
    if fillbadpix == False:
        model = gaussian2d(*params)(*np.indices(data.shape))
        data_img[nanpix] = model[nanpix]
    else:
        pass
    errorfunction = lambda p: np.ravel(gaussian2d(*p)(*np.indices(data_img.shape)) - data_img)
    p, success = opt.leastsq(errorfunction, params)
    amp,centy,centx,widy,widx,bgfit,rot = p
    # Make sure rot is between -180 and 180
    while(rot<-180.): rot+=360.
    while(rot>180.): rot-=360.
    # print(" Fit Values: returning h=%f, x=%f y=%f wx=%f wy=%f bg=%f rot=%f" % (amp, centx, centy, widx, widy, bgfit, rot))
    # return amp,centy,centx,widy,widx,bgfit,rot,success
    return widy,widx, rot, amp, bgfit

def gaussian(x, y, xmu, ymu, sigma):

    X = x - xmu
    Y = y - ymu
    C =1 / (sigma**2 * 2*math.pi)
    E = np.exp( -(X**2 + Y**2) / (2 * sigma **2))
    return C * E

def display_fit(raw, gauss, qual, xc, yc, quality, cmap = 'virdis'):
    vmax = np.max ( [raw.max(), gauss.max()]) 
    vmin = np.min ( [raw.min(), gauss.min()]) 

    fig, ((ax1), (ax2), (ax3)) = plt.subplots(nrows=3, ncols=1)
    fig.set_size_inches((2160/550),(3240/400))

    im = ax1.imshow(raw, interpolation='nearest', cmap = cmap,
                origin='lower', vmin = vmin, vmax = vmax)
    ax1.scatter(xc, yc)
    fig.colorbar(im, ax=ax1, aspect=10)

    im = ax2.imshow(gauss, interpolation='nearest', cmap = cmap,
                origin='lower', vmin = vmin, vmax = vmax)
    ax2.scatter(xc, yc)
    plt.colorbar(im, ax = ax2, aspect=10)

    im = ax3.imshow(qual, interpolation='nearest', cmap = cmap, 
                origin='lower')#, vmin = 0, vmax = quality)
    ax3.scatter(xc, yc)
    plt.colorbar(im, ax = ax3, aspect=10)
    # fig.tight_layout()
    plt.show()

def norm(X):
    return (X - X.min()) / (X.max() - X.min())

def plot_values(xs, ys, xs_err = None, ys_err = None, ys_error = None, yerrors = [False], xerrors = [False], xlog = False, ylog = True, lines = [None],
               colors = ['steelblue', 'firebrick', 'k', 'steelblue', 'firebrick', 'thistle'],
               markersize = 25, linewidth = 2, legend = None, labels = None, ncol = 1,
               xlim = None, ylim = None, savename = None):
    """
    function plotting a given dataset
    
    Parameters:
    -----------------------------------------------------------
    xs, yx: lists of arrays
    xlog, ylog: booleans for the scale
    lines: list determining how to values are 
           plotted (None for points or '-', '--', ':', ':.' for lines)
    colors: list of matplotlib colors
    markersize: int
    linewidth: int
    legend: None for no legend or list of label for legend
    label: None for no axis labels or list of axis labels
    ncol: int
    
    Returns:
    -------------------------------------------------------------
    plot object
    """
    fig, ((ax1)) = plt.subplots(figsize = (10, 10), nrows=1, ncols=1)
    plt.subplots_adjust(wspace=.3)
    for i, x in enumerate(xs):
        y = ys[i]
        line = lines[i]
        color = colors[i]
        yerror = yerrors[i]
        xerror = xerrors[i]

        if legend:
            label = legend[i]
        else:
            label = None
        if line:
            plt.plot(x, y, color = color, linewidth = linewidth, linestyle = line, label = label)

        if not xerror and not yerror:
            plt.scatter(x, y, color = color, s = markersize, label = label)
        elif xerror and not yerror:
            plt.errorbar(x, y, xerr = xs_err[i], color = color, markersize = markersize, linestyle ='', marker = 'o')
        elif yerror and not xerror:
            plt.errorbar(x, y, yerr = ys_err[i], color = color, markersize = markersize, linestyle ='', marker = 'o')
        else:
            plt.errorbar(x, y, xerr = xs_err[i], yerr = ys_err[i], color = color, markersize = markersize, linestyle ='', marker = 'o')

    if ylog:
        plt.yscale('log')
    if xlog:
        plt.xscale('log')
    if labels:
        ax1.set_xlabel(labels[0], fontsize = 25)
        ax1.set_ylabel(labels[1], fontsize = 25)

    if ylim:
        ax1.set_ylim(ylim)
    if xlim:
        ax1.set_xlim(xlim)

    ax1.tick_params(labelsize = 15, which = 'major', length = 10, width = 2)
    ax1.tick_params(which = 'minor', length = 6, width = 1)
    if legend:
        plt.legend(frameon = False, fontsize = 15, ncol = ncol)
    plt.show()
    if savename:
        fig.savefig(savename+".png", dpi=100, box_inches='tight')
