'''
LINE INTEGRAL CONVOLUTION (LIC) BASED ON CABRAL ET AL. (1993) AND IDL CODE FROM
FALCETA-GONCALVES

CODE WAS TRANSLATED FROM IDL TO PYTHON

AUTHOR: J. MICHAIL

DATE: 6/1/2017
'''

#Import modules
import numpy as np
from scipy import ndimage

#LIC function definition
def lic(vx, vy, background, length = 20, niter = 10, normalize = True,
        amplitude = False, level = True, scalar = True, log = None):
    """ LIC function: Calculates the lic map for the given vx/vy direction.
    
        vx, vy: 2D array with x/y component of the direction at each point.
        
        log: a logging.Logger (or similar) object to send messages to.
             (is ignored if None)
        
    """
    nx, ny = vx.shape[1], vx.shape[0]

    uu = np.sqrt(vx**2 + vy**2)
    ii = np.where(uu == 0.0)
    uu[ii] = 1.0

    if normalize:
        ux = vx / uu
        uy = vy / uu

    else:
        ux = vx / np.nanmax(uu)
        uy = vy / np.nanmax(uu)

    for it in range(niter):
        texture = background
        
        vv = np.zeros((ny, nx))
        pi = np.tile(np.reshape(np.arange(nx), (1, nx)), (ny, 1))
        pj = np.tile(np.reshape(np.arange(ny), (ny, 1)), (1, nx))
        mi = pi
        mj = pj

        ppi = 1. * pi
        ppj = 1. * pj
        mmi = 1. * mi
        mmj = 1. * mj

        for l in range(length):
            if log:
                log.debug("LIC: Iteration=%d L=%d" % (niter, l) )
            dpi = ndimage.map_coordinates(ux.T, (ppi, ppj), order=1)
            dpj = ndimage.map_coordinates(uy.T, (ppi, ppj), order=1)
            dmi = ndimage.map_coordinates(ux.T, (ppi, ppj), order=1)
            dmj = ndimage.map_coordinates(uy.T, (ppi, ppj), order=1)
            
            ppi += 0.25 * dpi
            ppj += 0.25 * dpj
            mmi -= 0.25 * dmi
            mmj -= 0.25 * dmj

            pi = (np.round(ppi) + nx) % nx
            pj = (np.round(ppj) + ny) % ny
            mi = (np.round(mmi) + nx) % nx
            mj = (np.round(mmj) + ny) % ny

            ppi = pi + (ppi - np.round(ppi))
            ppj = pj + (ppj - np.round(ppj))
            mmi = mi + (mmi - np.round(mmi))
            mmj = mj + (mmj - np.round(mmj))
            
            vv = vv + ndimage.map_coordinates(texture.T, (ppi, ppj), order=1) + ndimage.map_coordinates(texture.T, (mmi, mmj), order=1)

        background = 0.25 * vv / length
	
        if amplitude:
	        if type(scalar) != bool:
		        uu = scalar
	        else:
		        uu = np.ones(background.shape)
			
	        if type(level) == bool:
		        level = 0.1
		
	        uu[ii] = 0.0
	        level = max([0.0, min([level, 1.0])])
	        uu = ((1.0 - level)/np.max(uu)) * uu + (np.ones(background.shape) * level)
	        background *= uu
	
        background[ii] = 0.0	
    return background
