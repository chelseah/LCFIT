#!/usr/bin/env python

import numpy as np
import pyfits as pyf
from scipy import ndimage


def main():
    """
    Example to read the specific intensity files.
    This is shown using both nearest neighbor and trilinear interpolation.
    If passing the array to C/C++ before knowing the exact values needed,
    nearest-neighbor will probably be preferable.
    """

    #filtername = 'TYCHO_B'
    filtername = '2MASS_J'
    specfile = 'spec_intens_fullgrid.' + filtername + '.fits'

    hdu = pyf.open(specfile)[0]

    # axes are log(g), log(T), mu (fits header should be right)
    nx, ny, nz = hdu.data.shape
    print 'Specific intensity for filter ' + hdu.header['FILT']
    print 'Axes are ' + hdu.header['CTYPE1'].strip() + ', ' + \
        hdu.header['CTYPE2'].strip() + ', ' + \
        hdu.header['CTYPE3'].strip() + '.'  

    log_g = hdu.header['CRVAL1'] + hdu.header['CDELT1']*np.arange(nx)
    log_T = hdu.header['CRVAL2'] + hdu.header['CDELT2']*np.arange(ny)
    mu = hdu.header['CRVAL3'] + hdu.header['CDELT3']*np.arange(nz)

    # Example points to extract
    #log_g_list = np.asarray([3.33, 4.198, 4.23])
    #log_T_list = np.log10(np.asarray([4560, 6050, 12120]))
    #mu_list = np.asarray([0.967, 0.768, 0.442])
    log_g_list = np.array([4.436])
    log_T_list = np.array([3.762]) 
    mu_list = np.array([0.])
    # Pixel IDs
    indx_g = (log_g_list - hdu.header['CRVAL1'])/hdu.header['CDELT1']
    indx_T = (log_T_list - hdu.header['CRVAL2'])/hdu.header['CDELT2']
    indx_mu = (mu_list - hdu.header['CRVAL3'])/hdu.header['CDELT3']

    # Nearest neighbor, explicit index calculation
    indx_g_int = (indx_g + 0.5).astype(int)
    indx_T_int = (indx_T + 0.5).astype(int)
    indx_mu_int = (indx_mu + 0.5).astype(int)
    print indx_g_int, indx_T_int, indx_mu_int, hdu.header['CRVAL2'],hdu.header['CDELT2']
    print hdu.data[indx_g_int, indx_T_int, indx_mu_int]
    
    # Nearest neighbor, using map_coordinates
    print ndimage.map_coordinates(hdu.data, [indx_g, indx_T, indx_mu], order=0)

    # Trilinear
    print ndimage.map_coordinates(hdu.data, [indx_g, indx_T, indx_mu], order=1)
         


if __name__ == "__main__":
    main()
