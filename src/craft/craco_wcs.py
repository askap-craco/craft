#!/usr/bin/env python
"""
WCS utilities for CRACO

Borrowing some of the idea from ORD et al 
"Interferometric Imaging with the 32 Element Murchison Wide-Field Array"


https://iopscience.iop.org/article/10.1086/657160/pdf

Copyright (C) CSIRO 2020
"""
import numpy as np
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation,AltAz
from astropy import units as u
#from astropy.utils import iers
#iers.conf.auto_download = False
from numpy import sin,cos,tan,arccos,arctan


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def calc_ab(H, delta, phi):
    '''
    Calculate a, and b coefficients according to equations 4--7 from Ord et al.
    
    :param: H - target hour angle (radians)
    :param: delta - target declination (radians)
    :param: phi - observatory latitude (radians)
    :return: (a,b) coeffiicients
    '''

    # chi is parallactic angle
    tan_chi = cos(phi) * sin(H) / (sin(phi)*cos(delta) - cos(phi)*sin(delta)*cos(H)) # equation 6
    chi = arctan(tan_chi)

    # Z is zenith angle
    cos_Z = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(H) # equation 7
    Z = arccos(cos_Z)
    
    a = tan(Z) * sin(chi) # Equation 4
    b = -tan(Z) * cos(chi) # Equation 5

    return (a, b)

class CracoWCS:
    def __init__(self, target: SkyCoord, time: Time, cellsize, npix:int, inttime, site: EarthLocation=None):
        '''
        :param: target - SkyCoord look direction phase center
        :time: Time to set WCS to
        :cellsize: (lcell, mcell) cell size 2 tuple as astropy Angle
        :npix: Numebr of pixels (int)
        :inttime: Integration tiem as astropy time Quanitty
        :site: if None defaults to EarthLocation.of_site('ASKAP') else the site
        '''
        if site is None:
            site = EarthLocation.of_site('ASKAP')
        
        lst = time.sidereal_time('apparent', site.lon)
        altaz = target.transform_to(AltAz(obstime=time, location=site))
        hour_angle = lst - target.ra
        a,b = calc_ab(hour_angle.rad, target.dec.rad, site.lat.rad)

        self.site = site
        self.time = time
        self.target = target
        self.lst = lst
        self.altaz = altaz
        self.hour_angle = hour_angle
        self.a = a
        self.b = b
        self.cellsize = cellsize
        self.npix = npix

        lcell, mcell = cellsize
        assert lcell > 0, 'Invalid LCELL'
        assert mcell > 0, 'Invalid MCELL'
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [npix/2 + 1,npix/2 + 1] # honestly, I dont' understand if we need to +1 or not
        wcs.wcs.crval = [target.ra.deg, target.dec.deg]
        wcs.wcs.ctype = ['RA---SIN','DEC--SIN']
        wcs.wcs.cunit = ['deg','deg']
        wcs.wcs.cdelt = [-lcell.deg, mcell.deg]
        # The SIN projection is extended with the use
        #of the PV2_1 and PV2_2 FITS keywords set to a and -b, re-spectively.
        wcs.wcs.set_pv(((2,1,a), (2,2,-b)))
        wcs.wcs.timesys='TAI'
        wcs.wcs.trefpos = 'TOPOCENTER'
        wcs.wcs.timeunit = 's'
        mjdref = time.tai.mjd
        # mjdref integer and fractional part
        wcs.wcs.mjdref[:] = (int(mjdref), mjdref - int(mjdref))

        self.wcs2 = wcs

        # Astropy makes a mess if we use the WCS
        # We can add keywords to the header though
        w3hdr = wcs.to_header()
        w3hdr['NAXIS'] = 3
        w3hdr['CUNIT3'] = 's'
        w3hdr['CRVAL3'] = 0
        w3hdr['CRPIX3'] = 0
        w3hdr['CDELT3'] = inttime.to(u.second).value
        w3hdr['CTYPE3'] = 'TAI'
        self.wcs3_hdr = w3hdr
        

        # Now we make a copy of the 2d WCS and add a TAI
        # Dimension, just in case people want to use that.
        w3 = WCS(naxis=3)
        # Copy data from 2D WCS
        w3.wcs.crpix[:2] = wcs.wcs.crpix
        w3.wcs.crval[:2] = wcs.wcs.crval
        w3.wcs.ctype = ['RA--SIN','DEC--SIN','UTC'] # string proxies don't like to be sliced. Grrr
        w3.wcs.cunit = ['deg','deg','s']
        w3.wcs.cdelt[:2] = wcs.wcs.cdelt
        wcs.wcs.set_pv(((2,1,a), (2,2,-b)))
        wcs.wcs.timesys='TAI'
        wcs.wcs.trefpos = 'TOPOCENTER'
        wcs.wcs.timeunit = 's'
        mjdref = time.tai.mjd
        # mjdref integer and fractional part
        wcs.wcs.mjdref[:] = (int(mjdref), mjdref - int(mjdref))

        # Add time axis
        w3.wcs.crpix[2] = 0 # beginning of integration
        w3.wcs.crval[2] = 0 # 0 w.r.t. MJDREF
        w3.wcs.cdelt[2] = inttime.to(u.second).value # integration time
        w3.wcs.set_pv(((2,1,a), (2,2,-b)))
        w3.wcs.timesys='TAI'
        w3.wcs.trefpos = 'TOPOCENTER'
        w3.wcs.timeunit = 's'
        # mjdref integer and fractional part
        w3.wcs.mjdref[:] = (int(mjdref), mjdref - int(mjdref))

        self.wcs3 = w3

       
def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
