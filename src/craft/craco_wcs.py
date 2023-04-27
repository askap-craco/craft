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

    # Steve's equation 5 has a minus sign in front of tan(Z)
    # I found I needed to remove it because otherwise my UVWs didn't fall on teh correct plane
    # see CRACO-115 for details.
    b = tan(Z) * cos(chi) # Equation 5

    return (a, b)

class CracoWCS:
    def __init__(self, target: SkyCoord, time: Time, cellsize, npix:int, inttime, site: EarthLocation=None):
        '''

        Has 2 useful attributes
        wcs2 is the 2D WCS for an image centered on the given phase center, which includes the generalised SIN slant orthographic projection described in Ord et al. 2010

        wcs3 is wcs2 but includes the time axis, which is set to start at the given time and increases by inttime every pixel. This is useful for encoding a movie of images
        
        :param: target - SkyCoord look direction phase center
        :time: Time to set WCS to
        :cellsize: (lcell, mcell) cell size 2 tuple as astropy Angle
        :npix: Numebr of pixels (int)
        :inttime: Integration tiem as astropy time Quanitty
        :site: if None defaults to EarthLocation.of_site('ASKAP') else the site
        '''
        if site is None:
            try:
                site = EarthLocation.of_site('ASKAP')
            except:
                # stupid old versions of stupid astropy don't stupid have askap
                # stupid
                site = EarthLocation.from_geocentric(-2556084.65961682,
                                                     5097398.3818179,
                                                     -2848424.06141933,
                                                     u.meter)

        time.location = site # needed to remove warning
        
        lst = time.sidereal_time('apparent', site.lon)
        altaz = target.transform_to(AltAz(obstime=time, location=site))
        assert altaz.alt.deg > 12, f'target is below ASKAP horizon {altaz} {lst}'
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
        wcs = WCS(naxis=3)
        wcs.wcs.crpix = [npix/2 + 1,npix/2 + 1, 1] # honestly, I dont' understand if we need to +0.5 or not, or 1. 
        wcs.wcs.crval = [target.ra.deg, target.dec.deg, 0]
        wcs.wcs.ctype = ['RA---SIN','DEC--SIN', 'TAI']
        wcs.wcs.cunit = ['deg','deg','s']
        wcs.wcs.cdelt = [-lcell.deg,
                         mcell.deg,
                         inttime.to(u.second).value]
        # From Ord 2010...
        # The SIN projection is extended with the use
        #of the PV2_1 and PV2_2 FITS keywords set to a and -b, re-spectively.
        wcs.wcs.set_pv(((2,1,a), (2,2,-b)))
        wcs.wcs.timesys='TAI'
        # wcs.wcs.trefpos = 'topocenter' astropy prints warning if we set trefpos anod don't set obsgeo - and I don't know how to do that. it's probably not too important
        wcs.wcs.timeunit = 's'
        mjdref = time.tai.mjd
        # mjdref integer and fractional part
        wcs.wcs.mjdref[:] = (int(mjdref), mjdref - int(mjdref))

        # wcs3 is 3D wcs
        self.wcs3 = wcs

        # drop time axis to make 2d wcs
        self.wcs2 = wcs.dropaxis(2)


    
    @staticmethod
    def from_plan(plan):
        return CracoWCS(plan.phase_center,
                        plan.tstart,
                        plan.lmcell,
                        plan.npix,
                        plan.tsamp_s)

       
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
