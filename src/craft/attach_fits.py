#!/usr/bin/env python
"""
Takes an ascii and attaches a fits header to it

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from  astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime
from aces.footprint_class import Footprint


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def skycoord(s):
    return SkyCoord(s, frame='icrs', unit='deg')

class RyanMultinestImage(object):
    def __init__(self, f):
        d = np.loadtxt(f)
        with open(f, 'rU') as fin:
            for iline, line in enumerate(fin):
                if iline == 0:
                    self.nl, self.nm = list(map(int, np.array(line.split())[[2,4]]))
                elif iline == 1:
                    self.lmin, self.dl, self.mmin, self.dm = list(map(float, np.array(line.split())[[2,4,6,8]]))

        img = np.zeros((self.nm, self.nl))
        for idx in range(d.shape[0]):
            l, m = list(map(int, d[idx, 0:2]))
            img[m, l] = d[idx, 4]

        self.img = img
        lvec = np.arange(self.lmin, self.lmin + self.dl*self.nl, self.dl)
        mvec = np.arange(self.mmin, self.mmin + self.dm*self.nm, self.dm)
        print(lvec.min(), lvec.max(), mvec.min(), mvec.max())
        self.gx, self.gy = np.meshgrid(lvec, mvec)


class GaussImage(object):
    def __init__(self):
        x = np.linspace(-2.5,2.5,100)
        self.gx, self.gy = np.meshgrid(x,x)
        self.img = np.exp(-(self.gx - 2.*0)**2/0.1 + -(self.gy - 0.5*0)**2/0.1)
        self.nl = len(x)
        self.nm = len(x)
        self.dl = x[1] - x[0]
        self.dm = self.dl
        self.lmin = min(x)
        self.mmin = min(x)


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Parses ryans multinest output and adds a fits header', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--refpos', help='Center pixel', type=skycoord)
    parser.add_argument('--resolution', help='Resolution in degrees', type=float, default=1.)
    parser.add_argument('--pa', type=float, help='PA degrees', default=0)
    parser.add_argument(dest='file', nargs='?')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.file is None:
        rf = GaussImage()
    else:
        rf = RyanMultinestImage(values.file)

    v = rf.img
    fp = Footprint.named('square_6x6',np.radians(0.9),np.radians(170+45.))
    lms = np.degrees(np.array(fp.offsetsRect))

    pylab.imshow(v, aspect='auto', origin='bottom', extent=(rf.gx.min(), rf.gx.max(), rf.gy.min(), rf.gy.max()))
    ax = pylab.gca()
    print(lms.shape)
    pylab.scatter(lms[:, 0], lms[:, 1])
    for ibeam, lm in enumerate(lms):
        pylab.text(lm[0], lm[1], '%d'%ibeam)
        cir = plt.Circle(lm, 0.9)
        #ax.add_artist(cir)

    pylab.xlabel('l (deg)')
    pylab.ylabel('m (deg)')


    hdu = fits.PrimaryHDU(v)
    logging.debug('Refpos is %s', values.refpos.to_string('hmsdms'))
    
    hdulist = fits.HDUList([hdu])
    if values.file is None:
        values.file = 'test'

    fout = values.file + '.fits'
    print('WRiting', fout)
    hdr = hdu.header
    hdr['infile'] = values.file
    hdr['object'] = 'FRB170107'
    hdr['imtype'] = 'Error region'
    hdr['telescop'] = 'ASKAP'
    hdr['instrume'] = 'CRAFT'
    
    hdr['history'] = 'Created with attached_fits.py on %s UT' % (datetime.utcnow().isoformat())
    xoff = rf.lmin/rf.dl
    yoff = rf.mmin/rf.dm
    hdr['ctype1'] = 'RA---TAN'
    hdr['crpix1'] =  -xoff - 0.5
    hdr['crval1'] = values.refpos.ra.deg
    hdr['cdelt1'] = rf.dl
    
    hdr['ctype2'] = 'DEC--TAN'
    hdr['crpix2'] = -yoff  - 0.5
    hdr['crval2'] = values.refpos.dec.deg
    hdr['cdelt2'] = rf.dm

    hdr['CROTA1'] = values.pa
    hdr['CROTA2'] = values.pa

    print(dir(hdr))

    for h in list(hdr.items()):
        print(h)

    hdulist.writeto(fout, overwrite=True)

    pylab.show()

        

if __name__ == '__main__':
    _main()
