#!/usr/bin/env python
"""
Makes dedispersed visibilities suitable for the CRACO imaging pipelien that does gridding, FFT and boxcar.


Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from astropy.io import fits
from scipy import constants

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def bl2ant(bl):
    '''
    Convert baseline to antena numbers according to UV fits convention
    Antenna numbers start at 1 and:

    baseline = 256*ant1 + ant2

    :see: http://parac.eu/AIPSMEM117.pdf

    :returns: (ant1, ant2) as integers

    >>> bl2ant(256*1 + 2)
    (1, 2)

    >> bl2ant(256*7 + 12)
    (7, 12)
    '''
    ibl = int(bl)
    a1 = ibl // 256
    a2 = ibl % 256

    assert a1 >= 1
    assert a2 >= 1

    return (a1, a2)

def runidxs(x):
    ''' 
    Return the indexes of the start an end of a list numbers that might be equal

    '''
    istart = 0
    for i in xrange(1, len(x)):
        if x[i] != x[istart]:
            yield (istart, i-1)
            istart = i
            
    yield (istart, i)

def arcsec2rad(strarcsec):
    return np.radians(float(strarcsec)/3600.)

class BaselineCell(object):
    def __init__(self, blid, uvpix, chan_start, chan_end, freqs):
        self.blid = blid
        self.uvpix = uvpix
        self.chan_start = chan_start
        self.chan_end = chan_end
        self.freqs = freqs
        self.a1, self.a2 = bl2ant(blid)

    @property
    def nchan(self):
        return self.chan_end - self.chan_start + 1

def grid(uvcells, values):
    Npix = values.npix
    np2 = int(float(Npix)/2.)
    g = np.zeros((Npix, Npix), dtype=np.complex)

    for b in uvcells:
        upix, vpix = b.uvpix
        g[vpix, upix] += b.nchan
        g[Npix-vpix, Npix-upix] += b.nchan

    return g


def image_fft(g):
    '''
    Do the complex-to-complex imaging FFT with the correct shifts and correct inverse thingy
    If g.shape = (Npix, Npix) then we assume the center of the UV plane is at
    Npix/2, Npix/2 = DC
    '''
    # The old version was incorrect!
    #cimg = np.fft.fftshift(np.fft.ifft2(g)).astype(np.complex64)
    
    cimg = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g))).astype(np.complex64)
    return cimg/np.prod(cimg.shape)
            

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Makes dedispersed visibilities', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--cell', help='Image cell size (arcsec)', default='10,10')
    parser.add_argument(dest='files', nargs='?')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Opening file %s', values.files)
    hdul = fits.open(values.files)

    vis = hdul[0].data

    # get data for first integration to work out UVW coordinatse
    d0 = vis[0]['DATE']
    baselines = {}
    for i in xrange(vis.size):
        row = vis[i]
        baselines[row['BASELINE']] = row
        if row['DATE'] != d0:
            break

    hdr = hdul[0].header
    assert hdr['CTYPE4'] == 'FREQ', 'Cant find freq axis'
    fch1 = hdr['CRVAL4']
    foff = hdr['CDELT4']
    ch1 = hdr['CRPIX4']
    assert ch1 == 1.0, 'Unexpected initial frequency'
    nchan = vis[0].data.shape[-3]
    freqs = np.arange(nchan)*foff + fch1 # Hz
    lcell, mcell = map(arcsec2rad, values.cell.split(','))
    Npix = values.npix
    lfov = lcell*Npix
    mfov = mcell*Npix
    ucell, vcell = 1./lfov, 1./mfov
    logging.info('Fch1=%f foff=%f uvcell=%s', fch1, foff, (ucell, vcell))


    uvcells = []

    for blid, bldata in baselines.iteritems():
        #UU, VV WW are in seconds
        ulam = bldata['UU'] * freqs
        vlam = bldata['VV'] * freqs
        upix = np.round(ulam/ucell + Npix/2).astype(int)
        vpix = np.round(vlam/vcell + Npix/2).astype(int)

        pylab.plot(ulam/1e3, vlam/1e3)

        uvpos = list(zip(upix, vpix))
        for istart, iend in runidxs(uvpos):
            uvpix = uvpos[istart]
            assert uvpos[istart] == uvpos[iend]
            b = BaselineCell(blid, uvpix, istart, iend, freqs[istart:iend+1])
            uvcells.append(b)

    uvcells = sorted(uvcells, key=lambda b:b.uvpix)
    d = np.array([(f.a1, f.a2, f.uvpix[0], f.uvpix[1], f.chan_start, f.chan_end) for f in uvcells], dtype=np.uint32)
    np.savetxt('grid.txt', d, fmt='%d',  header='ant1, ant2, u(pix), v(pix), chan1, chan2')

    pylab.xlabel('U(klambda)')
    pylab.ylabel('V(klambda)')

    pylab.figure()
    fix, ax = pylab.subplots(1,2)
    g = grid(uvcells, values)
    ax[0].imshow(abs(g), aspect='auto', origin='lower')
    ax[1].imshow(image_fft(g).real, aspect='auto', origin='lower')
    pylab.show()
    
        
    
    
        

if __name__ == '__main__':
    _main()
