#!/usr/bin/env python
"""
CRACO utilities

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

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

def time_blocks(vis, nt):
    '''
    Generator that returns nt time blocks of the given visiblity table


    :returns: Directionary, keyed by baseline ID, value = np.array size (nchan, nt) dtype=complex64
    :see: bl2ant()
    '''

    nrows = vis.size
    nchan = vis[0].data.shape[-3]
    d = {}
    t = 0
    d0 = vis[0]['DATE']
    for irow in xrange(nrows):
        row = vis[irow]
        if row['DATE'] > d0:
            t += 1
            if t == nt:
                yield d
                d = {}
                d0 = row['DATE']
                t = 0

        blid = row['BASELINE']
        if blid not in d.keys():
            d[blid] = np.zeros((nchan, nt), dtype=np.complex64)

        d[blid][:, t].real = row.data[...,0].reshape(nchan)
        d[blid][:, t].imag = row.data[...,1].reshape(nchan)

        
    if len(d) > 0:
        if t < nt:
            warnings.warn('Final integration only contained {} of {} samples'.format(t, nt))
        yield d

def grid(uvcells, Npix):
    '''
    Grid the data onto an Npix grid
    '''
    np2 = int(float(Npix)/2.)
    g = np.zeros((Npix, Npix), dtype=np.complex)

    for b in uvcells:
        upix, vpix = b.uvpix
        g[vpix, upix] += b.nchan
        g[Npix-vpix, Npix-upix] += b.nchan

    return g


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
