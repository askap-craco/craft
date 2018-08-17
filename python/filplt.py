#!/usr/bin/env python
"""
PLot a filterbank - very simple

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d','--dm', type=float, default=0)
    parser.add_argument('-m', '--mjd', type=float, help='MJD to plot')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples', default=1024)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        plot(f, values)

def plot(f, values):
    s = sigproc.SigprocFile(f)
    fch1 = s.header['fch1']
    foff = s.header['foff']
    nchan = s.header['nchans']
    tsamp = s.header['tsamp']
    tstart = s.header['tstart']
    nsamp = values.nsamp
    if values.mjd is None:
        samp_start = nsamp/2
    else:
        samp_start = int(np.round((values.mjd - tstart)*86400.0/tsamp)) - nsamp/2
        if samp_start < 0:
            raise ValueError, 'Start sample is before start of filterbank start={}'.format(samp_start)
        if samp_start > s.file_size_elements:
            raise ValueError, 'End sample is after end of filterbank start={}'.format(samp_start)

    print values.mjd, tstart, tsamp, samp_start

    d = s[samp_start:samp_start+nsamp]
    # rescale to roughly 0 mean and 1 variance
    d = d.astype(float)
    d -= 128
    d /= 18.
    
    assert s.header['nbits'] == 8
    assert s.header['nifs'] == 1
    assert d.shape == (nsamp, nchan)
    channels = np.arange(nchan)*foff + fch1
    refchan = channels.min()

    dd = roll_dedisperse(d, channels, refchan , tsamp, values.dm)
    sn = dd.mean(axis=1)*np.sqrt(nchan)
    offset = np.arange(nsamp) - nsamp/2

    fig, ax = pylab.subplots(2,1, sharex=True)
    ax[0].imshow(dd.T, aspect='auto', extent=(offset[0], offset[-1], channels[0], channels[-1]), origin='lower')
    
    ax[1].plot(offset, sn)
    ax[0].set_title(f)
    #fig.title(f)
    ax[1].set_ylabel('S/N')
    ax[1].set_xlabel('Offset (samples) from %0.15f'%values.mjd)
    fig.savefig(f+'.png')

    if values.show:
        pylab.show()


def roll_dedisperse(d, channels, refchan, tsamp, dm):
    dd = d.copy()
    for c, f in enumerate(channels):
        delayms = 4.15*dm*((refchan/1e3)**-2 - (f/1e3)**-2)
        delaysamp = int(np.round(delayms/1e3/tsamp))
        dd[:, c] = np.roll(d[:, c], delaysamp)


    return dd
    
    

if __name__ == '__main__':
    _main()
