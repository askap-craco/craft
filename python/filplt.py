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
    parser.add_argument('-n','--nsamp', type=int, help='Number of sample', default=2048)
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
        samp_start = int(np.round((values.mjd -tstart)*86400.0/tsamp)) - nsamp/2

        print values.mjd, tstart, tsamp, samp_start
            
        if samp_start < 0:
            raise ValueError, 'Start sample is before start of filterbank'
        if samp_start > s.nsamples:
            raise ValueError, 'End sample is after end of filterbank'



    d = s[samp_start:samp_start+nsamp]
    # rescale to roughly 0 mean and 1 variance
    d = d.astype(float)
    d -= 128
    d /= 18.
    
    assert s.header['nbits'] == 8
    assert s.header['nifs'] == 1
    assert d.shape == (nsamp, nchan)

    dd = roll_dedisperse(d, nchan, fch1, foff, tsamp, values.dm)
    sn = dd.mean(axis=1)*np.sqrt(nchan)

    fig, ax = pylab.subplots(2,1, sharex=True)
    ax[0].imshow(dd.T, aspect='auto')
    
    ax[1].plot(sn)
    ax[0].set_title(f)
    #fig.title(f)
    ax[1].set_ylabel('S/N')
    fig.savefig(f+'.png')

    if values.show:
        pylab.show()


    pylab.close(fig)

        

    

def roll_dedisperse(d, nchan, fch1, foff, tsamp, dm):
    dd = d.copy()
    for c in xrange(nchan):
        f = fch1 + foff*c
        delayms = 4.15*dm*((fch1/1e3)**-2 - (f/1e3)**-2)
        delaysamp = -int(abs(np.round(delayms/1e3/tsamp)))
        dd[:, c] = np.roll(d[:, c], delaysamp)


    return dd
    
    

if __name__ == '__main__':
    _main()
