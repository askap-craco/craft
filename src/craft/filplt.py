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
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d','--dm', type=float, default=0)
    parser.add_argument('-m', '--mjd', type=float, help='MJD to plot')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('-c','--cand-file', help='Candidate file')
    parser.add_argument('-a','--all', action='store_true', help='Plot all beams if a candidate file is specified, not just the detection beam')
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples', default=1024)
    parser.add_argument('-r','--rescale', action='store_true', help='Rescale')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.cand_file is None:
        for f in values.files:
            plot(f, values)
    else:
        candidates = np.loadtxt(values.cand_file, dtype=np.float64)
        if len(candidates) == 0:
            raise ValueError('No candidates in file {}'.format(values.cand_file))
        if len(candidates.shape) == 1:
            candidates.shape = (1, -1)
            
        ncand = candidates.shape[0]

        for c in range(ncand):
            # 103.02 70298 60.9747 11 685 475.36 15 58374.2622446670 694.09
            mjd=candidates[c,7]
            dm=candidates[c,5]
            beam=int(candidates[c,6])
            if values.all:
                beam_files = values.files
            else:
                beam_files = [f for f in values.files if int(f.split('.')[-2]) == beam]
            print(beam_files)

            for f in beam_files:
                print(f, mjd, dm, beam)
                plot(f, values, mjd, dm)

def plot(f, values, mjd=None, dm=None):
    s = sigproc.SigprocFile(f)
    fch1 = s.header['fch1']
    foff = s.header['foff']
    nchan = s.header['nchans']
    tsamp = s.header['tsamp']
    tstart = s.header['tstart']
    nsamp = values.nsamp
    if mjd is None:
        mjd = values.mjd

    if dm is None:
        dm = values.dm
        
    if mjd is None:
        samp_start = nsamp/2
        mjd = tstart + samp_start*tsamp/86400.0
    else:
        samp_start = int(np.round((values.mjd -tstart)*86400.0/tsamp)) - nsamp/2

        print(values.mjd, tstart, tsamp, samp_start)
            
        if samp_start < 0:
            raise ValueError('Start sample is before start of filterbank')
        if samp_start > s.nsamples:
            raise ValueError('End sample is after end of filterbank')

    d = s[samp_start:samp_start+nsamp]
    # rescale to roughly 0 mean and 1 variance
    d = d.astype(float)
    assert s.header['nifs'] == 1
    assert d.shape == (nsamp, nchan)
    channels = np.arange(nchan)*foff + fch1

    if values.cand_file and values.cand_file.endswith('.finf'):
        refchan = 1e300
    else:
        refchan = channels.min()
        
    if values.rescale:
        d -= d.mean(axis=0)
        d /= d.std(axis=0)

    dd = roll_dedisperse(d, channels, refchan , tsamp, dm)
    sn = dd.mean(axis=1)*np.sqrt(nchan)
    offset = np.arange(nsamp) - nsamp/2

    fig, ax = pylab.subplots(2,1, sharex=True)
    ax[0].imshow(dd.T, aspect='auto', extent=(offset[0], offset[-1], channels[0], channels[-1]), origin='lower')
    
    ax[1].plot(offset, sn)
    ax[0].set_title(f)
    #fig.title(f)
    ax[1].set_ylabel('S/N')
    
    if mjd is None:
        lbl = 'Offset (samples)'
    else:
        lbl = 'Offset (samples) from %0.9f'%mjd
        
    ax[1].set_xlabel(lbl)
    try:
        fig.savefig(f.replace('.fil','.png'))
    except:
        logging.exception('Could not save figure')

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
