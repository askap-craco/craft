#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import itertools
import re
from .craftobs import load_beams
from .plotutil import subplots as mysubplots
from . import plotutil

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-t', '--times', help='Integration range to plot')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.times:
        bits = list(map(int, values.times.split(',')))
        if len(bits) == 1:
            tstart = bits[0]
            ntimes = 128*8
        elif len(bits) == 2:
            tstart, ntimes = bits
        else:
            raise ValueError('Invalid times: %s', values.times)
    else:
        tstart = 0
        ntimes = 128*8

    nrows = 2
    ncols = 5
    
    fig1, axes1 = mysubplots(nrows, ncols, sharex=True, sharey=True)
    fig2, axes2 = mysubplots(nrows, ncols, sharex=True, sharey=True)
    fig3, axes3 = mysubplots(nrows, ncols, sharex=True, sharey=True) 
    fig4, axes4 = mysubplots(nrows, ncols, sharex=True, sharey=True)
    fig5, axes5 = mysubplots(nrows, ncols, sharex=True, sharey=True)
    fig6, axes6 = mysubplots(nrows, ncols, sharex=True, sharey=True)


    antbeams = {}

    axes1.flat[0].set_ylabel('Beam number')
    axes2.flat[0].set_ylabel('Bandpass mean')
    axes3.flat[0].set_ylabel('Bandpass std')
    axes4.flat[0].set_ylabel('FFT(DM0)')
    axes5.flat[0].set_ylabel('Channel')
    axes6.flat[0].set_ylabel('FFT(chans)')

    chan = 10
    beamno = 0

    fs = 1e6*32./27.
    tint = 1500./fs # seconds
    fint = 1./tint
    print('fs', fs, 'Tint', tint, fint)

    
    for i, antdir in enumerate(values.files):
        antpath = os.path.join(antdir, 'C000/')
        try:
            beams = load_beams(antpath, tstart, ntimes)
            print('Loaded beams', beams.shape, i, antdir)
        except:
            logging.exception('Could not load beams at path %s', antdir)
            continue
        
        # shape is time, beam, freq
        
        print(beams.shape)
        ntimes, nbeams, nfreqs = beams.shape
        #dm0 = beams.mean(axis=2)

        dm0 = beams[:, :, 280:336].mean(axis=2)
        cov = np.cov(dm0.T)
        cov /= np.diag(cov)
        antname = antdir.replace('/','')

        ax = axes1.flat[i]
        ax.imshow(np.log10(cov), aspect='auto')
        ax.set_title(antname)

        ax2 = axes2.flat[i]
        mean_bandpass = beams.mean(axis=0).T
        print('mean bp', mean_bandpass.shape)
        ax2.plot(mean_bandpass)
        ax2.set_title(antname)

        #ax2b = ax2.twinx()

        ax3 = axes3.flat[i]
        beam_std = beams.std(axis=0).T
        for b in range(nbeams):
            ax3.plot(beam_std[:, b], label='Beam %d' % b, picker=3)

        ax3.set_title(antname)

        ax4 = axes4.flat[i]
        df = np.abs(np.fft.rfft(dm0, axis=0))**2
        nft = df.shape[0]
        freqs = np.arange(nft)/float(nft)*fint/2.
        print('Freqs', freqs.shape, df.shape)
        for b in range(nbeams):
            ax4.loglog(freqs[1:], df[1:, b], label='Beam %d' % b, picker=3)

        ax4.set_title(antname)

        ax5 = axes5.flat[i]
        b0 = beams[:, beamno, :]
        print(b0.shape)
        ax5.imshow(b0.T, aspect='auto')
        ax5.set_title(antname)

        ax6 = axes6.flat[i]
        bf = np.abs(np.fft.rfft(b0, axis=0))**2
        ax6.imshow(np.log10(df), aspect='auto')
        ax6.set_title(antname)


    ax.set_xlabel('Beam number')
    ax.set_ylabel('Beam nunber')
    ax2.set_xlabel('Channel')
    ax3.set_xlabel('Channel')
    ax4.set_xlabel('Frequency (Hz)')
    ax5.set_xlabel('Integration number')
    ax6.set_xlabel('Channel number')

    #cb = ax.colorbar()
    #cb.set_label('Covariance (dB)')

    fig1.savefig('ant_cov_c%d.png' % chan)
    fig2.savefig('ant_bpmean.png')
    fig3.savefig('ant_bpstd.png')
    fig4.savefig('ant_dm0fft.png')

    for f in (fig2, fig3, fig4):
        f.canvas.mpl_connect('pick_event', plotutil.onpick)



    pylab.show()
    

if __name__ == '__main__':
    _main()
