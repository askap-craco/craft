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
from craftobs import load_beams
from plotutil import subplots as mysubplots

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
        bits = map(int, values.times.split(','))
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

    nrows = 1
    ncols = 1
    
    fig1, axes1 = mysubplots(nrows, ncols)
    fig2, axes2 = mysubplots(nrows, ncols)
    fig3, axes3 = mysubplots(nrows, ncols)
    fig4, axes4 = mysubplots(nrows, ncols)
    fig5, axes5 = mysubplots(nrows, ncols)


    antbeams = {}

    axes1.flat[0].set_ylabel('Beam number')
    axes2.flat[0].set_ylabel('Bandpass mean')
    axes3.flat[0].set_ylabel('Bandpass std')
    axes4.flat[0].set_ylabel('FFT(DM0)')
    axes5.flat[0].set_ylabel('Channel')

    chan = 300
    beamno = 0 
    
    
    for i, antdir in enumerate(values.files):
        antpath = os.path.join(antdir, 'C000/')
        try:
            beams = load_beams(antpath, tstart, ntimes)
            print 'Loaded beams', beams.shape, i, antdir
        except:
            logging.exception('Could not load beams at path %s', antdir)
            continue
        
        # shape is time, beam, freq
        
        print beams.shape
        dm0 = beams.mean(axis=2)

        #dm0 = beams[:, :, chan]
        cov = np.cov(dm0.T)
        cov /= np.diag(cov)
        antname = antdir.replace('/','')

        ax = axes1.flat[i]
        ax.imshow(np.log10(cov), aspect='auto')
        ax.set_title(antname)

        ax2 = axes2.flat[i]
        mean_bandpass = beams.mean(axis=0).T
        print 'mean bp', mean_bandpass.shape
        ax2.plot(mean_bandpass)
        ax2.set_title(antname)

        #ax2b = ax2.twinx()

        ax3 = axes3.flat[i]
        ax3.plot(beams.std(axis=0).T)
        ax3.set_title(antname)

        ax4 = axes4.flat[i]
        dm0 = beams.mean(axis=2)
        df = np.abs(np.fft.rfft(dm0.T))**2
        ax4.semilogy(df[1:].T)
        #ax4.plot(dm0)
        ax4.set_title(antname)

        ax5 = axes5.flat[i]
        b0 = beams[:, beamno, :]
        print b0.shape
        ax5.imshow(beams[:,0,:].T, aspect='auto')
        ax5.set_title(antname)



    ax.set_xlabel('Beam number')
    ax.set_ylabel('Beam nunber')
    ax2.set_xlabel('Channel')
    ax3.set_xlabel('Channel')

    #cb = ax.colorbar()
    #cb.set_label('Covariance (dB)')

    fig1.savefig('ant_cov_c%d.png' % chan)
    fig2.savefig('ant_bpmean.png')
    fig3.savefig('ant_bpstd.png')
    fig4.savefig('ant_dm0fft.png')


    pylab.show()
    

if __name__ == '__main__':
    _main()
