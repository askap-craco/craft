#!/usr/bin/env python
"""
Plots a lot of beams

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
from craftobs import load_beams
from plotutil import subplots

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def onpick(event):
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print thisline.get_label(), xdata[ind], ydata[ind]

def annotate(fig, title, xlabel, ylable):
    fig.add_text(title, 0.5, 0.98, ha='center', va='top')
    fig.add_text(xlabel, 0.5, 0.02, ha='center', va='bottom')
    fig.add_text(ylabel, 0.02, 0.5, angle=90, ha='center', va='top')
    

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

    ffts = []
    nrows = 6
    ncols = 12

    fig1, axes1 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig2, axes2 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig3, axes3 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig4, axes4 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig5, axes5 = subplots(nrows, ncols, sharex=True, sharey=True)

    beams = load_beams(values.files[0], tstart, ntimes)
    print 'Loaded beams', beams.shape
    ntimes, nbeams, nfreq = beams.shape

    pylab.figure()
    chan = 150
    dm0 = beams.mean(axis=2)
    dm0 = beams[:, :, chan]
    dm0 -= dm0.mean(axis=0)
    
    cov = np.cov(dm0.T)
    #b0 = beams[:, 0, :]
    #cov = np.cov(b0.T)
    cov /= np.diag(cov)
    pylab.imshow(np.log10(cov), interpolation='none')
    pylab.xlabel('Beam number')
    pylab.ylabel('Beam number')
    pylab.colorbar()

    nplots = min(nrows*ncols, nbeams)

    for i in xrange(nplots):

        ax1 = axes1.flat[i]
        ax2 = axes2.flat[i]
        ax3 = axes3.flat[i]
        ax4 = axes4.flat[i]
        ax5 = axes5.flat[i]
        bi = beams[:, i, :]
        print 'bi', bi.shape
        ntimes, nfreq = bi.shape
        
        ax1.imshow(bi.T, aspect='auto')
        ax1.set_title('Beam {}'.format(i))

        ax2.plot(bi.mean(axis=0))
        ax3.plot(bi.std(axis=0))
        dm0 = bi.mean(axis=1)

        dm0f = abs(np.fft.rfft(dm0, axis=0))**2
        bf = abs(np.fft.rfft(bi, axis=0).T)**2
        print 'bf', bf.shape
        #ax4.plot(np.log10(bf[0::16, :]).T)
        ax4.imshow(np.log10(bf.T)[1:, :], aspect='auto')
        
        ax5.plot(np.log10(dm0f[1:]))

    

    annotate(fig1, 'Dynamic spectrum', 'Time', 'Channel')
    annotate(fig2, 'Mean bandpass', 'Channel','Mean bandpass')
    annotate(fig3, 'Bandpass stdDev', 'Channel','Bandpass StdDev')
    annotate(fig4, 'FFT of all channels', 'Digital frequency ', 'Channel')
    annotate(fig5, 'FFT of DM0', 'Digital Frequency', 'FFT (dB)')


    pylab.show()

if __name__ == '__main__':
    _main()
