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

def annotate(fig, title, xlabel, ylabel):
    fig.text( 0.5, 0.98,title, ha='center', va='top')
    fig.text(0.5, 0.02, xlabel, ha='center', va='bottom')
    fig.text(0.02, 0.5, ylabel, rotation=90, ha='center', va='top')
    fout='{}.png'.format(title.replace(' ', ''))
    print 'Saving', fout
    fig.savefig(fout)
    

def commasep(s):
    return map(int, s.split(','))

def floatcommasep(s):
    return map(float, s.split(','))


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-t', '--times', help='Integration range to plot', type=commasep)
    parser.add_argument('--nxy', help='number of rows,columns in plots', type=commasep)
    parser.add_argument('--imzrange', help='Z range for dynamic spectrum', type=floatcommasep)
    parser.set_defaults(verbose=False, nxy="3,3")
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.times:
        bits = values.times
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

    nrows, ncols = values.nxy

    ffts = []
    if values.imzrange is None:
        imzmin = None
        imzmax = None
    else:
        imzmin, imzmax = values.imzrange

    beams, files = load_beams(values.files, tstart, ntimes, return_files=True)
    print 'Loaded beams', beams.shape
    ntimes, nbeams, nfreq = beams.shape

    f0 = files[0]
    mjdstart = f0.tstart
    tsamp = f0.tsamp
    fch1 = f0.fch1
    foff = f0.foff
    src_raj = f0.src_raj
    src_dej = f0.src_dej
    freqs = np.linspace(fch1, fch1 + nfreq*foff, nfreq, endpoint=True)
    assert(len(freqs) == nfreq)

    if foff < 0:
        origin = 'upper'
        im_extent = (0, ntimes*tsamp, fch1 + foff*nfreq, fch1)
    else:
        origin = 'lower'
        im_extent = (0, ntimes*tsamp, fch1, fch1 + foff*nfreq)

    fig1, axes1 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig2, axes2 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig3, axes3 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig4, axes4 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig5, axes5 = subplots(nrows, ncols, sharex=True, sharey=True)
    fig6, axes6 = subplots(nrows, ncols, sharex=True, sharey=True)


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
        ax6 = axes6.flat[i]
        bi = beams[:, i, :]
        print 'bi', bi.shape
        ntimes, nfreq = bi.shape
        
        ax1.imshow(bi.T, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent)
        ax1.text(0.98, 0.98, 'B{:02d}'.format(i), va='top', ha='right', transform=ax1.transAxes)
        beam_mean = bi.mean(axis=0)
        beam_std = bi.std(axis=0)
        bmm = np.tile(beam_mean, (ntimes, 1))
        bsm = np.tile(beam_std, (ntimes, 1))
        bi_znorm = (bi - bmm)/bsm
        print 'Znorm', bi_znorm.shape
        beam_kurtosis = np.mean((bi_znorm)**4, axis=0)/np.mean((bi_znorm)**2, axis=0)**2
        print 'kurt shape', beam_kurtosis.shape

        ax2.plot(freqs, beam_mean)
        ax3.plot(freqs, beam_std)
        ax6.plot(freqs, beam_kurtosis)

        dm0 = bi.mean(axis=1)

        dm0f = abs(np.fft.rfft(dm0, axis=0))**2
        bf = abs(np.fft.rfft(bi, axis=0).T)**2
        print 'bf', bf.shape
        #ax4.plot(np.log10(bf[0::16, :]).T)
        ax4.imshow(np.log10(bf.T)[1:, :], aspect='auto')
        
        ax5.plot(np.log10(dm0f[1:]))

    annotate(fig1, 'Dynamic spectrum', 'Time (s)', 'Frequency (MHz)')
    annotate(fig2, 'Mean bandpass', 'Frequency (MHz)','Mean bandpass')
    annotate(fig3, 'Bandpass stdDev', 'Frequency (MHz)','Bandpass StdDev')
    annotate(fig4, 'FFT of all channels', 'Digital frequency ', 'Channel')
    annotate(fig5, 'FFT of DM0', 'Digital Frequency', 'FFT (dB)')
    annotate(fig6, 'Kurtosis', 'Frequency (MHz)', 'Kurtosis')


    pylab.show()

if __name__ == '__main__':
    _main()
