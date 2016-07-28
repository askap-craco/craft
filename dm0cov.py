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
import sigproc
import itertools
import re
import glob

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def onpick(event):
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print thisline.get_label(), xdata[ind], ydata[ind]

def load_beams(path, tstart, ntimes):
    files = sorted(glob.glob(os.path.join(path, '*.fil')))
    if len(files) == 0:
        raise ValueError('No files in path  %s' % path)
    data = None
    for fname in files:
        f = sigproc.SigprocFile(fname)
        tend = tstart + ntimes
        nelements = ntimes*f.nifs*f.nchans
        
        f.seek_data(f.bytes_per_element*tstart)
        if (f.nbits == 8):
            dtype = np.uint8
        elif (f.nbits == 32):
            dtype = np.float32
        
        v = np.fromfile(f.fin, dtype=dtype, count=nelements )
        v.shape = (ntimes, f.nifs, f.nchans)

        if f.nifs == 1:
            nifs = len(files)
            ifslice = slice(0, nifs)
            if data is None:
                data = np.zeros((ntimes, nifs, f.nchans))

            ifnum = int(fname.split('.')[-2])
            data[:, ifnum, :] = v[:, 0, :]

        else:
            data = v


        
    return data


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

    nrows = 3
    ncols = 4

    fig1, axes1 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
    fig2, axes2 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
    fig3, axes3 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
    fig4, axes4 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)


    antbeams = {}

    axes1.flat[0].set_ylabel('Beam number')
    axes2.flat[0].set_ylabel('Bandpass mean')
    axes3.flat[0].set_ylabel('Bandpass std')
    axes4.flat[0].set_ylabel('FFT(DM0)')
    chan = 290
    
    
    for i, antdir in enumerate(values.files):
        antpath = os.path.join(antdir, 'C000')
        try:
            beams = load_beams(antpath, tstart, ntimes)
        except:
            logging.exception('Could not load beams at path %s', antdir)
        
        # shape is time, beam, freq
        dm0 = beams.mean(axis=2)
        dm0 = beams[:, :, chan]
        cov = np.cov(dm0.T)
        cov /= np.diag(cov)
        antname = antdir.replace('/','')

        ax = axes1.flat[i]
        ax.imshow(np.log10(cov), aspect='auto')
        ax.set_title(antname)

        ax2 = axes2.flat[i]
        ax2.plot(beams.mean(axis=0).T)
        ax2.set_title(antname)

        #ax2b = ax2.twinx()

        ax3 = axes3.flat[i]
        ax3.plot(beams.std(axis=0).T)
        ax3.set_title(antname)

        ax4 = axes4.flat[i]
        dm0 = beams.mean(axis=2)
        df = np.abs(np.fft.rfft(dm0.T))**2
        ax4.semilogy(df[1:].T)

        ax4.set_title(antname)

        


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
