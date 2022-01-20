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
from . import sigproc
import itertools

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def onpick(event):
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print(thisline.get_label(), xdata[ind], ydata[ind])

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

    ffts = []

    antnos = [int(f.split('ak')[1][0:2]) for f in values.files]
    
    unique_antnos = sorted(set(antnos))
    print(len(unique_antnos), 'unique Antenna numbers', antnos)
    nrows = 3
    ncols = 4

    fig1, axes1 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
    fig2, axes2 = pylab.subplots(nrows, ncols, sharex=True, sharey=True)

    ant_ffts = {}
    for a in unique_antnos:
        ant_ffts[a] = []


    all_files = [sigproc.SigprocFile(fname) for fname in values.files]

        
    for f, antno in zip(all_files, antnos):
        tend = tstart + ntimes
        nelements = ntimes*f.nifs*f.nchans
        
        f.seek_data(f.bytes_per_element*tstart)
        if (f.nbits == 8):
            dtype = np.uint8
        elif (f.nbits == 32):
            dtype = np.float32
        v = np.fromfile(f.fin, dtype=dtype, count=nelements )
        v = v.astype(np.float)
        print('Nelements', nelements, 'Ntimes', ntimes, 'nchans', f.nchans, 'nifs', f.nifs, dtype, 'Actual length', len(v), 'Tsamp', f.tsamp)
    
        v.shape = (ntimes, f.nifs, f.nchans)
        beam = 0
        bd = v[:, beam, :]
        dm0 = bd.mean(axis=1)
        assert len(dm0) == ntimes
        dm0f = abs(np.fft.rfft(dm0))**2
        ffts.append(dm0f)
        
        aidx = unique_antnos.index(antno)
        ant_ffts[antno].append((f, dm0f))

    times = np.arange(ntimes)*f.tsamp
    fs = 1.0/f.tsamp
    ffts = np.array(ffts)
    nfftbins = ffts.shape[1]
    freqs = np.arange(nfftbins)/float(nfftbins)*fs/2.0
    #pylab.semilogy(freqs[1:], ffts.T[1:,:], label=fname, picker=1)
    #cid = pylab.gcf().canvas.mpl_connect('pick_event', onpick)

    #pylab.xlabel('Frequency (Hz)')
    #pylab.ylabel('abs(DM0)**2')

    pylab.figure()
    pylab.imshow(np.log10(ffts[:, 1:]), aspect='auto')

    pylab.figure()

    t0 = ant_ffts[antno][0][0].tstart
    marker = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))

    for ant in unique_antnos:
        bits = ant_ffts[ant]
        files = [b[0] for b in bits]
        times = np.array([f.tstart for f in files])
        ffts = np.array([b[1] for b in bits])
        iant = unique_antnos.index(ant)
        ax1 = axes1.flat[iant]
        ax2 = axes2.flat[iant]
        ax1.imshow(np.log10(ffts), aspect='auto')
        ax1.set_title('AK%02d' % ant)

        ax2.semilogy(freqs[1:], ffts[:, 1:].T)
        ax2.set_title('AK%02d' % ant)

        maxfreqs = [freqs[np.argmax(ffts[i, 1:])] for i in range(ffts.shape[0])]

        t = (times - t0)*24*60.

        pylab.plot(t, maxfreqs, label='AK%02d'%ant, marker=next(marker), alpha=0.5, ms=10)

        
    pylab.legend()
    pylab.xlabel('Time from start (mins)')
    pylab.ylabel('Frequency of FFT peak (Hz)')

    pylab.show()
    
        

    

if __name__ == '__main__':
    _main()
