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
import matplotlib.gridspec as gridspec
import sigproc

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

def commasep(s):
    return map(int, s.split(','))

def floatcommasep(s):
    return map(float, s.split(','))

def next_divisor(x, n):
    for t in xrange(x, n/2):
        if n % t == 0:
            return t

def divisors(n):
    d = [t for t in xrange(1, n/2 + 1) if n % t == 0]
    return d

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--outfile', help='Output file')
    parser.add_argument('-A', '--suffix', help='OUtput file suffix')
    parser.add_argument('-s', '--show',  action='store_true', help='show plots')

    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False, nxy="1,1")
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        filclean(f, values)

def filclean(fname, values):
    f = sigproc.SigprocFile(fname)
    tstart = 0
    ntimes = 4096
    tend = tstart + ntimes
    nelements = ntimes*f.nifs*f.nchans
    assert f.nifs == 1
    #f.seek_data(f.bytes_per_element*tstart)
    spout = None
    outf = None
    if values.outfile:
        outf = values.outfile
    elif values.suffix:
        outf = fname + '.' + values.suffix
        
    if outf:
        spout = sigproc.SigprocFile(outf, 'w', header=f.header)

    while True:
        dtype = np.uint8
        v = np.fromfile(f.fin, dtype=dtype, count=nelements )
        if len(v) != ntimes*f.nchans:
            # finished
            break
        v.shape = (ntimes, f.nchans)
        #v -= 128.
        #v /= 18.

        dm0 = v.mean(axis=1)
        offset = dm0.mean()
        # shape = time, channel
        vs = v - dm0[:, None] + offset # broadcast DM0
        print 'SHapes', v.shape, dm0[:, None].shape, vs.dtype


        if spout:
            vs.astype(np.uint8).flatten().tofile(spout.fin)

        dm0_vs = vs.mean(axis=1) 

        if values.show:
            
            f, (ax1, ax2) = pylab.subplots(1,2)
            ax1.imshow(v.T, aspect='auto')
            ax2.imshow(vs.T, aspect='auto')

            #bins = np.arange(-6, 6, 0.01)
            bins =100
            pylab.figure()

            pylab.hist(v.flatten(), histtype='step', bins=bins)
            pylab.hist(vs.flatten(), histtype='step', bins=bins)

            f, (ax1, ax2) = pylab.subplots(1,2)
            vf = abs(np.fft.rfft(v, axis=0))**2
            vsf = abs(np.fft.rfft(vs, axis=0))**2
            ax1.imshow(np.log10(vf[1:, 1:].T), aspect='auto')
            ax2.imshow(np.log10(vsf[1:, 1:].T), aspect='auto')
            

            pylab.figure()
            pylab.plot(v.std(axis=0))
            pylab.plot(vs.std(axis=0))
    
            pylab.figure()
            pylab.plot(dm0)
            pylab.plot(dm0_vs)
            pylab.figure()
            pylab.plot(abs(np.fft.rfft(dm0))**2)
            pylab.plot(abs(np.fft.rfft(dm0_vs))**2)
            
            
            pylab.show()

    if spout:
        spout.fin.flush()
        spout.fin.close()
    


if __name__ == '__main__':
    _main()
