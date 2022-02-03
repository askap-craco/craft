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
from .craftobs import load_beams
from .plotutil import subplots, addpick
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def annotate(fig, title, xlabel, ylabel):
    fig.text( 0.5, 0.98,title, ha='center', va='top')
    fig.text(0.5, 0.02, xlabel, ha='center', va='bottom')
    fig.text(0.02, 0.5, ylabel, rotation=90, ha='center', va='top')
    fout='{}.png'.format(title.replace(' ', ''))
    print('Saving', fout)
    fig.savefig(fout)
    

def commasep(s):
    return list(map(int, s.split(',')))

def floatcommasep(s):
    return list(map(float, s.split(',')))


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-t', '--times', help='Integration range to plot', type=commasep)
    parser.add_argument('--nxy', help='number of rows,columns in plots', type=commasep)
    parser.add_argument('--imzrange', help='Z range for dynamic spectrum', type=floatcommasep)
    parser.add_argument('-c','--channels', help='Channel number', type=commasep)
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

    
    nbeams = len(values.files)
    sig0 = sigproc.SigprocFile(values.files[0])
    assert sig0.nifs == 1
    nfreqs = sig0.nchans
    nsamps = sig0.nsamples - tstart
    all_data = np.zeros((nbeams, len(values.channels), nsamps))
    mean_data = np.zeros((nbeams, len(values.channels), nsamps/ntimes))
    std_data = np.zeros((nbeams, len(values.channels), nsamps/ntimes))
    files = [sigproc.SigprocFile(fname) for fname in values.files]
    
    print('Memory', all_data.shape, mean_data.shape)

    for ifname, f in enumerate(files):
        print('Loading file', f.filename)
        nelements = ntimes*f.nifs*f.nchans
        f.seek_data(f.bytes_per_element*tstart)
        t = 0
        if (f.nbits == 8):
            dtype = np.uint8
        elif (f.nbits == 32):
            dtype = np.float32

        tsum = 0
        while t < nsamps:
            print(f.filename, 'reading', nelements)
            v = np.fromfile(f.fin, dtype=dtype, count=nelements )
            if len(v) != nelements:
                break

            v.shape = (ntimes, f.nifs, f.nchans)
            assert f.nifs == 1
            d = v[:, 0, values.channels].T
            all_data[ifname, :, t:t+ntimes] = d
            print(d.shape, d.mean().shape, d.mean())
            mean_data[ifname, :, tsum] = d.mean(axis=1)
            std_data[ifname, :, tsum] = d.std(axis=1)
            t += ntimes
            tsum += 1


    print('Loaded', all_data.shape, mean_data.shape)
    
    f0 = files[0]
    mjdstart = f0.tstart
    tsamp = f0.tsamp
    fch1 = f0.fch1
    foff = f0.foff
    src_raj = f0.src_raj
    src_dej = f0.src_dej

    for c in range(len(values.channels)):
        fig = pylab.figure()
        for b in range(nbeams):
            pylab.plot(mean_data[b, c, :].T, label=files[b].filename, picker=3)
            pylab.xlabel('Sample')
            pylab.ylabel('Amplitude')

        addpick(fig)

    fig, axes = pylab.subplots(*values.nxy, sharex=True, sharey=True)
    axes = axes.flatten()
    for b in range(len(axes)):
        ax = axes[b]
        ax.plot(mean_data[b, :, :].T)
        ax.set_title(files[b].filename[-6:-4])
        

    pylab.show()

if __name__ == '__main__':
    _main()
