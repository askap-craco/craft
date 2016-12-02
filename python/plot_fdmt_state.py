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
import subprocess
from craftsim import dispersed_voltage, dispersed_stft
from FDMT import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}'.format(x, y, z)

def myimshow(ax, *args, **kwargs):
    kwargs['aspect'] = 'auto'
    kwargs['interpolation'] = 'nearest'
    im = ax.imshow(*args, **kwargs)
    ax.format_coord = Formatter(im)
    return im

def make_signal2():
    freq = np.linspace(fmin, fmax, nf)
    dm = 150
    assert(len(freq) == nf)
    #d = np.ones((nf, nt), dtype=np.float32)
    d = np.zeros((nf, nt), dtype=np.float32)
    d += np.random.randn(nf, nt)
    return d


    N_total = N_f*N_t*N_bins
    PulseLength = N_f*N_bins

def load4d(fname, dtype=np.float32):
    fin = open(fname, 'r')
    theshape =  np.fromfile(fin, dtype=np.int32, count=4)
    d = np.fromfile(fin, dtype=dtype, count=theshape.prod())
    d.shape = theshape
    fin.close()
    print 'load4d', fname, d.shape
    return d

def file_series(prefix):
    i = 0
    while True:
        fname = prefix % i
        if os.path.exists(fname):
            yield fname
        else:
            break
        i += 1

def show_series(prefix, theslice):
    for fname in file_series(prefix):
        ostate = load4d(fname)
        pylab.figure()
        v = ostate[theslice]
        print fname, ostate.shape, 'zeros?', np.all(ostate == 0), 'max', v.max(), np.unravel_index(v.argmax(), v.shape)
        myimshow(pylab.gca(), v, aspect='auto', origin='lower')
        pylab.title(fname)

def plot_series(prefix, theslice):
    for fname in file_series(prefix):
        ostate = load4d(fname)
        pylab.figure()
        v = ostate[theslice]
        print fname, ostate.shape, 'zeros?', np.all(ostate == 0), 'max', v.max(), np.unravel_index(v.argmax(), v.shape)
        pylab.plot(v.T)
        pylab.title(fname)



def plot_stats():
    for fname in file_series('mean_e%d.dat'):
        fig, ax = pylab.subplots(2,2)
        ax = ax.flatten()
        bslice = [0,0, slice(None), slice(None)]
        bmean = load4d(fname)[bslice].T
        bstd = load4d(fname.replace('mean','std'))[bslice].T
        bkur = load4d(fname.replace('mean','kurt'))[bslice].T
        bdm0 = load4d(fname.replace('mean','dm0'))[bslice].T
        ax[0].plot(bmean)
        ax[0].set_ylabel('Mean')
        ax[1].plot(bstd)
        ax[1].set_ylabel('Std')
        ax[2].plot(bkur)
        ax[2].set_ylabel('Kurtosis')
        ax[3].plot(bdm0)
        ax[3].set_ylabel('Dm0')
        fig.text(0.5, 0.98, fname.replace('mean_','').replace('.dat',''), ha='center', va='top')
        print fname


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--beam', type=float, help='beam number')
    parser.set_defaults(verbose=False, beam=0)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    beam = values.beam

    #show_series('ostate_e%d.dat', [beam, 0, slice(None), slice(None)])
    #show_series('finalstate_e%d.dat', [beam, 0, slice(None), slice(None)])
    #show_series('initstate_e%d.dat',[beam, slice(None),0, slice(None)])
    #plot_series('mean_e%d.dat',[0,0, slice(None), slice(None)])
    #plot_series('std_e%d.dat',[0,0, slice(None), slice(None)])
    #plot_series('kurt_e%d.dat',[0,0, slice(None), slice(None)])
    #plot_series('dm0_e%d.dat',[0,0, slice(None), slice(None)])
    plot_stats()
    #plot_series('nsamps_e%d.dat',[0,0, slice(None), slice(None)])


    for i, fname in enumerate(file_series('state_s%d.dat')):
        if 0 < i < 9:
            continue

        d = load4d(fname)
        fig, ax = pylab.subplots(1,2)
        myimshow(ax[0], d[beam, :, 0, :], aspect='auto', origin='lower')
        myimshow(ax[1], d[beam, 0, :, :], aspect='auto', origin='lower')

        ax[0].set_title(fname)

        print fname, d.shape, np.prod(d.shape)

    pylab.show()
    



if __name__ == '__main__':
    _main()
