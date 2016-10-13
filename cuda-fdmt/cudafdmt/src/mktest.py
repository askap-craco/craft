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
    kwargs['interpolation'] = 'None'
    im = ax.imshow(*args, **kwargs)
    ax.format_coord = Formatter(im)
    return im

def stft(X, N_f, N_t, N_bins):
    ''' Copied from FDMT.py'''
    XX = np.abs(np.fft.fft(X.reshape([N_f,N_t*N_bins]),axis = 1))**2
    XX = np.transpose(XX)
    XX = np.sum(XX.reshape(N_f,N_bins,N_t),axis=1)
    
    E = np.mean(XX[:,:10])
    XX -= E
    V = np.var(XX[:,:10])
    
    XX /= (0.25*np.sqrt(V))
    V = np.var(XX[:,:10])

    return XX

def make_signal(f_min = 1200., f_max = 1600, N_bins = 40, N_t = 512, N_f = 512, D=5):
    ''' COpied from FDMT.py'''
    N_total = N_f*N_t*N_bins
    PulseLength = N_f*N_bins
    PulseSig = 0.4
    PulsePosition = 4134567/4
    maxDT = N_t
    
    PDB("Signal preparation")
    practicalD = DispersionConstant * D
    I = np.random.normal(0,1,N_total)
    I[PulsePosition:PulsePosition+PulseLength] += np.random.normal(0,PulseSig,PulseLength)*10
    print "MAX Thoretical SNR:", np.sum(np.abs(I[PulsePosition:PulsePosition+PulseLength])**2 - np.mean(abs(I)**2)) / (np.sqrt(PulseLength*np.var(abs(I)**2)))
    
    X = CoherentDedispersion(I, -D, f_min,f_max,False)    

    return stft(X, N_f, N_t, N_bins)


def make_signal2():
    freq = np.linspace(fmin, fmax, nf)
    dm = 150
    assert(len(freq) == nf)
    #d = np.ones((nf, nt), dtype=np.float32)
    d = np.zeros((nf, nt), dtype=np.float32)
    d += np.random.randn(nf, nt)
    return d


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    tint = 0.1
    nf = 512
    nt_sim = 512
    nt = 256
    tstart = nt
    fmax = 1440.
    fmin = fmax - float(nf)

    if os.path.exists('test.in'):
        d = np.fromfile('test.in', dtype=np.float32, count=nf*nt)
        d.shape = (nf, nt)
    else:
        with open('test.in', 'w') as fout:
            d = make_signal(fmin, fmax, 40, nt_sim, nf, D=5)
            d = d[:, tstart:tstart+nt]
            d.astype(np.float32).tofile(fout)

    print 'signal shape', d.shape

    nf, nt = d.shape

    with open('test.out', 'r') as fin:
        dout = np.fromfile(fin, dtype=np.float32, count=nf*nt)
        dout.shape = (nf, nt)

    fig, ax  = pylab.subplots(1,2)
    ex = (0, nt*tint, fmin, fmax)
    ex = None
    myimshow(ax[0], d, origin='lower', extent=ex)
    myimshow(ax[1], dout, origin='lower', extent=ex)
    ax[0].set_title('test.out')
    print 'Dout max', dout.max()


    i = 0
    while True:
        fig, ax = pylab.subplots(1,2)
        fname = 'state_s%d.dat' % i
        if not os.path.exists(fname):
            break

        fin = open(fname, 'r')
        theshape =  np.fromfile(fin, dtype=np.int32, count=4)
        print theshape
        
        d = np.fromfile(fin, dtype=np.float32, count=theshape.prod())
        d.shape = theshape
        beam = 0
        myimshow(ax[0], d[beam, :, 0, :], aspect='auto', origin='lower')
        myimshow(ax[1], d[beam, 0, :, :], aspect='auto', origin='lower')

        ax[0].set_title(fname)

        print fname, d.shape, np.prod(d.shape)
        fin.close()
        i += 1

    pylab.show()
    



if __name__ == '__main__':
    _main()
