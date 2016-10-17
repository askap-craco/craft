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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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

    
    nf = 128
    nt = 128
    fmin = 1.200
    fmax = 1.600
    freq = np.linspace(fmin, fmax, nf)
    dm = 30.
    assert(len(freq) == nf)
    tint = 0.1
    #d = np.ones((nf, nt), dtype=np.float32)
    d = np.zeros((nf, nt), dtype=np.float32)
    d += np.random.randn(nf, nt)

    for iff, f in enumerate(freq):
        dt = int(np.round(4.15*dm*(f**-2 - fmax**-2)))
        d[iff, dt] += 10.
    
    with open('test.in', 'w') as fout:
        d.tofile(fout)
    
    try:
        subprocess.check_call('fdmt test.in test.out'.split())
    except:
        print 'PRoblem processing, will try anyway'

    with open('test.out', 'r') as fin:
        dout = np.fromfile(fin, dtype=np.float32, count=nf*nt)
        dout.shape = (nf, nt)

    fig, ax  = pylab.subplots(1,2)
    ex = (0, nt*tint, fmin, fmax)
    ex = None
    ax[0].imshow(d, origin='lower', extent=ex)
    ax[1].imshow(dout, origin='lower', extent=ex)


    for i in xrange(7):
        pylab.figure()
        fname = 'state_s%d.dat' % i
        fin = open(fname, 'r')
        theshape =  np.fromfile(fin, dtype=np.int32, count=3)
        d = np.fromfile(fin, dtype=np.float32, count=theshape.prod())
        d.shape = theshape
        pylab.imshow(d[:, 0, :], aspect='auto', origin='lower')
        print fname, d.shape
        fin.close()

    pylab.show()
    



if __name__ == '__main__':
    _main()
