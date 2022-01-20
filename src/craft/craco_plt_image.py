#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def plot(f, values):
    nd = values.nd
    nt = values.nt
    npix = values.npix
    xpix, ypix = list(map(int, values.pix.split(',')))
    
    d = np.fromfile(f, dtype=np.complex64).reshape(nd, nt/2, npix, npix)
    imgs = np.zeros((nd, nt, npix, npix))
    imgs[:, 0::2, :, :] = d.real
    imgs[:, 1::2, :, :] = d.imag
    maxidx = np.argmax(imgs)
    dmax, tmax, ymax, xmax = np.unravel_index(maxidx, imgs.shape)
    print('max at', maxidx, dmax, tmax, ymax, xmax, imgs.flat[maxidx])
    fig,ax = pylab.subplots(1,2)
    ax[0].imshow(imgs[dmax, tmax, :, :], aspect='auto', origin='lower')
    ax[1].plot(imgs[:, :, ymax, xmax].T)
    ax[1].set_xlabel('Time')
    ax[1].set_ylabel('Amplitude')

    
    pylab.show()

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--nd', type=int, default=4)
    parser.add_argument('--nt', type=int, default=16)
    parser.add_argument('--npix', type=int, default=256)
    parser.add_argument('--pix', default='128,128')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        plot(f, values)
    

if __name__ == '__main__':
    _main()
