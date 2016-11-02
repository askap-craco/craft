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
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"



def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Makes a filterbank with a dispersed pulse in it')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-f', '--fch1', type=float, help='First channel frequency MHz')
    parser.add_argument('-c', '--foff', type=float, help='Frequency offset MHz')
    parser.add_argument('-n', '--nchan', type=int, help='Number of channels')
    parser.add_argument('-b', '--nbins', type=int, help='Number of bins to integrate')
    parser.add_argument('-D', '--dispersion', type=float, help='Dispersion Measure')
    parser.add_argument('-t','--ntimes', type=int, help='Number of samples')
    parser.add_argument('-s','--pulsesig', type=float, help='Pulse S/N')
    parser.add_argument(dest='file', help='Output sigproc filterbank file name')

    parser.set_defaults(verbose=False, fch1=1440., foff=-1., nchan=256, nbins=40, dispersion=5.0, pulsesig=40.0, ntimes=2048)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    nf = values.nchan
    nt_sim = values.ntimes
    tstart = 0
    fmax = values.fch1
    fmin = fmax + values.nchan * values.foff
    assert fmin < fmax
    tint = values.nbins/((fmax - fmin)*1e6)
    hdr = {'fch1': fmax, 'foff':values.foff, 'nchans': values.nchan, 'nifs':1,'tstart':55000.0, 'tsamp':tint, 'nbits':8,
           'data_type':1}
    spfile = sigproc.SigprocFile(values.file, 'w', hdr)

    np.random.seed(42)

    d = dispersed_stft(fmin, fmax, nt_sim, nf, N_bins=values.nbins, D=values.dispersion, PulseSig=values.pulsesig)
    #d = d[:, tstart:tstart+nt]
    d *= 2
    d += 128
    pylab.imshow(d)
    pylab.show()
    d[0:nf, :] = d[nf::-1, :]

    print 'D mean=', d.flatten().mean(), 'D stdev', d.flatten().std(), 'shape', d.shape
    d.astype(np.uint8).tofile(spfile.fin)




if __name__ == '__main__':
    _main()
