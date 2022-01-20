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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose', default='sum.fil')
    parser.add_argument('-o','--output', help='output file', default='sum.fil')
    #parser.add_argument('-c','--channel-average', help='Average this many channels', type=int, default=0)
    parser.add_argument('-n','--nsamps-per-block', help='Number of samples per block for rescaling', type=int, default=1)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    infiles = [sigproc.SigprocFile(f) for f in values.files]
    s0 = infiles[0]
    tstarts = np.array([f.tstart for f in infiles])
    tstart_max = max(tstarts)
    sampoffs = np.round((tstart_max - tstarts)*86400.0/s0.tsamp).astype(int)
    print('Sample offsets', sampoffs)
    nin = len(infiles)
    foff_out = s0.foff 
    # If you average N channels together, the center frequency moves.
    fch1_out = s0.fch1 
    #assert s0.nchans % values.channel_average == 0
    nchan_out = s0.nchans
    
    mjd = s0.tstart

    hdr = {'fch1':fch1_out, 'foff':foff_out,'tsamp':s0.tsamp, 'tstart':tstart_max, 'nbits':32, 'nifs':1, 'nchans':nchan_out, 'src_raj':s0.src_raj, 'src_dej':s0.src_dej}
    fout = sigproc.SigprocFile(values.output, 'w', hdr)
    ntimes = values.nsamps_per_block
    sampno = 0
    nsamples = min([f.nsamples - soff for (f, soff) in zip(infiles, sampoffs)])
    while sampno < nsamples - ntimes:
        d = np.zeros((ntimes, nchan_out), dtype=np.float32)
        for ifile, f in enumerate(infiles):
            sampoff = sampoffs[ifile]
            i = sampno + sampoff
            din = f[i:i + ntimes].T
            (nc, ntimes) = din.shape
            din.shape = (nchan_out, -1, ntimes)
            davg = din.mean(axis=1).T
            d += davg

        d /= float(nin)
        d.tofile(fout.fin)
        sampno += ntimes

    fout.fin.close()

    

if __name__ == '__main__':
    _main()
