#!/usr/bin/env python
"""
Cleans multibeam filterbanks using eigenflagging/subtraction.

See Kocz 2012 and Wax and Kailath 1985

Copyright (C) CSIRO 2018
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples per block to calculate covariance over and subtract', default=256)
    parser.add_argument('-k', '--nsub', type=int, default=1, help='Fixed Number of eigenvalues to subtract (todo: Dynamic)')
    parser.add_argument('-D','--outdir', required=True, help='Output directory for cleaned filterbanks')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    insfs = [sigproc.SigprocFile(f) for f in values.files]
    outfile_names = [os.path.join(values.outdir, os.path.basename(f)) for f in values.files]
    outsfs = [sigproc.SigprocFile(f, mode='w', header=s.header) for f, s in zip(outfile_names, insfs)]
    endsamp = min([s.nsamples for s in insfs])
    nsamp = values.nsamp
    nblock = endsamp // nsamp
    k = values.nsub

    for b in xrange(nblock):
        d = np.array([s[b*nsamp:(b+1)*nsamp] for s in insfs])
        nbeam, nint, nchan = d.shape
        dout = np.empty_like(d)
        for ic in xrange(nchan):
            r = d[:, :, ic]
            thecov = np.cov(r)
            #dcovar[ic, :, :] = thecov
            u, v, vh = np.linalg.svd(thecov)
            for t in xrange(nint):
                rfi = np.dot(u.T, r[:, t])
                #rfiout[:, t, ic] = rfi[:]
                rfi[:k] = 0
                dout[:, t, ic] = np.dot(rfi, u.T)
                
        for iout, sout in enumerate(outsfs):
            dout[iout, :, :].astype(d.dtype).flatten().tofile(sout.fin)

    for s in outsfs:
        s.close()

if __name__ == '__main__':
    _main()
