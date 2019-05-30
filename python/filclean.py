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


def eigenflag(d, values):
    nbeam, nint, nchan = d.shape
    dout = np.empty_like(d)
    vic = np.empty((nchan, nbeam))
    k = values.nsub
    for ic in xrange(nchan):
        r = d[:, :, ic]
        thecov = np.cov(r)
        #dcovar[ic, :, :] = thecov
        u, v, vh = np.linalg.svd(thecov)
        vic[ic, :] = v
        for t in xrange(nint):
            rfi = np.dot(u.T, r[:, t])
            #rfiout[:, t, ic] = rfi[:]
            rfi[:k] = 0
            dout[:, t, ic] = np.dot(rfi, u.T)


    if values.show:
        fig, ax = pylab.subplots(2,1)
        ax = ax.flatten()
        ax[0].plot(vic)
        ax[0].set_ylabel('Singular value')
        ax[0].set_xlabel('Channel number')
        pylab.show()


    return dout

def subtract_dm0(d, values):
    dm0 = d.mean(axis=2)
    dout = d - dm0[:,:,np.newaxis]

    return dout

def eigenflag_dm0(d, values):
    nbeam, nint, nchan = d.shape
    dout = np.empty_like(d)
    dm0 = d.mean(axis=2)
    thecov = np.cov(dm0)
    u, v, vh = np.linalg.svd(thecov)
    k = values.nsub
    opt_dm0 = np.empty((nbeam, nint))

    for t in xrange(nint):
        r = dm0[:, t]
        rfi = np.dot(u.T, r)
        rfi[k:] = 0
        rfixed = np.dot(rfi, u.T)
        opt_dm0[:, t] = rfixed

    dout = d - opt_dm0[:,:,np.newaxis]
        
    if values.show:
        # find most correlated beams - peak off diagonal element
        maxidx = np.argmax(np.triu(thecov, k=1))
        maxbeams = np.unravel_index(maxidx, thecov.shape)

        print maxbeams
        
        fig, ax = pylab.subplots(2,5)
        #ax[0].imshow(dm0.T, aspect='auto')
        #ax[1].imshow(opt_dm0.T, aspect='auto')
        #ax[0,0].plot(dm0.T)
        #ax[0,1].plot(opt_dm0.T)
        ax[0,0].plot(dm0.T[:, maxbeams])
        ax[0,1].plot(opt_dm0.T[:, maxbeams])

        ax[0,2].plot(v)
        ax[0,3].imshow(thecov)
        ax[1,0].imshow(d[9,:,:].T, aspect='auto')
        ax[1,1].imshow(dout[9,:,:].T, aspect='auto')
        ax[1,2].imshow(u)
        ax[1,3].plot(u[:,0])# first eigenvector

        ax[0,4].scatter(dm0[maxbeams[0],:], dm0[maxbeams[1],:])
        ax[1,4].plot(dm0.std(axis=1), 'o')
        ax[1,4].plot(opt_dm0.std(axis=1), 'o')
        pylab.show()

    return dout
        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples per block to calculate covariance over and subtract', default=256)
    parser.add_argument('-k', '--nsub', type=int, default=1, help='Fixed Number of eigenvalues to subtract (todo: Dynamic)')
    parser.add_argument('-D','--outdir', required=True, help='Output directory for cleaned filterbanks')

    parser.add_argument('-e','--eigenflag', action='store_true', help='Apply eigenflagging')
    parser.add_argument('-d','--subtract-dm0', action='store_true', help='Subtract dm0')
    parser.add_argument('-f','--eigenflag-dm0', action='store_true', help='Eigenflag on dm0')
    parser.add_argument('--start-samp', type=int, help='Sample to start on', default=0)
    parser.add_argument('--nblock', type=int, help='Number of blocks to process')
    parser.add_argument('--show', action='store_true', help='Show plots')
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
    startblock = values.start_samp // nsamp
    nblock = values.nblock
    if nblock is None:
        nblock = endsamp // nsamp

    k = values.nsub

    for b in xrange(startblock, nblock):
        d = np.array([s[b*nsamp:(b+1)*nsamp] for s in insfs])
        dtype = d.dtype

        if values.eigenflag:
            d = eigenflag(d, values)
 
        if values.eigenflag_dm0:
            d = eigenflag_dm0(d, values)

        if values.subtract_dm0:
            d = subtract_dm0(d, values)

        for iout, sout in enumerate(outsfs):
            d[iout, :, :].astype(dtype).flatten().tofile(sout.fin)

    for s in outsfs:
        s.fin.close()

if __name__ == '__main__':
    _main()
