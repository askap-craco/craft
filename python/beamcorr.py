#!/usr/bin/env python
"""
Try out Wael's beamn crrelation

Copyright (C) CSIRO 2019
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
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    sfs = [sigproc.SigprocFile(f) for f in values.files]
    startsamp = 3176800
    nsamp = 256
    d = np.array([s[startsamp:startsamp+nsamp] for s in sfs])
    print d.shape
    nbeam, nint, nchan = d.shape
    dcovar = np.empty((nchan, nbeam, nbeam))
    dsub = np.empty_like(d)
    svd_u = []
    svd_v = []
    svd_vh = []
    k = 0
    
    for c in xrange(nchan):
        thecov = np.cov(d[:, :, c])
        dcovar[c, :, :] = thecov
        u, v, vh = np.linalg.svd(thecov)
        svd_u.append(u)
        svd_v.append(v)
        svd_vh.append(vh)
        for  t in xrange(nint):
            print d[:, t, c].shape, v.shape
            Y = np.dot(d[:,t,c], .T)
            Y[:k] = 0
            dsub[:, t, c] = np.dot(Y, v.T)


    svd_u = np.array(svd_u)
    svd_v = np.array(svd_v)
    svd_vh = np.array(svd_vh)

    print 'u=',svd_u.shape, 'v=',svd_v.shape, 'vh=',svd_vh.shape
    f1 = sfs[0].fch1
    foff = sfs[0].foff
    fend = f1 + nchan*foff
    freqs = np.arange(nchan)*foff + f1

    extent = (0, nbeam*nbeam, fend, f1)# left right bottom top
    pylab.imshow(dcovar.reshape(nchan, -1), extent=extent)
    pylab.ylabel('Frequency')
    pylab.xlabel('Covariance matrix entry')

    pylab.figure()
    pylab.plot(freqs, svd_v)
    pylab.xlabel('Frequency')
    pylab.title('SVD V')

    pylab.figure()
    pylab.imshow(svd_u.reshape(nchan, -1))
    pylab.title('SVD U')

    pylab.figure()
    pylab.imshow(svd_vh.reshape(nchan, -1))
    pylab.title('SVD VH')

    fig, ax = pylab.subplots(1,2)
    b = 9
    ax[0].imshow(d[b,:,:])
    ax[1].imshow(dsub[b, :, :])

    pylab.show()

    

if __name__ == '__main__':
    _main()
