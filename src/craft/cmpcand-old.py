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

def mergecand(dall, dt, ddm, dm):
    d0 = dall[0]
    d1 = dall[1]
    t0 = d0[:, 7]
    t1 = d1[:, 7]
    dm0 = d0[:, 5]
    dm1 = d1[:, 5]
    dmerge = []

    for it, t in enumerate(t0):
        candidx = np.argmin(np.abs(t - t1))
        if abs(t1[candidx] - t0[it]) < dt and abs(dm0[it] - dm1[candidx]) < ddm:
            dmerge.append((d0[it, :], d1[candidx, :]))

    return np.array(dmerge)
    

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--spfile', help='Sigproc file to overplot')
    parser.add_argument('-d','--dm', type=float, help='DM')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    dall = []
    fig, (ax1, ax2, ax3) = pylab.subplots(3,1)
    mjd0 = None
    
    for f in values.files:
        d = np.loadtxt(f)
        sn = d[:, 0]
        sec = d[:, 2]
        mjd = d[:, 7]
        if mjd0 is None:
            mjd0 = mjd[0]
        print(d.shape)
        secoff = (mjd - mjd0)*86400
        ax1.plot(sec,sn, label=f)
        dall.append(d)

    ax1.legend(frameon=False)
    ax1.set_ylabel('Pulsar Single pulse S/N')

    pylab.show()

    if values.spfile:
        spfile = sigproc.SigprocFile(values.spfile)
        d = np.fromfile(spfile.fin, dtype=np.uint8)
        nchan = spfile.nchans
        nsamp = len(d)/nchan
        d.shape = (nsamp, nchan)
        #pylab.figure()
        #pylab.imshow(d, aspect='auto')
        rfipower = d[:, 250:260].mean(axis=1)
        norfipower = d[:, 50:60].mean(axis=1)
        times = np.arange(len(rfipower))*spfile.tsamp
        ax3.plot(times, rfipower, label='GPS channel')
        ax3.plot(times, norfipower, label='Non-RFI channel')

        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Power')
        ax3.legend(frameon=False, loc='upper left')

    pylab.figure()
    dt = 0.001/86400.
    merged = mergecand(dall,dt, 10, dm=values.dm)
    pylab.scatter(merged[:, 0, 0], merged[:, 1, 0])
    pylab.xlabel('S/N (%s)' % values.files[0])
    pylab.ylabel('S/N (%s)' % values.files[1])
    pylab.plot([merged[:, 0, 0].min(), merged[:, 0,0].max()], [merged[:, 1, 0].min(), merged[:, 1, 0].max()])

    times = merged[:, 0, 2]
    ax2.plot(times, merged[:, 0, 0] - merged[:, 1, 0])
    ax2.set_ylabel('S/N difference'.format(values.files[0], values.files[1]))

    pylab.show()


    
    

if __name__ == '__main__':
    _main()
