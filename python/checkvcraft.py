#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import vcraft
from craftcor import print_delay

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = vcraft.VcraftFile(values.files[0])
    ndelay = 32
    nint = 1024*32
    nfft = 64
    nchan = len(f.freqs)
    nchan = 2
    #data = np.zeros((ntimes, f.nchans, ndelay), dtype=np.complex)

    rdata = np.random.randn((nint+ndelay)*nfft) + 1j*np.random.randn((nint+ndelay)*nfft)
    fig, (ax1, ax2) = pylab.subplots(2,1)
    delays = np.zeros((ndelay, nchan))
    phases = np.zeros((ndelay, nchan))
    df0 = np.zeros((nint, nfft, nchan), dtype=np.complex)


    for delay in xrange(ndelay):
        d = f.read(delay, nfft*nint)
        print d.shape
        for c in xrange(nchan):
            dc = d[:, c] # chanel 0
            #dc = rdata[delay:delay+nint*nfft]
            dc.shape = (-1, nfft)
            df = np.fft.fftshift(np.fft.fft(dc, axis=1), axes=1)
            if delay == 0:
                df0[:,:,c] = df

            dfx = (df*np.conj(df0[:,:,c])).mean(axis=0)
            ax1.plot(np.abs(dfx))
            ax2.plot(np.angle(dfx))
            print 'DELAY = {} sample'.format(delay)
            delays[delay, c], phases[delay, c] = print_delay(dfx)

    fig, (ax1, ax2) = pylab.subplots(2,1, sharex=True)
    ax1.plot(delays)
    ax1.set_ylabel('Delays')
    ax2.plot(phases)
    ax2.set_ylabel('Phase (deg)')
    ax2.set_xlabel('Delay (samp)')

    pylab.show()






if __name__ == '__main__':
    _main()
