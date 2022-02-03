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
import matplotlib.gridspec as gridspec
from .plot_fredda import load4d, file_series, comma_list, Formatter

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--beam', type=int, help='beam number')
    parser.add_argument('-d','--dm', type=int, help='iDm to show in time series plot')
    parser.add_argument('-t','--time', type=int, help='Time sample to print stats of')
    parser.add_argument('-s','--start', type=int, help='Start block')
    parser.add_argument('-w','--boxcar', type=int, help='Boxcar to print stats of')
    parser.add_argument('-n','--maxn', type=int, help='max number of blocks ot plot')
    parser.add_argument('--tsamp', type=float, help='Sample duration for FFT plot (milliseconds)', default=1.0)
    parser.set_defaults(verbose=False, beam=0, start=4, maxn=10, dm=0, boxcar=0, time=0)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for ifname, fname in enumerate(file_series('boxcar_e%d.dat', values.start)):
        alld = load4d(fname)
        
        # shape = beam, idm, t, boxcar
        d = alld[values.beam, values.dm, :, :]
        fig, axes = pylab.subplots(2,3)
        ax = axes.flatten()
        ax[0].plot(d)
        ax[0].set_xlabel('Sample')
        ax[0].set_ylabel('S/N')

        ax[1].plot(d.std(axis=0))
        ax[1].set_xlabel('Boxcar')
        ax[1].set_ylabel('Std')
        
        ax2 = ax[2]
        ax2.plot(d.mean(axis=0), 'r')
        ax2.set_ylabel('mean')
        ax2.set_xlabel('Boxcar')

        print('S/N per beam idm=%d t=%d bc=%d: %s' % (values.dm, values.time, values.boxcar, alld[:, values.dm, values.time, values.boxcar]))

        ax3 = ax[3]
        ax3.imshow(alld[values.beam, :, :, values.boxcar], aspect='auto', origin='lower')
        ax[3].set_xlabel('Sample')
        ax[3].set_ylabel('idt')

        ax4 = ax[4]
        for w in (0,1,30,31):
            ax4.hist(alld[values.beam,:,:,w].flatten(), bins=100, histtype='step', label='w={}'.format(w))
            
        ax4.legend()
        ax4.set_xlabel('S/N')
        ax4.set_ylabel('Number of points')

        ax5 = ax[5]
        x = alld[values.beam,values.dm,:,values.boxcar]
        xamp = abs(np.fft.rfft(x))**2
        tsamp = values.tsamp
        fs = 1.0/tsamp
        fsfreq = np.arange(len(xamp))*fs/2.0
        ax5.loglog(fsfreq, xamp, label='b={} idt={} w={}'.format(values.beam, values.dm, values.boxcar))
        ax5.set_ylabel('FFT^2')
        ax5.set_xlabel('Hz')
                    
        pylab.show()
        
  

if __name__ == '__main__':
    _main()
