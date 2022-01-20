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
from . import crafthdr
import itertools

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def load_headers(hdr_files):
    hdrs = [crafthdr.DadaHeader.fromfile(f) for f in hdr_files]
    hdrs.sort(key=lambda h: (int(h.get_value('ANTENNA_NO'), int(h.get_value('CARD_NO')), int(h.get_value('FPGA_ID')))))
    group_hdrs = itertools.groupby(hdrs, key=lambda h:h.get_value('ANT'))
    ant_hdrs = {ant:hdrs for ant, hdrs in group_hdrs}
    print(list(ant_hdrs.keys()))
    return ant_hdrs

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', default=False, help='Show plots')
    parser.add_argument('-n','--nant', help='Number of antennas', type=int)
    parser.add_argument('-c','--nchan', help='Number of channels per MHz', type=int, default=54)
    parser.add_argument('--headers', help='Header files', nargs='*')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    nant = values.nant
    nchan = values.nchan

    for f in values.files:
        rfile = f
        imagfile = f.replace('.real','.imag')
        rdata = np.loadtxt(rfile)
        idata = np.loadtxt(imagfile)
        assert rdata.shape == idata.shape
        xaxis_all = rdata[:, 0]
        xaxis_all.shape = (nant, -1) # should be the same for all antennas
        for a in range(nant):
            xaxis = xaxis_all[a, :]
            assert np.all(xaxis == xaxis_all[0, :])


        
        assert np.all(rdata[:, 0] == idata[:, 0])
        fulld = np.zeros(len(rdata), dtype=np.complex64)
        fulld.real = rdata[:, 1]
        fulld.imag = idata[:, 1]
        #fulld.shape = (nant, -1)
        #pylab.plot(xaxis, np.angle(fulld.T))
        print('FULLDSHAPE', fulld.shape)
        fulld.shape = (nant, -1, nchan)
        delayspec = np.fft.fftshift(np.fft.fft(fulld, axis=2), axes=2)
        print('DElayspec', delayspec.shape, delayspec.dtype)
        ncol = int(np.ceil(np.sqrt(nant)))
        nrow = (nant)//ncol + 1
        print('NROW, NCOL', nrow, ncol)

        fig, ax = pylab.subplots(nrow, ncol)
        ax = ax.flatten()
        #fig2, ax2 = pylab.subplots(2,3)
        
        #ax2 = ax2.flatten()
        xaxis_reshape = xaxis.reshape(-1, nchan).mean(axis=1)
        extent = (xaxis[0], xaxis[-1], -nchan/2 + 0.5, nchan/2 + 0.5)
        for a in range(nant):
            ax[a].imshow(np.log10(abs(delayspec[a, :, :].T)), aspect='auto', extent=extent)
            if a == 0:
                ax[a].set_ylabel('delay (samples)')
            if a == 3:
                ax[a].set_xlabel('Frequency (GHz)')

        rootfile = f.replace('.real','')
        fig.suptitle(rootfile)
        fig.savefig(rootfile+'.delay.png')
        if values.show:
            pylab.show()

    subprocess.check_call('convert -delay 500 -loop 0 *.delay.png delays.gif', shell=True)
        
            
        
        
    

if __name__ == '__main__':
    _main()
