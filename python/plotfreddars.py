#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2019
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from rtdata import FreddaRescaleData


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def showplot(ax, d, values, name, x=None):
    if values.image:
        ax.imshow(d, aspect='auto')
    else:
        if x is None:
            ax.plot(d)
        else:
            ax.plot(x, d)

    ax.set_ylabel(name)

def plot(f, values):
    d = FreddaRescaleData(f)
    print d.hdr
    print d.dada_files
    print d.dada_files[0].nblocks
    print d.nblocks
    print d.antennas

    fig, ax = pylab.subplots(3,3)
    ax = ax.flatten()

    fig2, ax2 = pylab.subplots(3,3)
    ax2 = ax2.flatten()

    fig3, ax3 = pylab.subplots(3,3)
    ax3 = ax3.flatten()
    
    blkidx = values.blkidx
    iant = values.iant
    ichan = values.ichan
    ibeam = values.ibeam
    bd = d[blkidx]
    iax = 0

    for iname, name in enumerate(['mean','std','kurt','scale','offset', 'decay_offset', 'nsamps']):

        bdn = bd[name]
        if values.log:
            bdn= 10*np.log10(bdn)
        
        print name, bdn.shape # shape is (ant, beam, chan)
        showplot(ax[iax],  bdn[iant, :, :].T, values, name, x=d.freqs)
        showplot(ax2[iax], bdn[:,:,ichan], values, name)
        showplot(ax3[iax], bdn[:, ibeam, :].T, values, name, x=d.freqs)

        iax += 1

    for iname, name in enumerate(['dm0','dm0count']):
        bdn = bd[name]
        showplot(ax[iax],  bdn[iant, :, :].T, values, name)
        showplot(ax2[iax], bdn[:,:,ichan] , values, name)
        showplot(ax3[iax], bdn[:,ibeam,:].T, values, name)
        iax += 1

    iax -= 1
    fig.suptitle('Antenna = {}={}'.format(iant, d.antennas[iant]))
    
    fig2.suptitle('Channel = {}={} MHz'.format(ichan, d.freqs[ichan]))
    ax2[6].set_xlabel('Beam')
    ax2[6].set_ylabel('Antenna')

    fig3.suptitle('Beam = {}'.format(ibeam))
    ax3[6].set_xlabel('Channel')
    ax3[6].set_ylabel('Antenna')


    pylab.show()

    

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    #parser.add_argument('-b','--beam', type=int, antenna='bea
    parser.add_argument('-a','--iant', type=int, help='Antenna number', default=0)
    parser.add_argument('-c','--ichan', type=int, help='Channel number', default=0)
    parser.add_argument('-b','--ibeam', type=int, help='Beam number', default=0)
    parser.add_argument('-i','--blkidx', type=int, help='Block index', default=0)
    parser.add_argument('-l','--log', action='store_true', default=False, help='do log on imshow')
    parser.add_argument('--image', action='store_true', default=False, help='Show images rathe rthan lines')
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
