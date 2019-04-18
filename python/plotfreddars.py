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


def plot(f, values):
    d = FreddaRescaleData(f)
    print d.hdr
    print d.dada_files
    print d.dada_files[0].nblocks
    print d.nblocks

    fig, ax = pylab.subplots(3,3)
    ax = ax.flatten()
    blkidx = 4
    iant = 1
    bd = d[blkidx]
    iax = 0

    for iname, name in enumerate(['mean','std','kurt','scale','offset', 'decay_offset', 'nsamps']):

        bdn = bd[name][iant, :, :]
        print name, bdn.shape
        ax[iax].plot(d.freqs, bdn.T)
        ax[iax].set_title(name)
        iax += 1

    axstart = iname + 1
    for iname, name in enumerate(['dm0','dm0count']):
        bdn = bd[name][iant, :, :]
        print name, bdn.shape
        ax[iax].plot(bdn.T)
        ax[iax].set_title(name)
        iax += 1


    pylab.show()

    

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

    for f in values.files:
        plot(f, values)
    

if __name__ == '__main__':
    _main()
