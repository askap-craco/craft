#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2018
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import vcraft

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-n','--nsamp', type=int, default=1024, help='Number of samples to check')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    muxes = vcraft.mux_by_antenna(values.files)
    fix, allax = pylab.subplots(4,5, sharex=True, sharey=True)
    allax = allax.flatten()
    for mux, ax in zip(muxes, allax):
        samps = mux.read(0, values.nsamp)
        ax.set_title(mux.ant)
        print(samps.shape)
        m = samps.real.mean(axis=0)
        s = samps.real.std(axis=0)
        print(m.shape, s.shape, mux.freqs.shape)
        ax.plot(mux.freqs, m, 'o')
        ax2 = ax.twinx()
        ax.plot(mux.freqs, s, 'x')
        

    pylab.show()

    

if __name__ == '__main__':
    _main()
