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
from plot_fredda import load4d, file_series, comma_list, Formatter

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--beam', type=float, help='beam number')
    parser.add_argument('-d','--dmrange', type=comma_list, help='Dm range to show in time series plot')
    parser.add_argument('-s','--start', type=int, help='Start block')
    parser.add_argument('-n','--maxn', type=int, help='max number of blocks ot plot')
    parser.set_defaults(verbose=False, beam=0, start=4, maxn=10)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for ifname, fname in enumerate(file_series('boxcar_e%d.dat', values.start)):
        alld = load4d(fname)
        
        d = alld[0, 97, :, :]
        fig, axes = pylab.subplots(1,3)
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
                    
        pylab.show()
        
  

if __name__ == '__main__':
    _main()
