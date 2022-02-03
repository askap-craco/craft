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
from . import plotutil
import fnmatch

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False, show=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for fin in values.files:
        dechirp_file(fin, values)

    if values.show:
        pylab.show()



def dechirp_file(fin, values):
    vin = np.loadtxt(fin)
    sn = vin[:, 0]
    sampno = vin[:, 1]
    time = vin[:, 2]
    boxcar = vin[:, 3]
    dm = vin[:, 4]
    beamno = vin[:, 5]

    ubeams = set(beamno)
    fig, axs = plotutil.subplots(1,1)
    for b in sorted(ubeams):
        bmask = beamno == b
        axs[0].plot(time[bmask], dm[bmask]+1, marker='x', ls='None', label='Beam %d' %b, picker=3)


    return vin

if __name__ == '__main__':
    _main()
