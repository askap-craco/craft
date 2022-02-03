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
import pandas as pd

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def loadcand(fin):
    print('Loading', fin)
    icscand = np.loadtxt(fin, usecols=(0,1,2,3,4,5,6,7))
    icsmask = (icscand[:,0] >= 7) & (65 < icscand[:,5]) & (icscand[:, 5] < 70) 
    icscand = icscand[icsmask, :]
    idxs = np.argsort(icscand[:, 1]) # sort bty sample number
    icscand = icscand[idxs, :]
    return icscand


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d','--dm', type=float, help='DM')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    mjd0 = None
    for f in values.files:
        d = loadcand(f)
        print(d.shape)
        sn = d[:, 0]
        sampno = d[:, 1]
        secs = d[:, 2]
        mjd = d[:, 7]
        if mjd0 is None:
            mjd0 = mjd[0]

        secoff = (mjd - mjd0)*86400.0
        pylab.plot(secoff, sn, label=f)

    pylab.legend()
    pylab.show()

        
    

if __name__ == '__main__':
    _main()
