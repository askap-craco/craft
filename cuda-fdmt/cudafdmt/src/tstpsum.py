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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    NBOX = 32
    history = np.arange(NBOX)
    n = NBOX/2
    offset = 1

    state = 0
    outhist = history.copy()
    for t in xrange(NBOX-1,-1,-1):
        outhist[t] += state
        state = outhist[t]

    print 'Reverse prefix sum', history, outhist
        

    print 'Init history', history
    while n != 0:
        for t in xrange(NBOX):
            if t < n:
                oidx = t*2
                inidx = oidx + 1
                print 'n={} off={} t={} oidx={} inidx={} his[oidx]={} hist[inidx]={}'.format(n, offset, t, oidx, inidx, history[oidx], history[inidx])
                history[oidx] += history[inidx]


        n /= 2
        offset *= 2
        print 'History', n, offset, history
        

    

if __name__ == '__main__':
    _main()
