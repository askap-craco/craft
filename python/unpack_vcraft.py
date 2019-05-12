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
import numba

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def unpack_mode2(dwords, d):
    for samp in xrange(4):
        # convert from 4-bit two's complement to int8
        d[samp::4, :, 0] = (dwords & 0xf) - (dwords & 0x8)*2 # real
        dwords >>= 4
        d[samp::4, :, 1] = (dwords & 0xf) - (dwords & 0x8)*2 # imag
        dwords >>= 4
        
    return d
    
    
@numba.jit(nopython=True)
def unpack_mode2_jit(dwords, d):
    for samp in xrange(4):
        # convert from 4-bit two's complement to int8
        d[samp::4, :, 0] = (dwords & 0xf) - (dwords & 0x8)*2 # real
        dwords >>= 4 # numba doesn't support this - so it fails.
        d[samp::4, :, 1] = (dwords & 0xf) - (dwords & 0x8)*2 # imag
        dwords >>= 4
        
    return d 

@numba.jit(nopython=True)
def unpack_mode2_jit_v2(dwords, d):
    (nsampwords, nchan) = dwords.shape
    for s in xrange(nsampwords):
        for c in xrange(nchan):
            word = dwords[s, c]
            for samp in xrange(4):
                d[4*s + samp, c, 0] = (word & 0xf) - (word & 0x8)*2 # real
                word >>= 4
                d[4*s + samp, c, 1] = (word & 0xf) - (word & 0x8)*2 # imag
                word >>= 4
    return d


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
