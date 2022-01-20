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
from . import crafthdr

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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

    for f in values.files:
        hdr = crafthdr.DadaHeader.fromfile(f)
        sz = int(hdr.get_value('HDR_SIZE'))
        with open(f, 'r') as fin:
            fin.seek(sz)
            d = np.fromfile(fin, dtype=np.uint8)
            print(('{} contains {} zeros and {} nonzeros'.format(f, (d==0).sum(), (d!=0).sum())))
            
    

if __name__ == '__main__':
    _main()
