#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from astropy.io import fits
import shutil
import pylab

__author__ = "Keith Bannister <keith.bannister@csiro.au>"



def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--output', help='output file')
    parser.add_argument(dest='files', nargs=2)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    shutil.copy(values.files[0], values.output)        
    f2 = fits.open(values.files[1])
    fout = fits.open(values.output, mode='update')
    dout = fout['PRIMARY'].data
    d2 = f2['PRIMARY'].data
    assert dout.shape == d2.shape
    assert np.all(dout['BASELINE'] == d2['BASELINE'])
    fout['PRIMARY'].data['DATA'] -= f2['PRIMARY'].data['DATA']
    fout.flush()


if __name__ == '__main__':
    _main()
