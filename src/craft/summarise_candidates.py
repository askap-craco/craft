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
from .plot_fredda_cand import  find_files
from .crafthdr import DadaHeader
import glob

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def get_dada_header(f):
    thedir = os.path.dirname(f)
    files = glob.glob(os.path.join(thedir, 'ak*.hdr.fixed'))
    assert len(files) == 1, 'Too many headers in {}'.format(thedir)
    hdr = DadaHeader.fromfile(files[0])
    return hdr

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument('-w','--max-boxcar', help='max width to plot', default=32, type=int)
    parser.add_argument('-d','--min-dm', help='minimum dm to plot', default=0, type=float)
    parser.add_argument('-b','--min-sn', help='minimum S/N to plot', default=0, type=float)

    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    npass_all = 0
    for d in values.files:
        npass_d = 0
        for f in find_files(d, values.fname):
            hdr = get_dada_header(f)
            source = hdr['SOURCE'][0]
            try:
                ncand, npass = count_ncand(f, values)

                print(d, f, source, ncand, npass)
                if 'PSR' not in source:
                    npass_d += npass
                    npass_all += npass
            except:
                logging.exception('Could not count ncand %s', f)
            
        print(d, npass_d)

    print('TOTAL', npass_all)

def count_ncand(f, values):
    vin = np.loadtxt(f)
    if len(vin.shape) == 1:
        vin.shape = (1, len(vin))

    boxcar = vin[:, 3]
    dm = vin[:, 5]
    sn = vin[:, 0]
    mask = (boxcar < values.max_boxcar) & (dm >= values.min_dm) & (sn >= values.min_sn)
    ncand = vin.shape[0]
    npass = sum(mask)

    return ncand, npass 





    

if __name__ == '__main__':
    _main()
