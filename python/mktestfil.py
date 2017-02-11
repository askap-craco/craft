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
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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

    hdr = {'source_name':'ones',
           'tsamp': 0.001265625,
           'pulsarcentric':0,
           'az_start':1.,
           'nbits':8,
           'foff':-1.0,
           'fch1':1448.0,
           'nchans': 336,
           'telescope_id':1,
           'src_dej': -30724.3618819,
           'src_raj': 164015.762746,
           'tstart':57681.5338261,
           'nifs':1
       }
    
    sfile = sigproc.SigprocFile(values.files[0], 'w', hdr)
    sfile.seek_data(0)
    nchans = hdr['nchans']
    nifs = hdr['nifs']
    ntimes = 16384
    d = np.ones((ntimes, nifs, nchans), dtype=np.uint8)*128
    print 'Writing ', d.shape, d.dtype, 'to', values.files[0]
    d.tofile(sfile.fin)
    sfile.fin.flush()
    sfile.fin.close()
    

if __name__ == '__main__':
    _main()
