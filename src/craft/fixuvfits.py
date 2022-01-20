#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
import pyfits
import numpy as np
from astropy.time import Time

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

    fix_bzero = False

    for f in values.files:
        hdu = pyfits.open(f)
        h = hdu[0]
        twhole = h.data['DATE'].astype(float)
        if '_DATE' in h.columns.names:
            tpart = h.data['_DATE'].astype(float)
        else:
            tpart = np.zeros(1)
        tall = tpart + twhole
        t0 = Time(np.round(tall.min()) - 0.5, format='jd')
        datec = h.columns['DATE']
        twholec = datec.array
        tnew = tall - t0.jd
        print('Column = {} twhole={} twholec={} tpart={} tall={} t0={} tnew={}'.format(datec, twhole[0], twholec[0], tpart[0], tall[0], t0, tnew[0]))
        if datec.bzero == 0 and fix_bzero:
            print('Fixing DATE bzero')
            t0.format = 'fits'
            #h.header.set('DATE-OBS', t0.value, 'Fixed by fixuvfits.py')
            #h.header.set('DATE_OBS', t0.value, 'Fixed by fixuvfits.py')
            h.data['DATE'] = tnew.astype(np.float32)
            #datec.array = tnew.astype(np.float32)*0

            datec.bzero = t0.jd
            if '_DATE' in h.columns.names:
                h.data['_DATE'] = 0
                

            print('New data', h.data['DATE'].min(), h.data['_DATE'].min(), 'writing to', fout, h.header['DATE_OBS'])
            

        del hdu[0].header['DATE_OBS']
        del hdu[2].header['RDATE']
        fout = f+'.fixed'
        hdu.writeto(fout)

if __name__ == '__main__':
    _main()
