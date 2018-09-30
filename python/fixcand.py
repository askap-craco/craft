#!/usr/bin/env python
"""
Fix candidate files - make MJD referenced to infinity

Copyright (C) CSIRO 2018
"""
import matplotlib as mpl
import numpy as np
import fredfof
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-f', '--ref-freq', type=float, help='Low end of the band that s the reference frequency of the MJD (MHz)', default=(1488.-336.))
    parser.add_argument('--suffix', default='.finf', help='Suffix for the output files')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        fix(f, values)

def fix(f, values):
    ref_freq = values.ref_freq/1e3
    new_file = f +'.finf'
    fout = open(new_file, 'w')
    for line in open(f, 'rU'):
        if line.startswith('#'):
            fout.write(line)
            continue
            
        line = line.strip()
        bits = line.split()
        dm = float(bits[5])
        mjd = float(bits[7])
        offset_ms = 4.15*dm*(0 - ref_freq**-2)
        offset_days = offset_ms/ 86400. / 1e3
        bits[7] = '%0.15f'%(mjd + offset_days)
        fout.write(' '.join(bits) + '\n')

    fout.close()

if __name__ == '__main__':
    _main()
