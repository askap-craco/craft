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
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', help="beam offsets, one per file, comma separated")
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    nfiles = len(values.files)
    offsets = list(map(int, values.offset.split(',')))
    assert nfiles == len(offsets)
    for fin_name, boff in zip(values.files, offsets):
        fin = open(fin_name, 'rU')
        print('# new file %s beam offset %d' % (fin_name, boff))
        for line in fin:
            if line.strip().startswith('#'):
                print(line)
            else:
                bits = line.strip().split()
                if len(bits) == 0:
                    continue
                else:
                    new_beamno = int(bits[-1]) + boff
                    bits[-1] = str(new_beamno)
                    print(' '.join(bits))

        fin.close()

if __name__ == '__main__':
    _main()
