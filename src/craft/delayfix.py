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

def loadfile(f):
    d = {}
    for line in open(f, 'rU'):
        bits = line.split()
        d[bits[0]] = int(bits[1])

    return d

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

    

    f1 = loadfile(values.files[0])
    f2 = loadfile(values.files[1])
    allkeys = list(f1.keys()) + list(f2.keys())
    for k in sorted(allkeys):
        print(k, f1.get(k, 0) + f2.get(k,0))
        
        
    

if __name__ == '__main__':
    _main()
