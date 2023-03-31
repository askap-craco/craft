#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from .calc11 import CalcFile

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

    antd, antnos = load_parset(values.files[0], prefix='common.antenna.ant')
    foutname = values.files[0]+'.calc'
    calcfile = CalcFile()
    calcfile.add_antdata(antd, antnos)
    calcfile.writeto(foutname)

def load_parset(fcmfile, prefix='common.antenna.ant'):
    f = open(fcmfile, 'rU')
    ant_data = {}
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue

        bits = line.split('=')
        if len(bits) != 2:
            continue
        
        key, value = line.split('=')
        key = key.strip()
        value = value.strip()

        if line.startswith(prefix):
            bits = key.split('.')
            antname = bits[2]
            #print(line, antname)
            if len(antname) <= 3:
                continue

            antno = int(antname[3:])
            smallkey = '.'.join(bits[3:])
            d = ant_data.get(antno, {})
            d[smallkey] = value
            ant_data[antno] = d
        if key == 'common.antennas':
            antnos = [int(x.replace('ant','')) for x in value.replace('[','').replace(']','').split(',')]



    return ant_data, antnos






if __name__ == '__main__':
    _main()
