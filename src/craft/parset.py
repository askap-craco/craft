#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Parset:
    def __init__(self, s):
        self.d = {}
        for line in s.split('\n'):
            if '=' not in line:
                continue
            line = line.split('#')[0]
            
            name, value = line.strip().split('=')
            name = name.strip()
            value = value.strip()
            if value.startswith('[') and value.endswith(']'):
                value = value[1:-1].split(',')
            self.d[name] = value

    def keys(self):
        return self.d.keys()

    def __getitem__(self, name):
        return self.d[name]
            
    @classmethod
    def from_file(self, fname):
        with open(fname, 'rt') as f:
            s = f.read()
            return Parset(s)
    


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


