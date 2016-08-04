#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def load_candidates(d):
    '''
    Loads all candidate files at or below directory 
    :d: Root directory for candidate files
    '''
    all_data = []
    for dirpath, dirnames, filenames in os.walk(d):
        for f in filenames:
            if f.endswith('.cand'):
                fname = os.path.join(dirpath, f)
                if os.path.getsize(fname) != 0:
                    #all_files.append(fname)
                    cand = np.loadtxt(fname)
                    if cand.ndim == 1:
                        cand.shape = (1, len(cand))

                    all_data.extend(cand)

    all_data = np.array(all_data)

    return all_data


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
