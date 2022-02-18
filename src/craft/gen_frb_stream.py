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

def frb_data_stream(yaml_file):
    '''
    Generates a stream of FRB data specified by the given YAML file

    @param yaml_file path to yaml file
    @returns generator of numpy arrays. Each array is [nbl, nc, nt], dtype=np.complex64, orrrrr, [nbl, nc, nt, 2] dtype=np.int16
    '''

    # this is not working code - just an idea of how you'd do it.

    frb_config = parse_file(yaml_file)
    nt = frb_config.nt
    sampnum = 0
    current_frb = None
    for iblck in range(nblocks):
        if current_frb is None:
            current_frb, start_sample = make_new_frb(frb_config, iblk)
            # FRB number of times might be >> the block size

        data = make_blank_data_with_noise(frb_config)
        is_finished = add_frb_to_data_at_position(data, curent_frb, samp_num)
        if is_finished:
            current_frb = None
        
        sampnum += nt
        
        yield data
    

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
