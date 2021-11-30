#!/usr/bin/env python
"""
UV fits reader class and utility

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from astropy.io import fits
import craco

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class UvFits(object):

    def __init__(self, hdulist):
        self.hdulist = hdulist

    @property
    def vis(self):
        '''Returns visbility table
        '''
        return self.hdulist[0].data

    @property
    def start_date(self):
        d0 = self.vis[0]['DATE']
        return d0

    @property
    def channel_frequencies(self):
        return craco.get_freqs(self.hdulist)

    @property
    def baselines(self):
        '''
        Returns all data from first integration
        
        :returns: dictionary, keyed by baseline ID of all basesline data with a timestamp
        equal to the first timestamp in the file
        '''
            
        d0 = self.start_date
        baselines = {}
        vis = self.vis
        for i in range(self.vis.size):
            row = vis[i]
            baselines[row['BASELINE']] = row
            if row['DATE'] != d0:
                break

        return baselines

    def get_max_uv(self):
        ''' 
        Return the largest absolute values of UU and VV in lambdas
        '''
        fmax = self.channel_frequencies.max()
        baselines = self.baselines
        ulam_max = max([abs(bldata['UU'])*fmax for bldata in baselines.values()])
        vlam_max = max([abs(bldata['VV'])*fmax for bldata in baselines.values()])
        return (ulam_max, vlam_max)

    
    def plot_baselines(self):
        baselines = self.baselines
        freqs = self.channel_frequencies
        for blid, bldata in baselines.iteritems():
            ulam = bldata['UU'] * freqs
            vlam = bldata['VV'] * freqs
            
            pylab.plot(ulam/1e3, vlam/1e3)
        
        pylab.xlabel('U (klambda)')
        pylab.ylabel('V (klambda)')
        #pylab.show()

    def time_blocks(self, nt):
        '''
        Returns a sequence of baseline data in blocks of nt
        '''
        return craco.time_blocks(self.vis, nt)

def open(*args, **kwargs):
    logging.info('Opening file %s', args[0])
    return UvFits(fits.open(*args, **kwargs))

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
