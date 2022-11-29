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
from . import craco
from .craco import bl2ant
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class UvFits(object):

    def __init__(self, hdulist, max_nbl=None, mask=True):
        '''
        @param hdulist FITS HDLISt typically got from pyfits.open
        @param max_nbl - if not none, only return this many baselines
        '''
        self.hdulist = hdulist
        self.max_nbl = max_nbl
        self.flagant = []
        self.ignore_autos = True
        self.mask = mask

    def set_flagants(self, flagant):
        '''
        set list of 1-based antennas to flag
        '''
        self.flagant = flagant
        return self

    @property
    def header(self):
        return self.hdulist[0].header

    @property
    def vis(self):
        '''Returns visbility table
        '''
        return self.hdulist[0].data

    @property
    def start_date(self):
        row = self.vis[0]
        d0 = row['DATE']
        try:
            d0 += row['_DATE'] # FITS standard says add these two columns together
        except KeyError:
            pass

        return d0

    @property
    def channel_frequencies(self):
        return craco.get_freqs(self.hdulist)

    @property
    def baselines(self):
        '''
        Returns all data from first integration
        Doesn't include baselines containing flagant
        
        :returns: dictionary, keyed by baseline ID of all basesline data with a timestamp
        equal to the first timestamp in the file
        '''
            
        d0 = self.start_date
        baselines = {}
        vis = self.vis
        for i in range(self.vis.size):
            row = vis[i]
            blid = row['BASELINE']
            a1, a2 = bl2ant(blid)
            if a1 in self.flagant or a2 in self.flagant:
                continue

            if self.ignore_autos and a1 == a2:
                continue
            
            baselines[blid] = row
            if row['DATE'] != d0 or (self.max_nbl is not None and i > self.max_nbl):
                break

        return baselines

    def get_max_uv(self):
        ''' 
        Return the largest absolute values of UU and VV in lambdas
        '''
        fmax = self.channel_frequencies.max()
        baselines = self.baselines
        ulam_max = max([abs(bldata['UU'])*fmax for bldata in list(baselines.values())])
        vlam_max = max([abs(bldata['VV'])*fmax for bldata in list(baselines.values())])
        return (ulam_max, vlam_max)

    
    def plot_baselines(self):
        baselines = self.baselines
        freqs = self.channel_frequencies
        for blid, bldata in baselines.items():
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
        # WARNING TODO: ONLY RETURN BASELINES THAT HAVE BEEN RETURNED in .baselines
        # IF max_nbl has been set
        return craco.time_blocks(self.vis, nt, self.flagant, self.ignore_autos, mask=self.mask)
    
    def get_tstart(self):
        '''
        return tstart as astropy Time from header otherwise first row of fits table
        '''
        f = self
        if 'DATE-OBS' in f.header:
            tstart = Time(f.header['DATE-OBS'], format='isot', scale='utc')
        else:
            jdfloat = f.start_date
            tstart = Time(jdfloat, format='jd', scale='utc')
        
        return tstart
        
    def get_target_position(self, targidx=0):
        '''
        return (ra,dec) degrees from header if available, otherwise source table
        '''
        f = self
        if 'OBSRA' in f.header:
            ra = f.header['OBSRA'] * u.degree
            dec = f.header['OBSDEC'] * u.degree
            log.info('Got radec=(%s/%s) from OBSRA header', ra, dec)
        else:
            source_table = f.hdulist[3].data
            assert len(source_table)==1, f'Dont yet support multiple source files: {len(source_table)}'
            row = source_table[targidx]
            src = row['SOURCE']
            ra = row['RAEPO']*u.degree
            dec = row['DECEPO']*u.degree
            log.info('Got radec=(%s/%s) from source table for %s', ra, dec, src)

        return (ra, dec)
            


    def close(self):
        return self.hdulist.close()

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
