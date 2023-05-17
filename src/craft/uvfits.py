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
import warnings
from astropy.io import fits
from . import craco
from .craco import bl2ant,get_max_uv
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord


log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class UvFits(object):

    def __init__(self, hdulist, max_nbl=None, mask=True, skip_blocks=0):
        '''
        @param hdulist FITS HDLISt typically got from pyfits.open
        @param max_nbl - if not none, only return this many baselines
        @param skip_blocks - number of time blocks to skip
        '''
        self.hdulist = hdulist
        self.max_nbl = max_nbl
        self.flagant = []
        self.ignore_autos = True
        self.mask = mask

        # first we calculate the number of baselines in a block
        assert skip_blocks >= 0, f'Invalid skip_blocks={skip_blocks}'
        self.__nstart = 0
        startbl = self.baselines
        self.nbl = len(startbl)
        # next we set the start block for the rest of time
        self.__nstart = self.nbl*skip_blocks
        self.skip_blocks = skip_blocks
        nrows = len(self.hdulist[0].data)
        log.debug('File contains %d baselines. Skipping %d blocks with nstart=%d', self.nbl, skip_blocks, self.__nstart)

        if self.__nstart >= nrows:
            raise ValueError(f'Requested skip {skip_blocks} larger than file. nrows={nrows} nbl={self.nbl} nstart={self.__nstart}')

        

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
        '''
        Returns visbility table
        if skip_blocks is > in teh constructor, then that number of blocks will have been skipped and you won't see them

        '''
        return self.hdulist[0].data[self.__nstart:]

    @property
    def start_date(self):
        '''
        Returns MJD float of first (skipped) sample in the file
        '''
        row = self.vis[0]
        d0 = row['DATE']
        try:
            d0 += row['_DATE'] # FITS standard says add these two columns together
        except KeyError:
            pass

        return d0

    @property
    def tsamp(self):
        '''
        Return sample time in seconds
        '''
        try:
            ts = self.vis[0]['INTTIM']*u.second # seconds
        except KeyError:
            warnings.warn('Unknown int time in file. returning 1ms')
            ts = 1e-3*u.second # seconds

        return ts

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
        return get_max_uv(baselines, fmax)

    
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

    def time_blocks_with_uvws(self, nt):
        '''
        Returns a sequence of baseline data in blocks of nt
        '''
        # WARNING TODO: ONLY RETURN BASELINES THAT HAVE BEEN RETURNED in .baselines
        # IF max_nbl has been set
        return craco.time_blocks_with_uvws(self.vis, nt, self.flagant, self.ignore_autos, mask=self.mask, fetch_uvws = True)
    
    def time_blocks(self, nt):
        '''
        Returns a sequence of baseline data in blocks of nt
        '''
        # WARNING TODO: ONLY RETURN BASELINES THAT HAVE BEEN RETURNED in .baselines
        # IF max_nbl has been set
        return craco.time_blocks(self.vis, nt, self.flagant, self.ignore_autos, mask=self.mask)
    
    def time_block_with_uvw_range(self, trange):
        """
        return a block of data and uvw within a given index range
        """
        return craco.time_block_with_uvw_range(
            vis=self.vis, trange=trange, flagant=self.flagant,
            flag_autos=self.ignore_autos, mask=self.mask
        )

    def get_tstart(self):
        '''
        return tstart as astropy Time from header otherwise first row of fits table
        '''
        f = self
        first_sample_date = Time(f.start_date, format='jd', scale='utc')

        if 'DATE-OBS' in f.header:
            tstart = Time(f.header['DATE-OBS'], format='isot', scale='utc')
            assert self.skip_blocks == 0, f'This date will be be incorrect if skip_blocks={self.skip_blocks} is nonzero - header={tstart} but first integration is {first_sample_date}. This code needs ot be fixed if you neeed it'
        else:
            tstart = first_sample_date
        
        return tstart

    @property
    def tstart(self):
        return self.get_tstart()

    @property
    def source_table_entry(self):
        f = self
        source_table = f.hdulist[3].data
        first_datarow = next(iter(self.baselines.values()))
        first_targetidx = int(first_datarow['SOURCE'])
        row = source_table[first_targetidx-1] # FITS convention is 1 based
        
        if len(source_table) > 1:
            warnings.warn(f'Dont yet support multiple source files: {len(source_table)} - using source at {first_targetidx} which is {row}')
            
        return row

    
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
            row = self.source_table_entry
            src = row['SOURCE']
            ra = row['RAEPO']*u.degree
            dec = row['DECEPO']*u.degree
            log.info('Got radec=(%s/%s) from source table for %s', ra, dec, src)

        return (ra, dec)

    def get_target_skycoord(self):
        '''
        Returns skycoord of phase center
        Leave this for posterity
        '''
        (ra, dec) = self.get_target_position()
        coord = SkyCoord(ra, dec, frame='icrs')
        return coord

    @property
    def target_skycoord(self):
        '''
        Return skycorod of phase center. Nicer version.
        '''
        return self.get_target_skycoord()

    @property
    def target_name(self, targidx=0):
        f = self
        if 'OBJECT' in f.header:
            src =f.header['OBJECT']
        else:
            row = self.source_table_entry
            src = row['SOURCE']
            
        return src

    def close(self):
        return self.hdulist.close()

def open(*args, **kwargs):
    logging.info('Opening file %s', args[0])
    return UvFits(fits.open(*args, **kwargs), **kwargs)

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--skip-blocks', type=int, default=0, help='Skip this many bllocks in teh UV file before usign it for UVWs and data')

    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = open(values.files[0], skip_blocks=values.skip_blocks)
    f.plot_baselines()
    pylab.show()
    
    

if __name__ == '__main__':
    _main()
