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
from .freq_config import FrequencyConfig
from .craco import bl2ant,get_max_uv
from craft.vis_metadata import VisMetadata
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import builtins

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def parse_beam_id_from_filename(fname):
    '''
    return beamID as an int from a filebame that looks like bXX.uvfits
    '''
    fname = os.path.basename(fname)
    if fname.startswith('b') and fname.endswith('.uvfits'):
        b  = int(fname.split('.uvfits')[0][1:])
    else:
        raise ValueError(f'Filename {fname} doesnt contain beam')

    return b

class VisView:
    def __init__(self, uvfitsfile, start_idx):
        '''
        Allows you to read and seek rather than memory map
        '''
        self.uvfitsfile = uvfitsfile
        self.dtype = uvfitsfile.hdulist[0].data.dtype
        self.hdrsize = len(str(self.uvfitsfile.hdulist[0].header))
        assert self.hdrsize % 2880 == 0
        self.start_idx = start_idx
        self.fin = self.uvfitsfile.raw_fin

    @property
    def size(self):
        sz = self.uvfitsfile.hdulist[0].header['GCOUNT'] - self.start_idx
        assert sz >= 0
        return sz
        
    def __getitem__(self, sidx: slice):
        '''
        Slices data from UV fitsfile
        '''
        
        sidx_ = sidx # keep a copy of the input slice...
        if isinstance(sidx, int):
            sidx = slice(sidx, sidx+1, 1)
            
        assert sidx.step is None or sidx.step == 1, 'np.fromfile cant to strided access'
        assert sidx.stop >= sidx.start, 'Cant read backwards'
        
        start_byte  = self.hdrsize + (sidx.start + self.start_idx)*self.dtype.itemsize
        nelements = sidx.stop - sidx.start
        assert nelements >= 0
        
        self.fin.seek(start_byte)
        dout = np.fromfile(self.fin, count=nelements, dtype=self.dtype)

        ### update date to float64
        dtype_store = [
            (field, np.dtype(">f8")) if field == "DATE" 
            else (field, self.dtype.fields[field][0])
            for field in self.dtype.fields
        ]

        dout = dout.astype(dtype_store)
        dout["DATE"] += self.uvfitsfile.hdulist[0].header["PZERO4"]

        if isinstance(sidx_, int): return dout[0]
        return dout

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
        self.raw_fin = builtins.open(self.hdulist.filename(), 'rb')

        # first we calculate the number of baselines in a block
        assert skip_blocks >= 0, f'Invalid skip_blocks={skip_blocks}'
        self.__nstart = 0
        startbl = self._find_baseline_order()
        self.nbl = len(startbl)
        # next we set the start block for the rest of time
        self.__nstart = self.nbl*skip_blocks
        self.skip_blocks = skip_blocks
        nrows = len(self.hdulist[0].data)
        log.debug('File contains %d baselines. Skipping %d blocks with nstart=%d', self.nbl, skip_blocks, self.__nstart)
        self.nblocks = nrows // self.nbl
        self._freq_config = FrequencyConfig.from_hdu(self.hdulist[0])

        if skip_blocks >= self.nblocks:
            raise ValueError(f'Requested skip {skip_blocks} larger than file. nblocks = {self.nblocks} nrows={nrows} nbl={self.nbl} nstart={self.__nstart}')


    def set_flagants(self, flagant):
        '''
        set list of 1-based antennas to flag
        '''
        self.flagant = flagant
        return self
        

    @property
    def valid_ants(self):
        '''
        Return list of valid antennas
        It's 36 antennas less the antennas that are listed in flagant
        1-based
        TODO: Load up antenna table
        '''
        all_ants = set([a+1 for a in range(36)])
        va = sorted(list(all_ants - set(self.flagant)))
        return va

    @property
    def valid_ants_0based(self):
        '''
        Returns np array of valid ants (0 based)
        '''
        return np.array(self.valid_ants) - 1
        

    @property
    def filename(self):
        return self.hdulist.filename()

    @property
    def header(self):
        return self.hdulist[0].header

    @property
    def vis(self):
        '''
        Returns visbility table
        if skip_blocks is > in teh constructor, then that number of blocks will have been skipped and you won't see them

        '''
        return VisView(self, self.__nstart)

    @property
    def start_date(self):
        '''
        Returns MJD float of first (skipped) sample in the file
        '''
        row = self.vis[0]
        d0 = row['DATE']
        # seems like if you read it from the raw data... don't do that!!!
        # try:
        #     d0 += row['DATE'] # FITS standard says add these two columns together
        # except KeyError:
        #     pass

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
        '''
        Returns channel frequencies in Hz
        '''
        return self.freq_config.channel_frequencies*1e6

    @property
    def freq_config(self):
        return self._freq_config

    
    def _find_baseline_order(self):
        '''
        Returns all data from first integration
        Doesn't include baselines containing flagant
        Set up internal baseline ordering thigns for fast_time_blocks
        
        :returns: dictionary, keyed by baseline ID of all basesline data with a timestamp
        equal to the first timestamp in the file
        '''
            
        d0 = self.start_date
        baselines = {}
        internal_baseline_order = []
        baseline_indices = []
        raw_nbl = 0
        vis = self.vis
        for i in range(self.vis.size):
            row = vis[i]
            if row['DATE'] != d0 or (self.max_nbl is not None and i > self.max_nbl):
                break
            blid = row['BASELINE']
            a1, a2 = bl2ant(blid)
            raw_nbl += 1
            if a1 in self.flagant or a2 in self.flagant:
                continue

            if self.ignore_autos and a1 == a2:
                continue
            
            baselines[blid] = row
            internal_baseline_order.append(blid)
            baseline_indices.append(i)


        self.internal_baseline_order = np.array(internal_baseline_order)
        self.baseline_indices = np.array(baseline_indices)
        self.raw_nbl = raw_nbl
        return baselines

    @property
    def baselines(self):
        return self.get_uvw_at_isamp(0)

    def get_uvw_at_isamp(self, isamp:int):
        d, uvw = next(self.fast_time_blocks(1, fetch_uvws=True, istart=isamp))
        return uvw[0]

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


    @property
    def nsamps(self):
        '''
        Returns number of samples in the file
        '''
        n =  self.vis.size // self.raw_nbl
        return n
    
    def fast_time_blocks(self, nt, fetch_uvws=False, istart=0):
        '''
        Reads raw data from uvfits file as an array and returns a block of nt samples, along with uvws as a list (which is empty if fetch_uvws is False)
        :istart: Sample number to start at. Doesnt have to be a multiple of nt

        Returns a tuple - (dout, uvws)
        '''
        assert istart >= 0, f'Invalid istart={istart}'
        samps_returned = istart
        samps_to_read = nt
        uvws = []
        if istart > self.nsamps - 1:
            raise ValueError(f'Asked to start reading past the end of the file. istart={istart} nsamps={self.nsamps}')

        while True:
            samps_left = self.nsamps - samps_returned

            if samps_left < nt:
                samps_to_read = samps_left

            if samps_left < 1:
                break

            byte_offset = self.vis.hdrsize + samps_returned * self.raw_nbl * self.vis.dtype.itemsize
            nbytes = samps_to_read*self.raw_nbl*self.vis.dtype.itemsize
            next_byte_offset = byte_offset+nbytes
            
            self.vis.fin.seek(byte_offset)

            log.debug('Reading %d bytes from %s at offset %d', nbytes, self.filename, byte_offset)
            dout = np.fromfile(self.vis.fin, count = samps_to_read * self.raw_nbl, dtype=self.vis.dtype).reshape(samps_to_read, -1)
            log.debug('read complete')

            # Tell Kernel we're going to need the next block. - doesnt seem to make much difference, but anyway.
            os.posix_fadvise(self.vis.fin.fileno(), next_byte_offset, nbytes, os.POSIX_FADV_WILLNEED)
            
            samps_returned += samps_to_read
            
            dout_data = dout['DATA'].transpose((*np.arange(1, dout['DATA'].ndim), 0))[self.baseline_indices]
            dout_complex_data = dout_data[..., 0, :] + 1j*dout_data[..., 1, :]

            if self.mask:
                mask_vals = (1 - dout_data[..., 2, :]).astype('int')
                dout_complex_data = np.ma.MaskedArray(data = dout_complex_data, mask = mask_vals)

            if fetch_uvws:
                uvws = []
                for it in range(samps_to_read):
                    this_uvw = {}

                    for ii, ibl in enumerate(self.baseline_indices):
                        blid = self.internal_baseline_order[ii]

                        # add the [0] index to get the scalar type
                        # so it exactly matches what you with the original
                        # .baselines numpy data
                        # There might be a better way to do it but I don't know how
                        this_uvw[blid] = np.array([(
                                            dout['UU'][it, ibl], 
                                            dout['VV'][it, ibl], 
                                            dout['WW'][it, ibl]
                                            )],
                                            dtype=[('UU', dout['UU'].dtype),
                                                   ('VV', dout['VV'].dtype),
                                                   ('WW', dout['WW'].dtype)
                                                ]
                                            )[0]

                    uvws.append(this_uvw)

            yield dout_complex_data, uvws


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

        # early version of MPIPIPELINE before 10 Nov wrote TAI
        # rather than UTC scale to the timestamps
        # now we write a 'TSCALE' header
        scale = f.header.get('TSCALE', 'TAI').lower()
        first_sample_date = Time(f.start_date, format='jd', scale=scale)

        if 'DATE-OBS' in f.header:
            tstart = Time(f.header['DATE-OBS'], format='isot', scale=scale)
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
        first_datarow = self.vis[0]
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

    @property
    def beamid(self):
        '''
        Returns beamId from header or from filename if filename is bXX.uvfits
        '''
        
        b = self.header.get('BEAMID', -1)
        if b == -1:
            try:
                b = parse_beam_id_from_filename(self.hdulist.filename())
            except ValueError:
                b = 0
                warnings.warn('No BeamID in header. Filename not useful. Returning 0')
        return b

    def vis_metadata(self, isamp:int):
        '''
        Return a vis info adapter for the given sample number
        '''
        tstart = self.tstart + self.tsamp*isamp
        uvw = self.get_uvw_at_isamp(isamp)
        m = VisMetadata(
            uvw,
            self.freq_config,
            self.target_name,
            self.target_skycoord,
            self.beamid,
            tstart,
            self.tsamp,
            tstart)
        m.isamp = isamp

        return m


    def close(self):
        self.raw_fin.close()
        self.raw_fin = None
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
