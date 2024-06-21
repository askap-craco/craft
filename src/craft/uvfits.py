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
from collections import OrderedDict

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
    def __init__(self, uvfitsfile, update_date=True):
        '''
        Allows you to read and seek rather than memory map
        '''
        start_idx = uvfitsfile._nstart
        self.uvfitsfile = uvfitsfile
        self.hdrsize = uvfitsfile.hdrsize
        self.dtype = uvfitsfile.dtype
        assert self.hdrsize % 2880 == 0
        self.start_idx = start_idx
        self.fin = self.uvfitsfile.raw_fin
        self.update_date = update_date

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
            
        assert sidx.step is None or sidx.step == 1, 'np.fromfile cant do strided access'
        assert sidx.stop >= sidx.start, 'Cant read backwards'
        if sidx.start is None:
            sidx.start = 0
        
        start_byte  = self.hdrsize + (sidx.start + self.start_idx)*self.dtype.itemsize
        nelements = sidx.stop - sidx.start
        assert nelements >= 0
        
        self.fin.seek(start_byte)
        dout = np.fromfile(self.fin, count=nelements, dtype=self.dtype)
        if len(dout) != nelements:
            raise EOFError(f'Read past the end of the file slice={sidx} start_byte={start_byte} count={nelements} actual={len(dout)} dtype={self.dtype}')

        ### update date to float64

        if self.update_date:
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

    def __init__(self, hdulist, max_nbl=None, mask=True, skip_blocks=0,
                 start_mjd=None, end_mjd=None):
        '''
        @param hdulist FITS HDLISt typically got from pyfits.open
        @param max_nbl - if not none, only return this many baselines
        @param skip_blocks - number of time blocks to skip
        @param start_mjd - if not None skip to start MJD - ignore skip_blocks
        @param end_mjd - don't iterate past this mjd
        '''
        self.hdulist = hdulist
        self.max_nbl = max_nbl
        self.flagant = []
        self.ignore_autos = True
        self.mask = mask
        self.raw_fin = builtins.open(self.hdulist.filename(), 'rb')
        self.hdrsize = len(str(self.hdulist[0].header))
        self.dtype = self.hdulist[0].data.dtype
        self._nstart = 0

        if start_mjd is not None:
            assert skip_blocks == 0, 'cant do both. Very tricky'
            skip_blocks = int(np.ceil(self.time_to_sample(start_mjd)))
            skip_blocks = max(skip_blocks, 0)
            
        assert skip_blocks >= 0, f'Invalid skip_blocks={skip_blocks}'

        if end_mjd is not None:
            iend = int(np.floor(self.time_to_sample(end_mjd)))
            assert 0 <= iend,'Invalid end sample for {iend} for time {end_mjd}'
            iend = min(iend, self.nblocks)
        else:
            iend = None

        self.iend = iend
        self.skip_blocks = skip_blocks

        self._freq_config = FrequencyConfig.from_hdu(self.hdulist[0])
        self._reset_baseline_setup()

    def _reset_baseline_setup(self):
        startbl = self._find_baseline_order()
        skip_blocks = self.skip_blocks
        self.nbl = len(startbl)
        # next we set the start block for the rest of time
        self._nstart = self.raw_nbl*skip_blocks
        self.skip_blocks = skip_blocks
        nrows = len(self.hdulist[0].data)
        self.nblocks_raw = nrows // self.raw_nbl
        if skip_blocks >= self.nblocks_raw:
            raise ValueError(f'Requested skip {skip_blocks} larger than file. nblocks = {self.nblocks_raw} nrows={nrows} nbl={self.nbl} nstart={self._nstart}')
        self.nblocks = self.nblocks_raw - skip_blocks


    def set_flagants(self, flagant):
        '''
        set list of 1-based antennas to flag
        '''
        self.flagant = flagant
        self._reset_baseline_setup()
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

    def gcount(self):
        return self.header['GCOUNT']

    @property
    def vis(self):
        '''
        Returns visbility table
        if skip_blocks is > in teh constructor, then that number of blocks will have been skipped and you won't see them

        '''
        return VisView(self)

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
    def end_date(self):
        '''
        Returns MJD float of final (skipped) sample in the file
        '''
        # Note: the final row int the table is occasionally corrupted for some reason
        # so we'll just calculate it from tstart + nsamp
        d0 = self.sample_to_time(self.nsamps)
        return d0

    @property
    def tsamp(self):
        '''
        Return sample time in seconds
        '''
        try:
            ts = self.vis[0]['INTTIM']*u.second # seconds
        except ValueError:
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
        raw_baseline_order = []
        baseline_indices = []
        raw_baseline_indices = []
        raw_nbl = 0
        vis = self.vis
        for i in range(self.vis.size):
            row = vis[i]
            if row['DATE'] != d0 or (self.max_nbl is not None and i > self.max_nbl):
                break
            blid = row['BASELINE']
            a1, a2 = bl2ant(blid)
            raw_nbl += 1
            raw_baseline_order.append(blid)
            raw_baseline_indices.append(i)
            if a1 in self.flagant or a2 in self.flagant:
                continue

            if self.ignore_autos and a1 == a2:
                continue
            
            baselines[blid] = row
            internal_baseline_order.append(blid)
            baseline_indices.append(i)


        self.internal_baseline_order = np.array(internal_baseline_order)
        self.raw_baseline_order = np.array(raw_baseline_order)
        self.baseline_indices = np.array(baseline_indices)
        self.raw_baseline_indices = np.array(raw_baseline_indices)
        self.raw_nbl = raw_nbl
        return baselines

    @property
    def baselines(self):
        return self.get_uvw_at_isamp(0)

    def sample_to_time(self, isamp:int):
        '''
        Returns sample number to astropy time
        '''
        assert isamp >= 0, 'Invalid isamp'
        t = self.tstart + isamp*self.tsamp
        return t


    def time_to_sample(self, time:Time)->float:
        '''
        Returns time in samples (can be float)
        between supplied time and the start time in units of TSAMP
        '''
        if not isinstance(time, Time):
            time = Time(time, format='mjd', scale='tai')
            #Andy asked me to use scale = tai here. If it is wrong -- kill Andy

        tdiff = time - self.tstart
        tdiff_samp = float(tdiff / self.tsamp)
        return tdiff_samp

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

    def _create_masked_data(self, dout_data, start_sampno):
        mask_vals = dout_data[..., 2, :] <= 0
        #mask_vals = (1 - dout_data[..., 2, :]).astype('int')
        dout_complex_data = dout_data[..., 0, :] + 1j*dout_data[..., 1, :]

        if self.mask:
            dout_complex_data = np.ma.MaskedArray(data = dout_complex_data, mask = mask_vals)
            
        return dout_complex_data

    def fast_raw_blocks(self, istart=0, nsamp=None, nt=1, raw_date=False):
        assert istart >= 0, f'Invalid istart={istart}'
        if nsamp is None:
            nsamp = self.nsamps

        assert nsamp > 0
        endsamp = min(istart + nsamp, self.nsamps - nt)

        if istart > self.nsamps - nt:
            raise ValueError(f'Asked to start reading past the end of the file. istart={istart} nsamps={self.nsamps}')

        vis = self.vis
        

        for startsamp in range(istart, endsamp, nt):           
            ix = startsamp*self.raw_nbl
            iy = ix + nt*self.raw_nbl
            byte_offset = self.hdrsize + ix * self.dtype.itemsize
            nbytes = nt*self.raw_nbl*self.dtype.itemsize
            next_byte_offset = byte_offset+nbytes

            log.debug('Reading %d bytes from %s at offset %d', nbytes, self.filename, byte_offset)            

            dout = vis[ix:iy].reshape(nt, -1)

            log.debug('read complete of visidx [%d:%d]. Pre-reading willneed %d for %d bytes', ix, iy, next_byte_offset, nbytes)

            # Tell Kernel we're going to need the next block. - doesnt seem to make much difference, but anyway.
            os.posix_fadvise(self.raw_fin.fileno(), next_byte_offset, nbytes, os.POSIX_FADV_WILLNEED)

            # dout has had it's ['DATE'] field converted to float64 and added ['PZERO4'] toit
            # so it can beused in baseline interpolation. Here we convert it back
            if raw_date:
                dtype_out = [
                    (field, np.dtype(">f4")) if field == "DATE" 
                    else (field, self.dtype.fields[field][0])
                    for field in self.dtype.fields
                ]
                dout["DATE"] -= self.hdulist[0].header["PZERO4"]
                dout = dout.astype(dtype_out)
            
            yield dout
    
    def convert_visrows_into_block(self, data, keep_all_baselines = False):
        if keep_all_baselines:
            bl_selector = np.arange(self.raw_nbl).astype('int')
        else:
            bl_selector = self.baseline_indices
        
        #TODO - I should take away the .reshape in fast_raw_blocks() and put it here instead

        dout_data = data['DATA'].transpose((*np.arange(1, data['DATA'].ndim), 0))[bl_selector]
        dout_complex_data = self._create_masked_data(dout_data, 0)
        return dout_complex_data
    
    def convert_block_into_visrows(self, block, add_back_removed_baslines = True):
        
        def compare_shapes(shape1, shape2):
            # Remove dimensions with length 1
            shape1_filtered = tuple(dim for dim in shape1 if dim != 1)
            shape2_filtered = tuple(dim for dim in shape2 if dim != 1)
    
            # Compare shapes
            return shape1_filtered == shape2_filtered
        
        if add_back_removed_baslines:
            nbl = self.raw_nbl
            block_nbl = block.shape[0]
            if block_nbl < nbl:
                bl_indices = self.baseline_indices
            else:
                bl_indices = np.arange(nbl)
        else:
            nbl = self.nbl
            bl_indices = np.arange(nbl)

        nt = block.shape[-1]
        nrows = nbl * nt
        
        dout = np.empty(nrows, dtype=self.dtype)
        
        expected_row_shape = dout.dtype['DATA'].shape[:-1]      #:-1 because the last axis is cmplx
        assert compare_shapes(expected_row_shape, block[0, ..., 0].shape), f"{expected_row_shape}, {block[0, ..., 0].shape}, {block.shape}"

        unity_weigths_array = np.ones(expected_row_shape)
        
        for it in range(nt):
            for ii, ibl in enumerate(bl_indices):
                if add_back_removed_baslines:
                    irow = it*nbl + ibl
                else:
                    irow = it*nbl + ii
                dout[irow]['DATA'][..., 0] = block[ii, ..., it].real.reshape(expected_row_shape)
                dout[irow]['DATA'][..., 1] = block[ii, ..., it].imag.reshape(expected_row_shape)

                if isinstance(block, np.ma.core.MaskedArray):
                    weights = 1 - block[ii, ..., it].mask.reshape(expected_row_shape)
                else:
                    weights = unity_weigths_array
            
                dout[irow]['DATA'][..., 2] = weights
        return dout

    def fast_time_blocks(self, nt, fetch_uvws=False, istart=0, keep_all_baselines = False):
        '''
        Reads raw data from uvfits file as an array and returns a block of nt samples
        :nt: number of times to get
        :fetch_uvws: True if you want to return uvws in teh second parameter
        
        :istart: Sample number to start at. Doesnt have to be a multiple of nt
        '''

        for dout in self.fast_raw_blocks(istart, nt=nt):
            dout_complex_data = self.convert_visrows_into_block(dout, keep_all_baselines = keep_all_baselines)
            #dout_data = dout['DATA'].transpose((*np.arange(1, dout['DATA'].ndim), 0))[self.baseline_indices]
            #dout_complex_data = self._create_masked_data(dout_data, 0)

            if keep_all_baselines:
                bl_selector = np.arange(self.raw_nbl).astype('int')
                bl_order = self.raw_baseline_order
            else:
                bl_selector = self.baseline_indices
                bl_order = self.internal_baseline_order
        
            uvws = []
            if fetch_uvws:
                for it in range(nt):
                    this_uvw = {}

                    for ii, ibl in enumerate(bl_selector):
                        blid = bl_order[ii]

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
        returns a numpy
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
        :trange: tuple (istart,istop) samples
        """
        return craco.time_block_with_uvw_range(
            vis=self.vis, trange=trange, flagant=self.flagant,
            flag_autos=self.ignore_autos, mask=self.mask
        )

    @property
    def tscale(self):
        scale = self.header.get('TSCALE', 'TAI').lower()
        return scale

    def get_tstart(self):
        '''
        return tstart as astropy Time from header otherwise first row of fits table
        '''
        f = self

        # early version of MPIPIPELINE before 10 Nov wrote TAI
        # rather than UTC scale to the timestamps
        # now we write a 'TSCALE' header
        scale = self.tscale
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
    def tend(self):
        last_sample_date = Time(self.end_date, format='jd', scale=self.tscale)
        return last_sample_date

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
        coord = SkyCoord(ra, dec, frame='icrs', equinox='J2000')
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
        try:
            row = self.source_table_entry
            src = row['SOURCE']
        except:
            if 'OBJECT' in f.header:
                src =f.header['OBJECT']
            else:
                src = 'UNKNOWN'
            
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
    
    @property
    def antenna_positions(self):
        '''
        Converts FITS antenna table to a collection suitable for UVW source
        
        :returns: Ordered dict keyed by antenna name and value is positions [3] vector in IERS coordinate frame in meters
        '''
        anttab = self.hdulist['AIPS AN']
        antdata = OrderedDict()
        nant = len(anttab.data)
        
        for iant, row in enumerate(anttab.data):
            antdata[row['ANNAME']] = row['STABXYZ']
        return antdata

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
