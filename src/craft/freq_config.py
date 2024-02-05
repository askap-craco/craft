#!/usr/bin/env python
"""
Class for managing frequency configuration

Copyright (C) CSIRO 2022
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import warnings

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class FrequencyConfig:
    '''
    Represents a frequency setup that handles masked channels. Computes bottom and center
    frequencies and handles fscrunching
    Works in MHz
    '''
    
    def __init__(self, fch1, foff, nchan, chan_mask=None):

        assert fch1 > 0
        assert foff > 0
        assert nchan >= 1

        if chan_mask is None:
            chan_mask = np.zeros(nchan, dtype=bool)

        assert len(chan_mask) == nchan
        
        self._fch1 = float(fch1)
        self._foff = float(foff)
        self._nchan = int(nchan)
        self._mask = chan_mask[:] # defensive copy

    @property
    def fch1(self) -> float:
        '''
        Returns center frequency of first channel
        '''
        return self._fch1

    @property
    def fchtop(self) -> float:
        '''
        Returns center frequency of top channel
        '''
        return self.freq_of_chan(self.nchan - 1)

    @property
    def foff(self) -> float:
        '''
        Returns offset between each channel
        '''
        return self._foff

    @property
    def fcent(self) -> float:
        '''
        Returns center of the band
        '''
        fm = (self.fch1 + self.fchtop)/2.0
        return fm

    @property
    def nchan(self) -> int:
        '''
        Returns total number of channels
        '''
        return self._nchan

    @property
    def channel_frequencies(self):
        f = np.arange(self.nchan)*self.foff + self.fch1
        return f

    @property
    def total_bandwidth(self) -> float:
        return self.nchan * self.foff

    @property
    def nvalid_channels(self) -> int:
        '''
        Returns the number of unmasked channels
        '''
        return self.nchan - self.nmasked_channels

    @property
    def nmasked_channels(self) -> int:
        return sum(self._mask)

    @property
    def channel_mask(self):
        '''
        Returns array of booleans. Nchan long. True=masked (i.e. is invalid)
        Returns a defensive copy
        '''
        return self._mask[:]

    @channel_mask.setter
    def set_channel_mask(self, mask):
        '''
        Set channel mask
        '''
        self._mask = mask[:]

    def mask_channels(self, indexs, maskval=True):
        '''
        Set the mask for the given channel numbers to True
        '''
        self._mask[indexs] = True

    def freq_of_chan(self, ichan):
        assert 0 <= ichan < self.nchan
        
        return self.fch1 + ichan*self.foff

    def __eq__(self, other):
        iseq =  self.fch1 == other.fch1 \
                and self.foff == other.foff \
                and self.nchan == other.nchan \
                and np.all(self.channel_mask == other.channel_mask) 
                
        return iseq


    def __str__(self):
        s = f'FreqConfig fch1={self.fch1} foff={self.foff} nc={self.nchan} nmasked={self.nmasked_channels}'
        return s

    __repr__ = __str__

    def fscrunch(self, factor:int) -> 'FrequencyConfig':
        '''
        Returns a new FreqConfig fscrunched by the given factor. 
        An output channel is marked as masked if *any* input channel that would be scrunched
        into that channel is masked
        '''
        assert factor >= 1, 'Invalid fscrunch factor'
        assert self.nchan % factor == 0, 'Factor must divide into nchan'
        

        new_nchan = self.nchan // factor
        new_foff = factor * self.foff
        new_fch1top = self.freq_of_chan(factor - 1)
        new_fch1 = (self.fch1 + new_fch1top)/2.0
        new_mask = self.channel_mask.reshape(-1, factor).any(axis=1)

        f = FrequencyConfig(new_fch1, new_foff, new_nchan, new_mask)
        return f

    @staticmethod
    def from_freqs_and_masks(freqs, card_masks) -> 'FrequencyConfig':
        '''
        Creates a frequency configuration from the card frequencies and masks from cardcap.py
        freqs should be (nfpga_in_this_header, NCHAN=4) array of channels
        card_masks should be an array of bool the same size as len(freqs)
        :freqs: array of (nfpga_per_file, 4) frequencies. Length = nfiles
        :card_masks: bool array of masks. True = masked. Lenght=nfiles
        '''
        
        freqs = np.array(freqs)
        card_masks = np.array(card_masks)
        assert freqs.shape[0] == card_masks.shape[0], f'Incongruent freqs of card masksfreqs= {freqs.shape} masks={card_masks.shape}'
        nfiles, nfpga_per_file, nchan = freqs.shape
        assert len(card_masks) == nfiles
        assert nchan == 4, 'Weird number of chanenls in a file'
        expected_nchan = freqs.size
        if np.all(card_masks):
            warnings.warn('All with freqs {freqs} were masked. Returning dummy frequency config')
            return FrequencyConfig(1,1,expected_nchan, np.ones(expected_nchan,dtype=bool))


        # reshape to make it have the FPGA shape so we can reshape it back
        # neatly
        if nfpga_per_file == 1:
            freqs = freqs.reshape(-1, 6, 4)
            card_masks = card_masks.reshape(-1,6)
        else:
            card_masks = card_masks[:, np.newaxis]


        mask = np.zeros(freqs.shape, dtype=bool)
        log.info('mask shape %s card_mask shape %s', mask.shape, card_masks.shape)
        mask[:,:,:] = card_masks[:, :, np.newaxis]
        fmask = np.ma.masked_array(freqs, mask=mask)

        fmask_ordered = fmask.transpose([0,2,1]).flatten()
        offsets = fmask_ordered[1:] - fmask_ordered[:-1]

        expected_foff = 0.16666666666666666
        assert np.all(abs(offsets - expected_foff) < 1e-6), f'Something wrong with channel offsets {offsets}'
        
        # You're just faffing around here keith - you know if you freqs and masks you can make
        # In the case that the first card is flagged you cant use that card to work out fch1
        # so you find the first good channnel and work backwards
        first_good_channel_idx = np.where(fmask_ordered.mask == False)[0][0]
        first_good_channel_freq = fmask_ordered[first_good_channel_idx]
        fch1 = first_good_channel_freq - first_good_channel_idx*expected_foff
        total_nchan = len(fmask_ordered)

        config = FrequencyConfig(fch1, expected_foff, total_nchan, fmask.flatten().mask)

        return config
        

    @staticmethod
    def from_cardcap_files(files:list) -> 'FrequencyConfig':
        '''
        Makes a freqconfig from a list of cardcap files
        '''
        all_freqs = [f.frequencies for f in files]
        all_masks = [not f.card_enabled for f in files]
        f =  FrequencyConfig.from_freqs_and_masks(all_freqs, all_masks)
        return f

    @staticmethod
    def from_fits_values(fcent, foff, nchan, chan_mask=None) -> 'FrequencyConfig':
        '''
        Creates a FreqConfig from center frequency, frequency offset and number of channels
        Which are values stored in a typical FITS header
        '''
        fch1 = fcent - (float(nchan) - 1)/2*foff
        return FrequencyConfig(fch1, foff, nchan, chan_mask)

    @staticmethod
    def from_hdu(hdu) -> 'FrequencyConfig':
        '''
        Creates a frequency config from a fits HDU
        Uses CRVAL4, CDELT4 and CRPIX4
        '''
        hdr = hdu.header
        fch1 = hdr['CRVAL4']
        foff = hdr['CDELT4']
        ch1 = hdr['CRPIX4']
        #assert ch1 == 1.0, f'Unexpected initial frequency: {ch1}'
        assert foff > 0, 'cant handle negative frequencies anywhere athe moment foff=%f' % foff
        vis = hdu.data
        nchan = vis[0]["DATA"].shape[-3]
        # Need to add 1 due to stupid FITS convention. Grr.

        freqs = (np.arange(nchan, dtype=np.float64) - ch1 + 1)*foff + fch1 # Hz
        actual_fch1 = fch1 - (ch1 - 1)*foff
        f = FrequencyConfig(actual_fch1/1e6, foff/1e6, nchan)
        assert np.all(abs(freqs/1e6 - f.channel_frequencies) < 1e-6)

        return f


    @staticmethod
    def from_filterbank_values(fch1:float, foff:float, nchan:int, chan_mask=None) -> 'FrequencyConfig':
        ''''
        Creates a FreqConfig from values in a filterbank
        '''
        return FrequencyConfig(fch1, foff, nchan, chan_mask)

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
