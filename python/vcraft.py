#!/usr/bin/env python
"""
Abstracts away VCRAFT files

Copyright (C) CSIRO 2017
"""
__author__ = "Keith Bannister <keith.bannister@csiro.au>"

import numpy as np
import os
import sys
import logging
from crafthdr import DadaHeader
import freqconfig
log = logging.getLogger(__name__)

bat_cards = ('START_WRITE_BAT40','STOP_WRITE_BAT40','TRIGGER_BAT40')
frame_cards = ('START_WRITE_FRAMEID','STOP_WRITE_FRAMEID','TRIGGER_FRAMEID')

# number of samples per 32 bit word indexed by mode
SAMPS_PER_WORD32 = [1,2,4,16,1,2,4,16]

def unpack_craft(fin, nsamp_per_word):
    dwords = np.fromfile(fin, dtype=np.uint32)
    nwords = len(dwords)/nchan
    nsamps = nwords*nsamp_per_word
    assert len(dwords) == nchan*nwords
    dwords.shape = (nwords, nchan)

    d = np.empty((nsamps, nchan, 2), dtype=np.int8)
    bitshift = 32/nsamp_per_word
    assert nsamp_per_word*bitshift == 32

    for samp in xrange(nsamp_per_word):
        # convert from 4-bit two's complement to int8
        d[samp::4, :, 0] = (dwords & 0xf) - (dwords & 0x8)*2 # real
        dwords >> 4
        d[samp::4, :, 1] = (dwords & 0xf) - (dwords & 0x8)*2 # imag
        dwords >> 4

class VcraftFile(object):
    def __init__(self, fname, mode=None):
        self.fname = fname
        f = fname
        hdr = DadaHeader.fromfile(f)
        self.hdr = hdr
        self.fin = open(f, 'r')
        self.hdrsize = int(hdr['HDR_SIZE'][0])
        self.fin.seek(self.hdrsize)
        if mode is None:
            self.mode = int(hdr['CRAFT_MODE'][0]) & 0x3
        else:
            self.mode = mode

        self.nbits = int(hdr['NBITS'][0])
        self.infsamp = float(hdr['SAMP_RATE'][0])
        # TODO: Calculate frequencies a bit better - tricky because individual
        # files have freqs with gaps, which makes life a little wierd in sigprocland
        self.freqs = np.array(map(float, hdr['FREQS'][0].split(',')))

        # TODO: Get tstart from BATs. This is the easy way
        self.tstart = float(hdr['ANT_MJD'][0])
        self.bats = np.array([int(hdr[c][0], base=16) for c in bat_cards])
        self.frames = np.array([int(hdr[c][0]) for c in frame_cards])
        self.bat0 = int(hdr['NOW_BAT'][0], base=16)
        self.bat0_40 = self.bat0 & 0xffffffffff


    @property
    def nsamps(self):
        ''' Return the number of complex samples in the file'''

        file_bytes= os.path.getsize(self.fname)
        data_bytes = file_bytes - self.hdrsize
        nchan = len(self.freqs)
        smp_per_word = SAMPS_PER_WORD32[self.mode]
        data_words = data_bytes/4

        nsamp = data_words*smp_per_word/nchan

        return nsamp

    def print_times(self):
        bats = self.bats
        frames = self.frames

        print 'BAT duration (s)', (bats[0] - bats[1])/1e6,'offset', (bats[2] - bats[0])/1e6, 'file offset', (bat0_40 - bats[0])/1e6
        # lowest 32 bits of bat from the time the file was written
        print 'FRAMES duration', (frames[0] - frames[1])/infsamp,'offset', (frames[2] - frames[0])/infsamp

    def read(self, startsamp=0, nsamp=None):
        ''' Reads everything and returns a numpy array
        with shape dtype=(nsamps,nchan) and dtype=np.complex64
        '''
        # See https://svn.atnf.csiro.au/askapsoft/Src/trunk/Code/Base/adbe/current/adbe/Beamformer.cc
        mode = self.mode
        fin = self.fin
        nchan = len(self.freqs)
        self.fin.seek(self.hdrsize)
        assert startsamp >= 0, 'Invalid startsamp {}'.format(startsamp)
        if nsamp is None:
            nsamp = self.nsamps - startsamp

        assert startsamp + nsamp <= self.nsamps, 'Asked for too many samples. startsamp {} nsamp {} nsamps {}'.format(startsamp, nsamp, self.nsamps)

        if mode == 0: # 16b+16b
            self.fin.seek(self.hdrsize + startsamp*4*nchan)
            d = np.fromfile(fin, dtype=np.int16, count=nsamp*nchan*2)
            # truncate if required
            #d = d[:2*nchan*nsamps]
            assert 2*nchan*nsamp == len(d), 'Not integral number of samples'
            d.shape = (nsamp, nchan, 2)
            nsamps = nsamp
        elif mode == 1: # 8b + 8b:
            if startsamp != 0:
                raise NotImplementedError('Give me a sec!')

            # This works perfectly, but we need some practice at the other method with non-native types
            #d = np.fromfile(fin, dtype=np.int8)
            #print 'mode1 = 8+8', d[0:6]
            #nsamps = len(d)/2/nchan
            #assert 2*nchan*nsamps == len(d), 'Not integral number of samples'
            # each 32 bit word contains 4x8 bit numbers imag/real/imag/real is two seuential samplesfrom a single channel
            #d.shape =(nsamps/2, nchan, 2, 2)
            # let's put those axes next to each other
            #d = d.transpose(0,2,1,3)
            # and reshape it to the orrect shape
            #d = d.reshape(nsamps, nchan, 2)
            dwords = np.fromfile(fin, dtype=np.uint32)
            # each word contains 4 8 bit numbers, (imag/real)*2
            nwords = len(dwords)/nchan
            assert len(dwords) == nchan*nwords
            dwords.shape = nwords, nchan

            nsamps = nwords*2
            d = np.empty((nsamps, nchan, 2), dtype=np.int8)
            for samp in xrange(2):
                # convert from 4-bit two's complement to int8
                d[samp::2, :, 0] = (dwords & 0xff) - (dwords & 0x80)*2 # real
                dwords >> 8
                d[samp::2, :, 1] = (dwords & 0xff) - (dwords & 0x80)*2 # imag
                dwords >> 8


        elif mode == 2: # 4b+4b
            if startsamp != 0:
                raise NotImplementedError('Give me a sec!')
            # each 32 bit word contais 8 4 bit numbers (imag/real)*4 for the same channel
            dwords = np.fromfile(fin, dtype=np.uint32)
            nwords = len(dwords)/nchan
            nsamps = nwords*4
            assert len(dwords) == nchan*nwords
            dwords.shape = (nwords, nchan)

            d = np.empty((nsamps, nchan, 2), dtype=np.int8)

            for samp in xrange(4):
                # convert from 4-bit two's complement to int8
                d[samp::4, :, 0] = (dwords & 0xf) - (dwords & 0x8)*2 # real
                dwords >> 4
                d[samp::4, :, 1] = (dwords & 0xf) - (dwords & 0x8)*2 # imag
                dwords >> 4

        elif mode == 3: # 1b+1b
            # each 32 bit word contais 32 1 bit numbers (imag/real)*16 for the same channel
            wordidx = startsamp / 16 # which 32 bit word the start sample is in
            sampoff = startsamp % 16 # how may samples into the first word the start sample is
            nwordsamps = (nsamp + 15 + sampoff) / 16 # how many times we need to read (in words)
            nwords = nwordsamps*nchan # total numberof words including channels
            seek_bytes = self.hdrsize + wordidx*nchan*4  # seek offset in bytes
            fin.seek(seek_bytes)
            dwords = np.fromfile(fin, dtype='<u4', count=nwords)
            #nwords = len(dwords)/nchan
            assert len(dwords) == nwords
            dwords.shape = (nwordsamps, nchan)
            nsamps = nsamp
            d = np.empty((nwordsamps*16, nchan, 2), dtype=np.int8)
            for samp in xrange(16):
                # Convert 0 into +1 and 1 into -1
                d[samp::16, :, 0] = 1 - 2*(dwords & 0x1) # real
                dwords >>= 1
                d[samp::16, :, 1] = 1 - 2*(dwords & 0x1) # imag
                dwords >>= 1
            log.debug('MODE3 startsamp=%s sampoff=%s wordidx=%s nwordssamps=%s nwords=%s seek_bytes=%s nsamps=%s dwords.shape=%s d.shape=%s',
                      startsamp, sampoff, wordidx, nwordsamps, nwords, seek_bytes, nsamps, dwords.shape, d.shape)
            d = d[sampoff:sampoff+nsamp, :, :]
            assert d.shape[0] == nsamp, 'Incorrect output shape {} expected {}'.format(d.shape, nsamp)
        else:
            raise ValueError('Unsupported mode (yet) {}'.format(mode))


        df = np.empty((nsamps, nchan), dtype=np.complex64)
        df.real = d[:, :, 0]
        df.imag = d[:, :, 1]
        return df

    def print_summary(self):
        bats = np.array([int(self.hdr[b][0], base=16) for b in bat_cards])
        frames = np.array([int(self.hdr[b][0], base=10) for b in frame_cards])
        for b, bc in zip(bat_cards, bats):
            print b,  self.hdr[b][0], bc, bc-bats[0], (bc-bats[0])/1e6, 's'

        for b, bc in zip(frame_cards, frames):
            print b,  self.hdr[b][0], bc, bc-frames[0], (bc-frames[0])/(1e6*32/27), 's'


class VcraftMux(object):
    '''
    Multiplexes together VCRAFT files to give a contiguous, monotonically increasing
    frequency axis
    '''

    def __init__(self, vcraft_files):
        '''
        :vcraft_files: A list of open Vcraft files
        '''
        self._files = vcraft_files
        self.ant = self.hdr_identical('ANT')
        self.beam = int(self.hdr_identical('BEAM'))
        self.nbits = int(self.hdr_identical('NBITS'))
        self.mode = int(self.hdr_identical('MODE'))
        self.data_type = self.hdr_identical('DATA_TYPE') # just a check
        assert self.data_type == 'CRAFT_VOLTAGES'
        self.samp_rate = float(self.hdr_identical('SAMP_RATE'))
        freq_str = self.allhdr('FREQS')
        freqs = np.array([map(float, flist.split(',')) for flist in freq_str])
        nchan_per_file = len(freqs[0])
        assert freqs.shape == (len(self._files), nchan_per_file)
        
        self.freqconfig = freqconfig.FreqConfig(freqs, reverse=True)
        self.freqs = np.arange(self.freqconfig.nchan_span)*self.freqconfig.bw + self.freqconfig.freq

        self.all_samps = [f.nsamps for f in self._files]
        self.nsamps = min(self.all_samps)
        self.trigger_frameids = np.array(map(int, self.allhdr('TRIGGER_FRAMEID')))
        self.trigger_mjds = np.array(map(float, self.allhdr('TRIGGER_MJD')))
        self.sample_offsets = self.trigger_frameids - min(self.trigger_frameids)
        assert np.all(self.sample_offsets >= 0)
        assert np.all(self.sample_offsets < self.nsamps)
        self.start_mjd = self.trigger_mjds[np.argmin(self.sample_offsets)]
        beam_ra = np.array(map(float, self.allhdr('BEAM_RA'))).mean() # TODO: make craft_vdump write same BEAM_RA for all card/fpgas
        beam_dec = np.array(map(float, self.allhdr('BEAM_DEC'))).mean()
        #self.beam_pos = SkyCoord(beam_ra, beam_dec, frame='icrs', unit=('deg','deg'))
        self.beam_pos = (beam_ra, beam_dec)


    def allhdr(self, cardname):
        '''
        Returns a list of header values for all files, with the given card name
        '''
        header_values = [f.hdr[cardname][0] for f in self._files]
        return header_values

    def hdr_identical(self, cardname):
        '''
        Returns the header value for the given cardname. 
        Checks the value is the same for all files. If not, it throws a ValueError
        '''
        hdr_values = set(self.allhdr(cardname))
        if len(hdr_values) != 1:
            raise ValueError('Exepected the same header value for {} for all files. Got these values {}'.format(cardname, hdr_values))

        return hdr_values.pop()

    def read(self, samp_start=0, nsamp=None):
        '''
        Read into a giant buffer - demuxing along the way.
        Probably not memory optimal, but it'll have to do for now.
        Something with 'yield' in it would probably make more sense - I should do that
        '''
        if nsamp is None:
            nsamp = self.nsamps

        assert samp_start + nsamp <= self.nsamps, 'Invalid read request. nsamp={} samp_start ={} differnce{}'.format(nsamp, samp_start, nsamp-samp_start)

        # allocate giant buffer
        d = np.empty((nsamp, len(self.freqs)), dtype=np.complex64)
        for ifile, f in enumerate(self._files):
            out_chans = self.freqconfig.chanmaps[ifile, :]
            fsamp_start = samp_start + self.sample_offsets[ifile]
            d[:, out_chans] = f.read(fsamp_start, nsamp)

        return d
         

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-n','--nsamps', help='Number of samples per integration', type=int, default=1500)
    parser.add_argument('-s','--show', help='Show plots', action='store_true', default=False)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    all_files = []
    for f in values.files:
        print f
        vf = VcraftFile(f)
        vf.print_summary()
        all_files.append(vf)

    mux = VcraftMux(all_files)
    print mux.freqconfig
    print mux.freqconfig.freqs
    print mux.freqconfig.chanmaps
    print mux.freqconfig.freqmaps
    print mux.freqs

        

    



if __name__ == '__main__':
    _main()
