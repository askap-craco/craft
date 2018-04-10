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
            log.debug('MODE3 %s', ' '.join([startsamp, sampoff, wordidx, nwordsamps, nwords, seek_bytes, nsamps, dwords.shape, d.shape]))
            d = d[sampoff:sampoff+nsamp, :, :]
            assert d.shape[0] == nsamp, 'Incorrect output shape {} expected {}'.format(d.shape, nsamp)
        else:
            raise ValueError('Unsupported mode (yet) {}'.format(mode))


        df = np.empty((nsamps, nchan), dtype=np.complex64)
        df.real = d[:, :, 0]
        df.imag = d[:, :, 1]
        return df

    def print_summary(self):
        d1 = self.read()
        print 'Data shape', d1.shape, 'freqs', self.freqs
        bats = np.array([int(self.hdr[b][0], base=16) for b in bat_cards])
        frames = np.array([int(self.hdr[b][0], base=10) for b in frame_cards])
        for b, bc in zip(bat_cards, bats):
            print b,  self.hdr[b][0], bc, bc-bats[0], (bc-bats[0])/1e6, 's'

        for b, bc in zip(frame_cards, frames):
            print b,  self.hdr[b][0], bc, bc-frames[0], (bc-frames[0])/(1e6*32/27), 's'


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

    for f in values.files:
        print f
        vf = VcraftFile(f)
        vf.print_summary()



if __name__ == '__main__':
    _main()
