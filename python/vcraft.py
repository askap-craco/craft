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

bat_cards = ('START_WRITE_BAT40','STOP_WRITE_BAT40','TRIGGER_BAT40')
frame_cards = ('START_WRITE_FRAMEID','STOP_WRITE_FRAMEID','TRIGGER_FRAMEID')

class VcraftFile(object):
    def __init__(self, fname):
        self.fname = fname
        f = fname
        hdr = DadaHeader.fromfile(f)
        self.hdr = hdr
        self.fin = open(f, 'r')
        self.hdrsize = int(hdr['HDR_SIZE'][0])
        self.fin.seek(self.hdrsize)
        self.mode = int(hdr['CRAFT_MODE'][0]) & 0x3
        self.nbits = int(hdr['NBITS'][0])
        self.infsamp = float(hdr['SAMP_RATE'][0])
        # TODO: Calculate frequencies a bit better - tricky because individual
        # files have freqs with gaps, which makes life a little wierd in sigprocland
        self.freqs = map(float, hdr['FREQS'][0].split(','))

        # TODO: Get tstart from BATs. This is the easy way
        self.tstart = float(hdr['ANT_MJD'][0])
        self.bats = np.array([int(hdr[c][0], base=16) for c in bat_cards])
        self.frames = np.array([int(hdr[c][0]) for c in frame_cards])
        self.bat0 = int(hdr['NOW_BAT'][0], base=16)
        self.bat0_40 = self.bat0 & 0xffffffffff

    def print_times(self):
        bats = self.bats
        frames = self.frames

        print 'BAT duration (s)', (bats[0] - bats[1])/1e6,'offset', (bats[2] - bats[0])/1e6, 'file offset', (bat0_40 - bats[0])/1e6
        # lowest 32 bits of bat from the time the file was written
        print 'FRAMES duration', (frames[0] - frames[1])/infsamp,'offset', (frames[2] - frames[0])/infsamp

    def read(self):
        ''' Reads everything and returns a numpy array
        with shape dtype=(nsamps,nchan) and dtype=np.complex64
        '''
        # See https://svn.atnf.csiro.au/askapsoft/Src/trunk/Code/Base/adbe/current/adbe/Beamformer.cc
        mode = self.mode
        fin = self.fin
        nchan = len(self.freqs)
        self.fin.seek(self.hdrsize)

        if mode == 0: # 16b+16b
            d = np.fromfile(fin, dtype=np.int16)
            nsamps = len(d)/2/nchan
            assert 2*nchan*nsamps == len(d), 'Not integral number of samples'
            d.shape = (nsamps, nchan, 2)
        elif mode == 1: # 8b + 8b:
            d = np.fromfile(fin, dtype=np.int8)
            nsamps = len(d)/2/nchan
            assert 2*nchan*nsamps == len(d), 'Not integral number of samples'
            d.shape = (nsamps, nchan, 2)
        elif mode == 2: # 4b+4b
            dbytes = np.fromfile(fin, dtype=np.int8)
            nsamps = len(dbytes)/2/nchan*2
            assert len(dbytes) == nsamps*nchan
            d = np.empty((nsamps, nchan, 2), dtype=np.int8)
            dbytes.shape = (nsamps, nchan)
            d[:, :, 0] = (dbytes & 0xf) - (dbytes & 0x8)*2 # real
            dbytes >>= 4
            d[:, :, 1] = (dbytes & 0xf) - (dbytes & 0x8)*2 # imag
        
        elif mode == 3: # 1b+1b
            dbytes = np.fromfile(fin, dtype=np.int8)
            nsamps = len(dbytes)/2/nchan*4
            assert len(dbytes) == nsamps*nchan*4
            d = np.empty((nsamps, nchan, 2), dtype=np.int8)
            dbytes.shape = (nsamps, nchan)
            for t in xrange(4):
                # Convert 0 into +1 and 1 into -1
                d[t::4, :, 0] = 1 - 2*(dbytes & 0x1) # real
                dbytes >>= 1
                d[t::4, :, 1] = 1 - 2*(dbytes & 0x1) # imag
                dbytes >>= 1

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
