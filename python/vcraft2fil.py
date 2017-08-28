#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from crafthdr import DadaHeader
from sigproc import SigprocFile

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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
        detect(f, values)

bat_cards = ('START_WRITE_BAT40','STOP_WRITE_BAT40','TRIGGER_BAT40')
frame_cards = ('START_WRITE_FRAMEID','STOP_WRITE_FRAMEID','TRIGGER_FRAMEID')

def detect(f, values):
    hdr = DadaHeader.fromfile(f)
    fin = open(f, 'r')
    hdrsize = int(hdr['HDR_SIZE'][0])
    fin.seek(hdrsize)
    mode = int(hdr['CRAFT_MODE'][0]) & 0x3
    nbits = int(hdr['NBITS'][0])
    infsamp = float(hdr['SAMP_RATE'][0])
    tsamp = float(values.nsamps)/infsamp
    # TODO: Calculate frequencies a bit better - tricky because individual
    # files have freqs with gaps, which makes life a little wierd in sigprocland
    freqs = map(float, hdr['FREQS'][0].split(','))
    fch1 = min(freqs)
    foff = freqs[1] - freqs[0]
    # TODO: Get tstart from BATs. This is the easy way
    tstart = float(hdr['ANT_MJD'][0])
    bats = np.array([int(hdr[c][0], base=16) for c in bat_cards])
    frames = np.array([int(hdr[c][0]) for c in frame_cards])
    bat0 = int(hdr['NOW_BAT'][0], base=16)
    bat0_40 = bat0 & 0xffffffffff

    print 'BAT duration (s)', (bats[0] - bats[1])/1e6,'offset', (bats[2] - bats[0])/1e6, 'file offset', (bat0_40 - bats[0])/1e6
    # lowest 32 bits of bat from the time the file was written
    print 'FRAMES duration', (frames[0] - frames[1])/infsamp,'offset', (frames[2] - frames[0])/infsamp


    hdr = {'data_type': 1,
           'tsamp': tsamp,
           'tstart': tstart,
           'fch1':fch1,
           'foff':foff,
           'nbits':32,
           'nifs':1,
           'nchans':8,
           'src_raj':0.0,
           'src_dej':0.0
           
    }
    foutname = f.replace('.vcraft','.fil')
    fout = SigprocFile(foutname, 'w', hdr)
    nchan = len(freqs)
           
    # See https://svn.atnf.csiro.au/askapsoft/Src/trunk/Code/Base/adbe/current/adbe/Beamformer.cc
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
    
    # detect
    dout = abs(df)**2
    dfil = np.add.reduceat(dout, np.arange(0, nsamps, values.nsamps), axis=0)
    # get rid of last sample
    dfil = dfil[0:-1:, :]
    print 'Input shape', d.shape, 'output shape', dfil.shape
    
    
    dfil.tofile(fout.fin)
    fout.fin.close()

    if values.show:
        pylab.figure()
        pylab.plot(np.real(df[0:4096, 0]), label='real')
        pylab.plot(np.imag(df[0:4096, 0]), label='imag')
        pylab.xlabel('Sample')
        pylab.ylabel('Voltage codeword')
        pylab.legend()
        bmax = 2**(nbits/2 - 1) # Factor of 2 for real+imag and -1 for full width
        bins = np.arange(-bmax-0.5, bmax+0.5, 1)
        print 'NBITS', nbits, 'bmax', bmax
        fig, axes = plt.subplots(1,2, sharex=True, sharey=True)
        for chan in xrange(nchan):
            axes[0].hist(np.real(df[:, chan]), label='chan %d' % chan, bins=bins, histtype='step')
            axes[1].hist(np.imag(df[:, chan]), label='chan %d' % chan, bins=bins, histtype='step')
        axes[0].legend(frameon=False)
        axes[0].set_xlabel('Real Codeword')
        axes[1].set_xlabel('Imag Codeword')
        axes[0].set_xlim(df.real.min(), df.real.max())

        pylab.figure()
        int_times = np.arange(dfil.shape[0])*tsamp
        pylab.plot(int_times, dfil)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Power')
        pylab.figure()
        df = np.fft.rfft(dfil, axis=0)
        fsamp = 1./tsamp
        fft_freqs = np.arange(len(df)) *fsamp/2.0/len(df)
        pylab.plot(fft_freqs[1:], abs(df[1:, :])**2)
        pylab.xlabel('Frequency (Hz)')
        pylab.show()


        
    
if __name__ == '__main__':
    _main()
