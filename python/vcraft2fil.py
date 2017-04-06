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

def detect(f, values):
    hdr = DadaHeader.fromfile(f)
    fin = open(f, 'r')
    hdrsize = int(hdr['HDR_SIZE'][0])
    fin.seek(hdrsize)
    mode = int(hdr['CRAFT_MODE'][0]) & 0x3
    infsamp = float(hdr['SAMP_RATE'][0])
    tsamp = float(values.nsamps)/infsamp
    # TODO: Calculate frequencies a bit better - tricky because individual
    # files have freqs with gaps, which makes life a little wierd in sigprocland
    freqs = map(float, hdr['FREQS'][0].split(','))
    fch1 = min(freqs)
    foff = freqs[1] - freqs[0]
    # TODO: Get tstart from BATs. This is the easy way
    tstart = float(hdr['ANT_MJD'][0])
    hdr = {'data_type': 1,
           'tsamp': tsamp,
           'tstart': tstart,
           'fch1':fch1,
           'foff':foff,
           'nbits':32
    }
    foutname = f.replace('.vcraft','.fil')
    fout = SigprocFile(foutname, 'w', hdr)
    nchan = len(freqs)
           
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
        pylab.ylabel('Voltage')
        pylab.legend()
        
    
        pylab.figure()
        pylab.hist(np.real(df[:, 0]), label='real', bins=128, histtype='step')
        pylab.hist(np.imag(df[:, 0]), label='imag', bins=128, histtype='step')
        pylab.legend()
        pylab.xlabel('Voltage')
        pylab.figure()
        int_times = np.arange(dfil.shape[0])*tsamp
        pylab.plot(int_times, dfil)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Power')
        pylab.figure()
        df = np.fft.rfft(dfil, axis=0)
        fsamp = 1./tsamp
        fft_freqs = np.arange(len(df)) *fsamp/2.0
        pylab.plot(fft_freqs[1:], abs(df[1:, :])**2)
        pylab.xlabel('Frequency (Hz)')
        pylab.show()


        
    
if __name__ == '__main__':
    _main()
