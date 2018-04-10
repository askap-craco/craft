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
import vcraft

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

    detect(values.files, values)

bat_cards = ('START_WRITE_BAT40','STOP_WRITE_BAT40','TRIGGER_BAT40')
frame_cards = ('START_WRITE_FRAMEID','STOP_WRITE_FRAMEID','TRIGGER_FRAMEID')

def detect(files, values):
    vfiles = [vcraft.VcraftFile(f) for f in files]
    mux = vcraft.VcraftMux(vfiles)
    infsamp = mux.samp_rate
    # TODO: Calculate frequencies a bit better - tricky because individual
    # files have freqs with gaps, which makes life a little wierd in sigprocland
    fch1 = mux.freqconfig.freq
    foff = mux.freqconfig.bw
    # TODO: Get tstart from BATs. This is the easy way
    tstart = mux.start_mjd
    
    nbits = 32
    nchan = mux.freqconfig.nchan
    nfft = 64
    nguard = 5
    nint = values.nsamps
    nchanout = nfft - 2*nguard
    tsamp = nint*nfft/infsamp

    hdr = {'data_type': 1,
           'tsamp': tsamp,
           'tstart': tstart,
           'fch1':fch1,
           'foff':foff/float(nchanout),
           'nbits':nbits,
           'nifs':1,
           'nchans':nchan*nchanout,
           'src_raj':0.0, # todo: add sigprog.deg2sex(mux.beam_pos[0])
           'src_dej':0.0

    }
    foutname = f.replace('.vcraft','.fil')
    fout = SigprocFile(foutname, 'w', hdr)

    nsamps = mux.nsamps
    nsampout = nsamps/(nint*nfft)
    nsampin = nsampout*nint*nfft

    logging.debug('Writing sigproc header %s', hdr)
    logging.debug(' nsamps=%s nsampin=%s nsapout=%s samp rate=%s',nsamps, nsampin, nsampout, infsamp)

    # reshape to integral number of integrations
    for s in xrange(nsampout):
        sampno = s*nint*nfft
        df = mux.read(sampno, nint*nfft)
        #df.shape = (nint, nfft, nchan)

        for c in xrange(nchan):
            #dc = df[:, :, c]
            dc = df[:, c]
            dc.shape = (-1, nfft)
            dfft = np.fft.fftshift(np.fft.fft(dc, axis=1), axes=1)
            dfft = dfft[:, nguard:nchanout+nguard]
            assert dfft.shape == (nint, nchanout)
            dabs = abs((dfft * np.conj(dfft))).mean(axis=0)
            assert len(dabs) == nchanout
            #dabs -= dabs.mean()
            #dabs /= dabs.std()
            #dabs *= 18*2
            #dabs += 128
            dabs = dabs[::-1]
            dabs= dabs.astype(np.float32)
            dabs.tofile(fout.fin)
            if values.show:
                pylab.plot(dabs)
                pylab.show()




#    dfil = np.add.reduceat(dout, np.arange(0, nsamps, values.nsamps), axis=0)
    # get rid of last sample

    #dfil.tofile(fout.fin)
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

        fig, ax = pylab.subplots(3,3, sharex=True)
        int_times = np.arange(dfil.shape[0])*tsamp
        ax = ax.flatten()
        for f in xrange(nchan):
            ax[f].plot(int_times*1e3, dfil[:, f])
            ax[f].set_xlabel('Time (ms)')
            ax[f].set_ylabel('Power')

        pylab.figure()
        pylab.imshow(dfil.T, aspect='auto', interpolation='None')

        pylab.figure()
        df = np.fft.rfft(dfil, axis=0)
        fsamp = 1./tsamp
        fft_freqs = np.arange(len(df)) *fsamp/2.0/len(df)
        pylab.plot(fft_freqs[1:], abs(df[1:, :])**2)
        pylab.xlabel('Frequency (Hz)')
        pylab.show()




if __name__ == '__main__':
    _main()
