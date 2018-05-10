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
    parser.add_argument('-i','--nsamps', help='Number of samples per integration', type=int, default=1500)
    parser.add_argument('-s','--show', help='Show plots', action='store_true', default=False)
    parser.add_argument('-d','--dm', help='Coherently dedisperse each channel to DM', type=float, default=None)
    parser.add_argument('-n','--nfft', help='FFT size / 64. I.e. for 128 point FFT specify 2', type=int, default=1)
    parser.add_argument('-o','--outfile', help='Outfile', default=None)
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
    infsamp = mux.samp_rate # samples/sec for a coarse channel
    fch1 = mux.freqconfig.freq
    foff = mux.freqconfig.bw
    tstart = mux.start_mjd
    nbits = 32
    nchan = len(mux.freqs)
    nfft = 64*values.nfft
    nguard = 5*values.nfft
    nint = values.nsamps
    nchanout = nfft - 2*nguard
    finebw = foff/float(nchanout)
    nsamps = mux.nsamps

    if values.dm is None: # just FFT, square and average over nint
        bw_file = finebw
        nchan_file = nchan*nchanout
        sample_block = nint*nfft
        tsamp = nint*nfft/infsamp
        
    else: # coherently dedisperse - FFT, phaseramp, IFFT, then average over nint
        bw_file = foff
        nchan_file = nchan
        sample_block = nfft
        tsamp = nint/infsamp
        assert nfft % nint == 0

    nsampout = nsamps/sample_block
    nsampin = nsampout*sample_block


    hdr = {'data_type': 1,
           'tsamp': tsamp,
           'tstart': tstart,
           'fch1':fch1,
           'foff':bw_file,
           'nbits':nbits,
           'nifs':1,
           'nchans':nchan_file,
           'src_raj':0.0, # todo: add sigprog.deg2sex(mux.beam_pos[0])
           'src_dej':0.0

    }
    if values.outfile:
        foutname = values.outfile
    else:
        foutname = f.replace('.vcraft','.fil')
        
    fout = SigprocFile(foutname, 'w', hdr)

    logging.debug('Writing sigproc header %s', hdr)
    logging.debug(' nsamps=%s nsampin=%s nsapout=%s samp rate=%s',nsamps, nsampin, nsampout, infsamp)

    # prepare phase ramp, if dm is specified
    if values.dm is not None:
        phaseramp = np.empty((nchan, nfft), dtype=np.complex64)
        foffset = (np.arange(nfft) - float(nfft)/2)*finebw # MHz
        print 'FOFFSET', foffset
        
        for ichan, chanfreq in enumerate(mux.freqs):
            freqs = (chanfreq + foffset)*1e-3 # GHz
            delay_ms = 4.15*values.dm*(freqs[0]**-2 - freqs[-1]**-2) # dispersion delay
            delay_us = delay_ms*1e3
            delay_samp = delay_ms *1e-3 * infsamp
            print values.dm, chanfreq, freqs[0], freqs[-1], delay_ms, delay_us,delay_samp
            phaseramp[ichan, :] = np.exp(np.pi*2j*delay_us*foffset)

        if values.show:
            pramps = np.degrees(np.angle(phaseramp.T))
            pylab.plot(pramps, 'o-')
            pylab.show()


    # reshape to integral number of integrations
    for s in xrange(nsampout):
        sampno = s*sample_block
        df = mux.read(sampno, sample_block)
        sampout = np.zeros((sample_block/nint, nchan), dtype=np.float32)

        for c in xrange(nchan):
            dc = df[:, c]
            if values.dm is None:
                dc.shape = (-1, nfft)
                dfft = np.fft.fftshift(np.fft.fft(dc, axis=1), axes=1)

                # discard guard channels
                dfft = dfft[:, nguard:nchanout+nguard]
                assert dfft.shape == (nint, nchanout)

                # just take the absolute values, average over time and be done with it
                dabs = abs((dfft * np.conj(dfft))).mean(axis=0)
                assert len(dabs) == nchanout
                dabs = dabs[::-1] # invert spectra
                dabs= dabs.astype(np.float32) # convert to float32
                dabs.tofile(fout.fin)

            else:
                # do FFT
                assert len(dc) == nfft
                dfft = np.fft.fftshift(np.fft.fft(dc))

                # apply phaseramp
                drotated = (dfft * phaseramp[c, :])
                assert drotated.shape == (nfft,)
                ddelayed = np.fft.fftshift(np.fft.ifft(drotated)) # should be roughly dedispersed
                if values.show:
                    pylab.plot(abs(ddelayed))
                    pylab.show()
                ddelayed.shape = (-1, nint)
                dabs = (abs(ddelayed)**2).mean(axis=1)**2
                sampout[:, c] = dabs


            if values.show:
                pylab.plot(dabs)
                pylab.show()

        if values.dm is not None:
            sampout.tofile(fout.fin)




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
