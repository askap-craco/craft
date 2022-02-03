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
from .crafthdr import DadaHeader
from .sigproc import SigprocFile
from . import vcraft


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-i','--nsamps', help='Number of samples per integration', type=int, default=1500)
    parser.add_argument('-s','--show', help='Show plots', action='store_true', default=False)
    #parser.add_argument('-n','--nfft', help='FFT size / 64. I.e. for 128 point FFT specify 2', type=int, default=1)
    parser.add_argument('-n','--nchan', help='FFT size / 64. I.e. for 128 point FFT specify 2. '+\
            'If 0 (default) no finer channels are produced', type=int, default=0)
    parser.add_argument('--tstart', type=float, default=0, help='Start time in seconds')
    parser.add_argument('--nsampout', type=int, help='Number of output samples')
    parser.add_argument('-d', '--npolout', help='1=XX+YY, 2=XX,XY,YX,YY', type=int, default=1)
    parser.add_argument('-o','--outfile', help='Outfile', default='beam.fil')
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
    #vfiles = [vcraft.VcraftFile(f) for f in files]
    #mux = vcraft.VcraftMux(vfiles)
    all_muxes = vcraft.mux_by_pol(files, read_ahead=values.nsamps*16)
    npolin = len(all_muxes)
    if npolin == 1:
        npolout = 1
    elif npolin == 2:
        if values.npolout == 1:
            npolout = 1
        else:
            raise NotimplementedError("Only XX+YY currently implemented")
    else:
        assert 'Unexpected number of input polarisations {}'.format(list(all_muxes.keys()))

    pols = sorted(all_muxes.keys())
    # pick random pol to get metadata - todo - check all muxes the same.
    mux = all_muxes[pols[0]] 
    infsamp = mux.samp_rate # samples/sec for a coarse channel
    fch1 = mux.freqconfig.freq
    foff = mux.freqconfig.bw
    tstart = mux.start_mjd
    nbits = 32
    nchan = len(mux.freqs)
 
    nint = values.nsamps
    nsamps_total = mux.nsamps
    nsamps_total = mux.overlap_nsamps

    # Performs a finer channelisation
    if values.nchan > 0:
        nfft = 64*values.nchan
        nguard = 5*values.nchan
        nchanout = nfft - 2*nguard
        bw_file = foff/float(nchanout)
        # Take care of first channel
        fch1 = fch1 - foff/2. + bw_file/2.

        nchan_file = nchan*nchanout
        sample_block = nint*nfft
        tsamp = nint*nfft/infsamp
        

    # Just square and integrate each channel
    elif values.nchan == 0:
        nfft = 0
        nguard = 0
        nchanout = 1 
        bw_file = foff

        nchan_file = nchan
        sample_block = nint
        tsamp = nint/infsamp

    else:
        raise RuntimeError("Wrong --nchans provided: %i. Should be >= 0", values.nchan)

    toff_sec = values.tstart
    toff_insamps = int(np.round(toff_sec * infsamp))
    if values.nsampout is None:
        nsampout = int(nsamps_total/sample_block)
    else:
        nsampout = values.nsampout
        
    nsampin = nsampout*sample_block - toff_insamps

    hdr = {'data_type': 1,
           'tsamp': tsamp,
           'tstart': tstart + toff_sec/86400.,
           'fch1':fch1,
           'foff':bw_file,
           'nbits':nbits,
           'nifs':npolout,
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
    logging.debug(' nsamps=%s nsampin=%s nsapout=%s samp rate=%s',nsamps_total, nsampin, nsampout, infsamp)

    din = np.empty((npolin, sample_block, nchan), dtype=np.complex64)
    dout = np.empty((npolout, nchan_file), dtype=np.float32) # XXX
    # reshape to integral number of integrations
    for s in range(nsampout):
        sampno = s*sample_block + toff_insamps
        logging.debug('block %i/%i, sampno: %i', s, nsampout,sampno)
        for ipol, pol in enumerate(pols):
            din[ipol, :, :]  = all_muxes[pol].read(sampno, sample_block)
            #din[ipol, :, :] = np.zeros((sample_block, nchan))

        # dout = np.empty((npolout, nchan_file), dtype=np.float32) #XXX
            

        for c in range(nchan):
            if values.nchan == 0:
                # We don't have to perform fft
                #for ipol in xrange(npolout):
                #    dout[ipol, c] = np.vdot(din[ipol,:, c], din[ipol, :, c]).real
                dfft = din[:, :, c]
                dfft = dfft[:,:,np.newaxis]


            else:
                dc = din[:, :, c]
                # make blocks of nfft samples
                dc.shape = (npolin, -1 ,nfft)
                dfft = np.fft.fftshift(np.fft.fft(dc, axis=1), axes=1)
                
                # discard guard channels
                dfft = dfft[:, :, nguard:nchanout+nguard]
                invert = True
                if invert:
                    dfft = dfft[:,:,::-1]
                    
                assert dfft.shape == (npolin, nint, nchanout)
            
            # just take the absolute values, average over time and be done with it
            if npolin == 1:
                #dabs = abs((dfft * np.conj(dfft))).mean(axis=0)
                #assert len(dabs) == nchanout
                #dabs = dabs[::-1] # invert spectra
                #dabs= dabs.astype(np.float32) # convert to float32
                #dout[0, :] = dabs
                for thechan in range(dfft.shape[2]):
                    p = dfft[0, :, thechan]
                    dout[0, c*nchanout+thechan] = np.vdot(p,p).real


            elif npolin == 2:
                if values.npolout == 1:
                    for thechan in range(dfft.shape[2]):
                        p1 = dfft[0, :, thechan]
                        p2 = dfft[1, :, thechan]
                        dout[0, c*nchanout+thechan] = np.vdot(p1,p1).real + np.vdot(p2,p2).real

                elif values.npolout == 2: 
                    for chanout in range(nchanout):
                        p1 = dfft[0, chanout, :]
                        p2 = dfft[1, chanout, :]
                        cross = np.vdot(p2, p1)
                        dout[0, chanout] = np.vdot(p1, p1).real
                        dout[1, chanout] = cross.real
                        dout[2, chanout] = cross.imag
                        dout[3, chanout] = np.vdot(p2, p2).real

            if values.show:
                pylab.plot(dabs)
                pylab.show()
                
        dout.flatten().tofile(fout.fin)




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
        print('NBITS', nbits, 'bmax', bmax)
        fig, axes = plt.subplots(1,2, sharex=True, sharey=True)
        for chan in range(nchan):
            axes[0].hist(np.real(df[:, chan]), label='chan %d' % chan, bins=bins, histtype='step')
            axes[1].hist(np.imag(df[:, chan]), label='chan %d' % chan, bins=bins, histtype='step')
        axes[0].legend(frameon=False)
        axes[0].set_xlabel('Real Codeword')
        axes[1].set_xlabel('Imag Codeword')
        axes[0].set_xlim(df.real.min(), df.real.max())

        fig, ax = pylab.subplots(3,3, sharex=True)
        int_times = np.arange(dfil.shape[0])*tsamp
        ax = ax.flatten()
        for f in range(nchan):
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
