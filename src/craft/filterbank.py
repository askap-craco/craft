#!/usr/bin/env python
"""
Pythonised version of ADE_analysis_syntehsis_WOLA.m from John Tuthill

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import scipy
import scipy.io
import scipy.signal

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def mreshape(x, shape):
    '''Matlab-compatible reshape command. Grrr'''
    return np.reshape(x, shape, order='F')

class PolyphaseFilterbank(object):
    def __init__(self, coeffs, num, den, paths):
        self.c = coeffs
        self.num = num
        self.den = den
        self.paths = paths
        self.K = len(coeffs)
        self.taps = self.K/self.paths # length of each polyphase FIR
        assert self.K == self.paths * self.taps, 'Invalid polyphase setup. K={} paths={} taps={}'.format(self.K, self.paths, self.taps)
        self.os = float(num)/float(den)
        #self.stride = int(self.paths/self.os)
        self.stride = self.paths * self.den / self.num
        self.overlap = paths - self.stride
        self.fft_shifts = np.mod(np.arange(self.num)*(self.paths-self.stride), self.paths);
        self.fft_shifts_rev = self.fft_shifts[::-1]
        
        self.h1 = self.c # analysis filter
        self.h2 = self.c # Syntehsis filter
        self.hh1 = mreshape(self.c[::-1], (self.paths, self.taps)) # reshape analysis filter into polyphase


    def __str__(self):
        s = 'Fiterbank os={}/{} K={} npaths={} ntaps={} stride={}'.format(self.num, self.den, self.K, 
                                                                       self.paths, self.taps, self.stride)
        return s

    __repr__ = __str__

    def analysis(self, x, dtype=None):
        '''
        Does an analysis filterbank with commutation
        :returns: [paths, nframes] array
        '''
        if dtype is None:
            dtype = x.dtype

        N = len(x)
        nframes = (N - self.K)/self.stride
        off = self.K - self.stride
        v1 = np.zeros(self.paths, dtype=dtype)
        afin = np.zeros(self.paths*self.taps, dtype=dtype)
        dout = np.zeros((self.paths, nframes), dtype=np.complex)

        afin[0:off] = x[off:0:-1]
        stride = self.stride
        K = self.K

        for f in range(nframes):
            afin[stride:K] = afin[0:off]
            afin[0:stride] = x[off+(f+1)*stride:off+f*stride:-1]
            reg1 = mreshape(afin, (self.paths, self.taps))
            
            for k in range(self.paths):
                v1[k] = np.dot(reg1[k, :], self.hh1[k, :])

            v1 = v1[::-1]
            shift_indx = self.fft_shifts[f % len(self.fft_shifts)]
            v1 = np.roll(v1, -shift_indx)
            v2 = np.fft.fft(v1)
            dout[:, f] = v2

        return dout

    def synthesis(self, x):
        raise NotImplementedError('Synthesis doesnt work well - you get spikes everywhere from the spike in the coefficients - weird')
        npaths, nframes = x.shape
        if npaths != self.paths:
            raise ValueError('Invalid input')
            
        reg = np.zeros(self.paths * self.taps)
        stride = self.stride
        sf_state = np.zeros(self.K + self.stride*(nframes - 1), dtype=np.complex)
        sf_out = np.zeros((nframes-1)*self.stride, dtype=np.complex)

        for f in range(nframes):
            ind = x[:, f]
            temp = np.fft.ifft(ind)
            shift_indx = self.fft_shifts_rev[f % len(self.fft_shifts_rev)]
            temp = np.roll(temp, -shift_indx)
            for lp1 in range(self.taps): # Could this be implemented as a tile & reshape?
                reg[lp1*self.paths:(lp1+1)*self.paths] = temp 

            reg *= self.h2
            sf_state[f*stride:self.K + f*stride] += reg

            if f > 0:
                sf_out[(f-1)*stride:f*stride] = sf_state[(f-1)*stride:f*stride]

                
            fig, axes = pylab.subplots(4, 1)
            axes = axes.flatten()

            axes[0].plot(sf_state)
            axes[1].plot(sf_out)
            #axes[1].set_xlim(0, 8192)
            axes[2].plot(reg)
            axes[3].plot(temp)

            #axes[0].set_xlim(0, 8192)
            fig.savefig('frame_{:02d}.png'.format(f))
        
            pylab.show()

        return sf_out
      
class AdeFilterbank(PolyphaseFilterbank):
    def __init__(self):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        coeff_file = os.path.join(this_dir, 'ADE_R6_OSFIR.mat')
        c = scipy.io.loadmat(coeff_file)['c'][0, :] 
        # From John Tuthill: the coefficients have rediculous values
        # at the start and end. This is probably a bug in how
        # the coefficients are calculated. Needs research.
        # He reckons it's harmless to set them to zero, so ...
        c[0] = 0
        c[-1] = 0
        PolyphaseFilterbank.__init__(self, c, 32, 27, 1536)

class FftFilterbank(object):
    def __init__(self, N):
        self.N = N

    def analysis(self, x):
        nframes = len(x)/self.N
        xout = np.reshape(x, (nframes, self.N))
        dout = np.fft.fft(xout, axis=1).T
        return dout

class ResampleFilterbank(object):
    def __init__(self, N, num, den, window=('kaiser', 0.5)):
        self.N = N
        self.num = num
        self.den = den
        self.window = window

    def analysis(self, x):
        nframes = len(x)/self.N
        xout = np.reshape(x, (nframes, self.N))
        dout = np.fft.fft(xout, axis=1).T
        for c in range(dout.shape[0]):
            df = scipy.signal.resample_poly(dout[c, :], self.num, self.den, window=self.window)
            start = len(df) - dout.shape[1]
            dout[c, :] = df[start:] # discard first few samples
            
        return dout


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    fb = AdeFilterbank()
    print(fb)
    pylab.figure()
    pylab.plot(fb.c)
    pylab.figure()
    fc = np.fft.fftshift(np.fft.fft(fb.c))
    np.savetxt('fc_fft.txt', abs(fc))
    np.savetxt('coeffs.txt', fb.c)
    pylab.plot(abs(fc))

    snr = 10.
    N = 10*fb.K
    x = np.random.randn(N)*0 + np.cos(2*np.pi*np.arange(N)*21./fb.paths)
    xfb = fb.analysis(x, x.dtype)
    spec = np.real(xfb*np.conj(xfb)).sum(axis=1)
    
    print(xfb.shape, xfb)
    print(spec.shape, spec)
    pylab.figure()
    pylab.plot(spec)
    fb2 = ResampleFilterbank(fb.paths, 32, 27)
    xfb2 = fb2.analysis(x)
    spec2 = np.real(xfb2*np.conj(xfb2)).sum(axis=1)
    pylab.plot(spec2)
    #pylab.xlim(0, 50)

    pylab.figure()
    #xfbb = fb.synthesis(xfb)

    pylab.plot(x)
    #pylab.plot(xfbb)

    pylab.show()
    

if __name__ == '__main__':
    _main()
