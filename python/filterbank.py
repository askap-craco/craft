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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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
        
        self.hh1 = np.reshape(self.c[::-1], (self.paths, self.taps)) # reshape analysis filter into polyphase


    def __str__(self):
        s = 'Fiterbank os={}/{} K={} npaths={} ntaps={} stride={}'.format(self.num, self.den, self.K, 
                                                                       self.paths, self.taps, self.stride)
        return s

    __repr__ = __str__

    def analysis(self, x, dtype=None):
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

        for f in xrange(nframes):
            afin[stride:K] = afin[0:off]
            afin[0:stride] = x[off+(f+1)*stride:off+f*stride:-1]
            reg1 = np.reshape(afin, (self.paths, self.taps))
            
            for k in xrange(self.paths):
                v1[k] = np.dot(reg1[k, :], self.hh1[k, :])

            v1 = v1[::-1]
            shift_indx = self.fft_shifts[f % len(self.fft_shifts)]
            v1 = np.roll(v1, shift_indx)
            v2 = np.fft.fft(v1)
            dout[:, f] = v2

        return dout



class AdeFilterbank(PolyphaseFilterbank):
    def __init__(self):
        this_dir = os.path.dirname(os.path.abspath(__file__))
        coeff_file = os.path.join(this_dir, 'ADE_R6_OSFIR.mat')
        c = scipy.io.loadmat(coeff_file)['c'][0, :] 
        c[-1]  = 0 # for some reason there's a dodgey sample on the end?
        PolyphaseFilterbank.__init__(self, c, 32, 27, 1536)


class FftFilterbank(object):
    def __init__(self, N):
        self.N = N

    def analysis(self, x):
        nframes = len(x)/self.N
        xout = np.reshape(x, (nframes, self.N))
        dout = np.fft.fft(xout, axis=1).T
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
    print fb
    pylab.figure()
    pylab.plot(fb.c)

    snr = 10.
    N = 10*fb.K
    x = np.random.randn(N)*0 + np.cos(2*np.pi*np.arange(N)*22./fb.paths)
    xfb = fb.analysis(x, x.dtype)
    spec = np.real(xfb*np.conj(xfb)).sum(axis=1)
    
    print xfb.shape, xfb
    print spec.shape, spec
    pylab.figure()
    pylab.plot(spec)
    fb2 = FftFilterbank(fb.paths)
    xfb2 = fb2.analysis(x)
    spec2 = np.real(xfb2*np.conj(xfb2)).sum(axis=1)
    pylab.plot(spec2)
    #pylab.xlim(0, 50)
    pylab.show()
    

if __name__ == '__main__':
    _main()
