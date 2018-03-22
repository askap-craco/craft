#!/usr/bin/env python
"""
Correlate vcraft files

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import vcraft


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', type=int, help='Offset samples')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot', default=0)
    parser.add_argument('-n','--fft-size', type=int, help='FFT size per coarse channel', default=128)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f1 = vcraft.VcraftFile(values.files[0])
    f2 = vcraft.VcraftFile(values.files[1])
    assert f1.freqs == f1.freqs

    d1 = f1.read()
    d2 = f2.read()

    # truncate to identical nsamp
    nsamp = min(d1.shape[0], d2.shape[0])
    print 'SHAPE BEFORE', d1.shape, d2.shape, nsamp

    d1 = d1[:nsamp, :]
    d2 = d2[:nsamp, :]
    print 'SHAPE AFTER', d1.shape, d2.shape

    assert d1.shape == d2.shape

    print 'Data shape', d1.shape, 'freqs', f1.freqs
    for s in ('START_WRITE','STOP_WRITE','TRIGGER'):
        for b in ('FRAMEID','BAT', 'MJD'):
            h = '{}_{}'.format(s,b)
            if b == 'BAT':
                func = lambda x: int(x, base=16)
            elif b == 'FRAMEID':
                func = lambda x: int(x, base=10)
            else:
                func = float
                
            if b == 'MJD':
                mul = 1./86400.
            elif b == 'BAT':
                mul = 1e-6
            elif b == 'FRAMEID':
                mul = 1e-6*32./27.
                
            d1t = func(f1.hdr[h][0])
            d2t = func(f2.hdr[h][0])
            print h, 'd1=',d1t, 'd2=',d2t, 'diff=',d2t - d1t,mul*(d2t-d1t)
            
    nsamp, nchan = d1.shape
    fig, axes = pylab.subplots(5,1)
    d1ax, d2ax,lagax,pax,lax = axes.flatten()
    N = 4096
    Nc = N*512
    if values.offset is None:
        h = 'TRIGGER_FRAMEID'
        offset = int(f2.hdr[h][0]) - int(f1.hdr[h][0])
    else:
        offset = values.offset


    if offset >= d1.shape[0] or offset >= d2.shape[0]:
        raise ValueError('Requested offset is larger than the number of samples. Offset={}, nsamp1={} nsamp2={}'.format(offset, d1.shape[0], d2.shape[0]))

    c = values.channel
    d1ax.plot(d1[:N, c].real, label='real')
    d1ax.plot(d1[:N, c].imag, label='imag')
    d1ax.legend()
    d2ax.plot(d2[:N, c].real)
    d2ax.plot(d2[:N, c].imag)

    Nf = values.fft_size
    shortsamp = ((nsamp-offset)/Nf)*Nf

    assert f1.freqs[c] == f2.freqs[c]
    x1 = d1[offset:shortsamp+offset, c].reshape(-1, Nf)
    xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)

    x2 = d2[:shortsamp, c].reshape(-1, Nf)
    xf2 = np.fft.fftshift(np.fft.fft(x2, axis=1), axes=1)
    xx12 = xf1 * np.conj(xf2)
    xx11 = xf1 * np.conj(xf1)
    xx22 = xf2 * np.conj(xf2)

    print 'PRODUCT SIZE', xx12.shape, xx12.shape[0]*Nf, nsamp, 'shortsamp', shortsamp, 'offset', offset
    punwrap = np.unwrap(np.angle(xx12.mean(axis=0)))
    xx = np.arange(len(punwrap))
    gradient, phase = np.polyfit(xx, punwrap, 1)
    delay = 32./27.*gradient/2./np.pi*len(punwrap)
    print 'Unwrapped phase = {} rad, graidnet={} rad per channel, delay={} us'.format(phase, gradient, delay)

    lagax.plot(abs(xx11.mean(axis=0)), label='auto0')
    lagax.plot(abs(xx22.mean(axis=0)), label='auto1')
    lagax.plot(abs(xx12.mean(axis=0)), label='crossamp')
    lagax.legend(frameon=False)
    pax.plot(np.degrees(np.angle(xx12.mean(axis=0))), 'o')
    pax.set_ylabel('Cross phase (deg)')
    pax.set_xlabel('Channel')
    lax.plot(np.fft.fftshift(abs(np.fft.fft(xx12.mean(axis=0)))), label='lag')

    Navg = 128
    Nout = xx12.shape[0]/Navg
    xx12 = xx12[:Nout*Navg, :]
    xx12.shape = [Nout, Navg, -1 ]
    xx12a = xx12.mean(axis=1)

    pylab.figure()
    pylab.imshow(np.angle(xx12a), aspect='auto')

    pylab.show()
    

if __name__ == '__main__':
    _main()
