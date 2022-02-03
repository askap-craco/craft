#!/usr/bin/env python2
"""
Checks vcraft files to see whether they have delays or not.

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import vcraft
import scipy.signal


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Check for sample delays in vcraft files', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', type=int, help='Offset samples')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot', default=0)
    parser.add_argument('-n','--fft-size', type=int, help='FFT size per coarse channel', default=128)
    parser.add_argument('-s','--show', action='store_true', help="plot where it's all gone wrong")
    parser.add_argument('-t','--nsamp', type=int, help='Number of samples to compare', default=4096)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    f1 = vcraft.VcraftFile(values.files[0])
    f2 = vcraft.VcraftFile(values.files[1])
    assert np.all(f1.freqs == f2.freqs)


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
            print(h, 'd1=',d1t, 'd2=',d2t, 'diff=',d2t - d1t,mul*(d2t-d1t))

    fig, axes = pylab.subplots(5,1)
    d1ax, d2ax,lagax,pax,lax = axes.flatten()
    N = 4096
    Nc = N*512
    if values.offset is None:
        h = 'TRIGGER_FRAMEID'
        offset = int(f2.hdr[h][0]) - int(f1.hdr[h][0])
    else:
        offset = values.offset


    print('OFFSET IS', offset)
    d1 = f1.read(offset)
    d2 = f2.read(0)
    nsamp, nchan = d1.shape
    # truncate to identical nsamp
    nsamp = min(d1.shape[0], d2.shape[0])
    print('SHAPE BEFORE', d1.shape, d2.shape, nsamp)

    d1 = d1[:nsamp, :]
    d2 = d2[:nsamp, :]
    print('SHAPE AFTER', d1.shape, d2.shape)

    assert d1.shape == d2.shape

    print('Data shape', d1.shape, 'freqs', f1.freqs)
    print('D1 channel0', d1[0:100, 0])
    print('D2 channel0', d2[0:100, 0])


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
    N = values.nsamp
    files = [vcraft.VcraftFile(f) for f in values.files]
    trigger_frameids = [f.trigger_frameid for f in files]
    fstart = max(trigger_frameids)
    f0 = fstart - trigger_frameids[0]
    d0 = files[0].read(f0, N)
    chan = 0
    ncrap = 1440
    for f in files:
        foff = fstart - f.trigger_frameid
        d = f.read(foff, N)
        is_equal = np.all(d[:, :] == d0[:, :])

        print('{} TRIGGER_FRAMEID={} FRAME_OFF={} equal ? {}'.format(f.fname, f.trigger_frameid, foff, is_equal))
        if not is_equal and values.show:
            wrap_frameid = int(f.hdr['START_WRITE_FRAMEID'][0])
            ngood = 96
            print('WRAP FRAME', wrap_frameid, fstart, f.trigger_frameid, f.trigger_frameid - wrap_frameid, d0.shape, 'ngood', ngood)
            d0c = d0[:, chan]
            dc = d[:, chan]

            #extra = d0c[0:ncrap]
            extra = dc[0:ngood]
            extra = np.conj(extra[::-1])
            c1 = scipy.signal.fftconvolve(extra, d0c)
            c2 = scipy.signal.fftconvolve(extra, dc)
            
            fig, ax = pylab.subplots(4,1)
            fig.suptitle(files[0].fname + '-' +f.fname)
            ax[0].plot(d0c.real)
            ax[0].plot(dc.real)
            ax[1].plot(d0c.imag)
            ax[1].plot(dc.imag)
            ax[2].plot((dc-d0c).real)
            ax[2].plot((dc-d0c).imag)
            ax[0].set_ylabel('Real part')
            ax[1].set_ylabel('Imag part')
            ax[2].set_ylabel('Difference (real/imag)')
            ax[3].plot(abs(c1), label=files[0].fname)
            ax[3].plot(abs(c2), label=f.fname)
            ax[3].set_ylabel('Correlation')
            ax[3].legend()

            pylab.show()

    #x1 = d1[offset:shortsamp+offset, c].reshape(-1, Nf)
    x1 = d1[0:shortsamp+0, c].reshape(-1, Nf)

    xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)

    x2 = d2[:shortsamp, c].reshape(-1, Nf)
    xf2 = np.fft.fftshift(np.fft.fft(x2, axis=1), axes=1)
    xx12 = xf1 * np.conj(xf2)
    xx11 = xf1 * np.conj(xf1)
    xx22 = xf2 * np.conj(xf2)

    print('PRODUCT SIZE', xx12.shape, xx12.shape[0]*Nf, nsamp, 'shortsamp', shortsamp, 'offset', offset)
    xx12m = xx12.mean(axis=0)
    punwrap = np.unwrap(np.angle(xx12m))
    xx = np.arange(len(punwrap))
    gradient, phase = np.polyfit(xx, punwrap, 1)
    delay = 32./27.*gradient/2./np.pi*len(punwrap)
    corramp =  abs(xx12m[5:-5]).mean()
    print('Unwrapped phase = {} deg, delay={} us cross amplitude={}'.format(np.degrees(phase), delay, corramp))

    lagax.plot(abs(xx11.mean(axis=0)), label='auto0')
    lagax.plot(abs(xx22.mean(axis=0)), label='auto1')
    lagax.plot(abs(xx12m), label='crossamp')
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
