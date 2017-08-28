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
from scipy import signal

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
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

    fig, axes = pylab.subplots(4,1)
    d1ax, d2ax,lagax,pax = axes.flatten()
    N = 4096
    Nc = N*512

    c = 0
    d1ax.plot(d1[:N, c].real)
    d1ax.plot(d1[:N, c].imag)
    d2ax.plot(d2[:N, c].real)
    d2ax.plot(d2[:N, c].imag)


    Nf = 128

    assert f1.freqs[c] == f2.freqs[c]
    x1 = d1[:, c].reshape(nsamp/Nf, Nf)
    xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)
    
    x2 = d2[:, c].reshape(nsamp/Nf, Nf)
    xf2 = np.fft.fftshift(np.fft.fft(x2, axis=1), axes=1)
    xx12 = xf1 * np.conj(xf2)
    xx11 = xf1 * np.conj(xf1)
    xx22 = xf2 * np.conj(xf2)

    lagax.plot(abs(xx11.mean(axis=0)), label='auto0')
    lagax.plot(abs(xx22.mean(axis=0)), label='auto1')
    lagax.plot(abs(xx12.mean(axis=0)), label='crossamp')
    lagax.legend(frameon=False)
    pax.plot(np.degrees(np.angle(xx12.mean(axis=0))), 'o')
    pax.set_ylabel('Cross phase (deg)')
    pax.set_xlabel('Channel')
    pylab.show()
    

if __name__ == '__main__':
    _main()
