#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import pyfits

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

    h = pyfits.open(values.files[0])

    uvdata = h[0].data
    print(h.info())
    info = []
    for irow in range(len(uvdata)):
        row = uvdata[irow]
        u,v,w,jd,jdfrac,baseline,inttim,freqsel,source,datareal = row
        u,v,w = np.array([u,v,w,])*1e9
        delay = np.sqrt(u**2 + v**2 + w**2)
        jd += jdfrac
        nchan = datareal.shape[3]
        data = np.zeros(nchan, dtype=np.complex)
        data.real = datareal[0,0,0,:,0,0]
        data.imag = datareal[0,0,0,:,0,1]
        weights = datareal[0,0,0,:,0,2]
        baseline = int(baseline)
        ant2 = baseline % 256
        ant1 = baseline / 256
        if ant1 == ant2:
            continue

        datahalf = data[0:nchan/2]
        delayspec = abs(np.fft.fftshift(np.fft.fft(datahalf)))
        peak_delay = np.argmax(delayspec) - float(len(delayspec))/2.
        ang = np.degrees(np.unwrap(np.angle(datahalf)))
        fx = np.arange(len(datahalf))
        phase, gradient = np.polyfit(fx, ang, 1)
        ch12diff = np.degrees(np.angle(datahalf[0:54].mean()) - np.angle(datahalf[54:2*54].mean()))
        info.append((ant1,ant2,u,v,w,delay,peak_delay, phase, gradient, ch12diff))
        plot = False

        if plot:
            fig, (ax1, ax2,ax3,ax4) = pylab.subplots(4,1)
            ax1.plot(abs(data))
            ax2.plot(np.polyval((phase, gradient), fx))
            ax2.plot(np.degrees(np.angle(data)))
            ax3.plot(abs(delayspec))

            fig.suptitle('{}-{} uvw={:.2f},{:.2f},{:.2f} delay={:.2f}'.format(ant1,ant2, u,v,w,delay))
            _info = np.array(info)
            pylab.scatter(_info[:, 5], _info[:, -1])

            pylab.show()

    info = np.array(info)
    pylab.scatter(info[:, 5], info[:, -1])

    pylab.show()



if __name__ == '__main__':
    _main()
