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
import vcraft
import scipy.signal


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Check for sample delays in vcraft files', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help="plot where it's all gone wrong")
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples to compare', default=4096)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    N = values.nsamp
    files = [vcraft.VcraftFile(f) for f in values.files]
    trigger_frameids = [f.trigger_frameid for f in files]
    fstart = max(trigger_frameids)
    f0 = fstart - trigger_frameids[0]
    d0 = files[0].read(f0, N)
    chan = 4
    ncrap = 1440
    for f in files:
        foff = fstart - f.trigger_frameid
        d = f.read(foff, N)
        is_equal = np.all(d[:, :] == d0[:, :])

        print '{} TRIGGER_FRAMEID={} FRAME_OFF={} equal ? {}'.format(f.fname, f.trigger_frameid, foff, is_equal)
        if not is_equal and values.show:
            wrap_frameid = int(f.hdr['START_WRITE_FRAMEID'][0])
            ngood = 96
            print 'WRAP FRAME', wrap_frameid, fstart, f.trigger_frameid, f.trigger_frameid - wrap_frameid, d0.shape, 'ngood', ngood
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



if __name__ == '__main__':
    _main()
