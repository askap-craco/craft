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


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Check for sample delays in vcraft files', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help="plot where it's all gone wrong")
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    N = 4096*4
    files = [vcraft.VcraftFile(f) for f in values.files]
    trigger_frameids = [f.trigger_frameid for f in files]
    fstart = max(trigger_frameids)
    f0 = fstart - trigger_frameids[0]
    d0 = files[0].read(f0, N)
    for f in files:
        foff = fstart - f.trigger_frameid
        d = f.read(foff, N)
        is_equal = np.all(d[:, :] == d0[:, :])

        print '{} TRIGGER_FRAMEID={} FRAME_OFF={} equal ? {}'.format(f.fname, f.trigger_frameid, foff, is_equal)
        if not is_equal and values.show:
            fig, ax = pylab.subplots(3,1)
            fig.suptitle(files[0].fname + '-' +f.fname)
            ax[0].plot(d0.real)
            ax[0].plot(d.real)
            ax[1].plot(d0.imag)
            ax[1].plot(d.imag)
            ax[2].plot((d-d0).real)
            ax[2].plot((d-d0).imag)
            ax[0].set_ylabel('Real part')
            ax[1].set_ylabel('Imag part')
            ax[2].set_ylabel('Difference (real/imag)')
            pylab.show()



if __name__ == '__main__':
    _main()
