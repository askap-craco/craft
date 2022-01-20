#!/usr/bin/env python
"""
Convert miriad uvspec output to filterbank

Write uvspec with options=avall,tstamp and axis=freq,amp

Copyright (C) CSIRO 2018
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import datetime
import dateutil
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def blockup(filename):
    with open(filename, 'rU') as f:
        d = []
        tstamp = None
        for line in f:
            bits = line.strip().split()
            if len(bits) == 1:
                t = line.split('.')[0].lower()
                tstamp = datetime.datetime.strptime(t, '%y%b%d:%H:%M:%S')
                if len(d) > 0:
                    yield tstamp, np.array(d)
                    d = []
                    tstamp = None
                
            elif len(bits) == 2:
                x, y = list(map(float, bits))
                d.append((x, y))
            else:
                assert False, 'Unepxected bits'
                
        if len(d) > 0:
            yield tstamp, np.array(d, dtype=np.float32)

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs=1)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    outfile = None
    for time, data in blockup(values.files[0]):
        freqs = data[:, 0]
        yvalues = data[:, 1]
        if outfile is None:
            fch1 = freqs[0]*1e3 # MHz # miriad freqs in GHz - sigproc MHz
            foff = (freqs[1] - freqs[0])*1e3 #MHz
            intsamp = 16. # Assume my correlator does 64 point FFTs and integrates 16 of them
            fsamp = 32./27.*1e6 # Hz
            fftsize = 64.
            inttime = 1.0/(fsamp/fftsize/intsamp)
            mjd = 0
            nbits = 32
            nchan = len(freqs)
            hdr = {'fch1':fch1, 'foff':foff,'tsamp':inttime, 'tstart':mjd, 'nbits':32, 'nifs':1, 'nchans':nchan, 'src_raj':0.0, 'src_dej':0.0}
            outf = values.files[0] + '.fil'
            outfile = sigproc.SigprocFile(outf, 'w', hdr)

        yvalues.astype(np.float32).tofile(outfile.fin)

if __name__ == '__main__':
    _main()
