#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-D', '--directory', help='Directory to write results to')
    parser.add_argument('-t','--times', help='samp-start,nsamps comma-separated')
    parser.add_argument('-m','--mjd', help='mjdmiddle,nsec comma-separated')
    
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    mjd_mid = None
    nsec = None
    if values.mjd:
        mjd_mid, nsec = list(map(float, values.mjd.split(',')))
    elif values.times:
        samp_start, nsamps = list(map(int, values.times.split(',')))

    block_nsamps = 1024

    for f in values.files:
        foutname = os.path.join(values.directory, os.path.basename(f))
        assert os.path.abspath(foutname) != os.path.abspath(f), 'Input and output paths are equal: {}={}'.format(os.path.abspath(foutname), os.path.abspath(f))
        fin = sigproc.SigprocFile(f)
        tint = fin.header['tsamp']
        if mjd_mid is not None:
            samp_mid = np.round((mjd_mid - fin.header['tstart'])*86400./tint)
            nsamps = np.round(nsec/tint)
            samp_start = samp_mid - nsamps/2
            if samp_start < 0:
                samp_start = 0

        hdr = fin.header.copy()
        fout_tstart = fin.tstart + fin.tsamp * samp_start / 86400.0
        hdr['tstart'] = fout_tstart
        del hdr['nsamples']
        del hdr['refdm']
        del hdr['fchannel']
        del hdr['period']
        hdr['rawdatafile'] = f
        fout = sigproc.SigprocFile(foutname, header=hdr, mode='w')

        fin.seek_sample(samp_start)
        nblocks = int((nsamps + block_nsamps +11)/block_nsamps)
        block_size = fin.nchans * fin.nifs * fin.nbits*block_nsamps/8
        for blk in range(nblocks):
            b = fin.fin.read(block_size)
            fout.fin.write(b)

        fout.fin.close()
        
        
    

if __name__ == '__main__':
    _main()
