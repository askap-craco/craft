#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-D', '--directory', help='Directory to write results to')
    parser.add_argument('-t','--times', help='samp-start,nsamps comma-separated')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    samp_start, nsamps = map(int, values.times.split(','))
    block_nsamps = 1024

    for f in values.files:
        foutname = os.path.join(values.directory, f)
        assert os.path.abspath(foutname) != os.path.abspath(f), 'Input and output paths are equal'
        fin = sigproc.SigprocFile(f)
        print dir(fin)
        print fin.header.keys()
        hdr = fin.header.copy()
        fout_tstart = fin.tstart + fin.tsamp * samp_start
        hdr['tstart'] = fout_tstart
        del hdr['nsamples']
        del hdr['refdm']
        del hdr['fchannel']
        del hdr['period']
        hdr['rawdatafile'] = f
        print hdr
        fout = sigproc.SigprocFile(foutname, header=hdr, mode='w')
        fin.seek_sample(samp_start)
        nblocks = (nsamps + block_nsamps -1)/block_nsamps
        block_size = fin.nchans * fin.nifs * fin.nbits*block_nsamps/8
        for blk in xrange(nblocks):
            b = fin.fin.read(block_size)
            fout.fin.write(b)

        fout.fin.close()
        
        
    

if __name__ == '__main__':
    _main()
