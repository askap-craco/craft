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

def parse_array(line):
    return list(map(float, line.split(' ')[1].split(',')))

def load_dada_beams(f):
    with open(f, 'rU') as fin:
        for line in fin:
            if line.startswith('BEAM_RA'):
                ras = parse_array(line)

            if line.startswith('BEAM_DEC'):
                decs = parse_array(line)

    return list(zip(ras, decs))

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-D', '--dada-header', help='Pck off beam positions from dada header')
    parser.add_argument('-b', '--block-nint', help='Number of integrations to read per block', type=int, default=1024)
    parser.add_argument('-o','--output-prefix', help='Output prefix to use. Default: input filname')
    parser.add_argument('-f','--fifo', help='Make FIFOs instead of files', action='store_true', default=False)
    parser.add_argument(dest='files', nargs=1)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    fname = values.files[0]
    fin = sigproc.SigprocFile(fname)
    nbeams = fin.nifs
    nchans = fin.nchans
    out_files = []
    header = fin.header
    
    if values.dada_header:
        radecs = load_dada_beams(values.dada_header)
    else:
        radecs = [(header['src_raj'], header['src_dej']) for i in range(nbeams)]

    if fin.nbits == 8:
        dtype = np.uint8
        sz = 1
    elif fin.nbits == 32:
        dtype = np.float32
        sz = 4
    elif fin.nbits == 2:
        dtype = np.uint8
        sz = 1
        assert nchans % 4 == 0, 'Cannot handle nbits=2 for non integer number of bytes'
        nchans = nchans / 4 # effective numebr of chnnels
    else:
        raise ValueError('Unknown nbits {}'.format(fin.nbits))
       
    assert len(radecs) == nbeams, 'Number of beams in dada_header {} doesnt match number of beams in filter file {}'.format(len(radec), nifs)
    out_files = []
    outf_names = []

    if values.output_prefix is None:
        prefix = fname
    else:
        prefix = values.output_prefix
    
    for b, radec in enumerate(radecs):
        outf_name = prefix.replace('.fil','') + '_sif{:02d}.fil'.format(b)
        outf_names.append(outf_name)
        if values.fifo:
            logging.debug('Making fifo %s', outf_name)
            outf = os.mkfifo(outf_name)

        outf = open(outf_name, 'w')

            
        # write stupid sigproc header

        header = fin.header
        header['src_raj'] = radec[0]
        header['src_dej'] = radec[1]
        header['nifs'] = 1

        # TODO: Write modified header
        outf.write(fin.hdr)
        out_files.append(outf)
        logging.debug('Writing beam %s radec %s to %s', b, radec, outf_name)

    count = values.block_nint * nchans * nbeams

    while True:
        v = np.fromfile(fin.fin, dtype=dtype, count=count)
        if len(v) == 0:
            logging.debug('EOF')
            break

        assert len(v) <= count

        actual_nelements = len(v) 
        actual_nint = actual_nelements/(nchans*nbeams)
        logging.debug('Read %s elements. Requested %s. nint=%s', len(v), count, actual_nint)

        if actual_nint != values.block_nint:
            new_size = actual_nint*nchans*nbeams
            logging.debug('Undersized block. Probably next EOF.  new size %s.', new_size)

            v = v[0:new_size]
            
        v.shape = (actual_nint, nbeams, nchans)

        for ibeam, outf in enumerate(out_files):
            v[:, ibeam, :].tofile(outf)

    for outf in out_files:
        outf.flush()
        outf.close()
        
    if values.fifo:
        for outf_name in outf_names:
            logging.debug('Removing fifo %s', outf_name)
            os.remove(outf_name)
        

if __name__ == '__main__':
    _main()
