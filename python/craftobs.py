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
import glob
import heimdall

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def load_beams(path, tstart, ntimes):
    files = sorted(glob.glob(os.path.join(path, '*.fil')))
    if len(files) == 0:
        raise ValueError('No files in path  %s' % path)
    data = None
    for fname in files:
        f = sigproc.SigprocFile(fname)
        tend = tstart + ntimes
        nelements = ntimes*f.nifs*f.nchans
        
        f.seek_data(f.bytes_per_element*tstart)
        if (f.nbits == 8):
            dtype = np.uint8
        elif (f.nbits == 32):
            dtype = np.float32
        
        v = np.fromfile(f.fin, dtype=dtype, count=nelements )
        v.shape = (ntimes, f.nifs, f.nchans)

        if f.nifs == 1:
            nifs = len(files)
            ifslice = slice(0, nifs)
            if data is None:
                data = np.zeros((ntimes, nifs, f.nchans))

            ifnum = int(fname.split('.')[-2])
            data[:, ifnum, :] = v[:, 0, :]

        else:
            data = v


        
    return data


def find_hdrs(root):
    hdrs = []
    for dirpath, dirnames, filenames in os.walk(root):
        h = [os.path.join(dirpath, f) for f in filenames if f.endswith('.hdr')]
        hdrs.extend(h)

    return hdrs


def find_scans(root):
    return [Scan(h) for h in sorted(find_hdrs(root))]

class Scan(object):
    def __init__(self, hdr):
        self.hdrname = hdr
        assert hdr.endswith('.hdr')
        self.root = os.path.dirname(hdr)
    
    def load_candidates(self):
        beam_cands = {}
        logging.debug('Loading candidates from %s', self.root)

        for hpath in glob.glob(os.path.join(self.root, '*.fil.heim')):
            ibeam = int(hpath.split('.')[-3])
            cands = heimdall.load_candidates(hpath)
            logging.debug("loaded candidates %s from path %s", cands.shape, hpath)
            beam_cands[ibeam] = cands

        return beam_cands

    def load_beams(self, start, ntimes):
        return load_beams(self.root, start, ntimes)


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()