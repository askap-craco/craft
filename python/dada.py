#!/usr/bin/env python
"""
Dada file utilities

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import logging
import numpy as np
from crafthdr import DadaHeader

class DadaFile(object):
    def __init__(self, filename):
        self.filename = filename
        assert os.path.isfile(self.filename)
        self.hdr = DadaHeader.fromfile(self.filename)
        self.hdr_size = int(self.hdr.get_value('HDR_SIZE'))
        self.shape = self.hdr.get_value('SHAPE') # throws exception if not defined
        self.shape = np.array(map(int, self.shape.split(',')))
        self.dtype = np.dtype(self.hdr.get_value('DTYPE', '<f4'))
        self.block_size_bytes = np.prod(self.shape)*self.dtype.itemsize
        self.fin = open(self.filename, 'r')

    @property
    def nblocks(self):
        nbytes = os.path.getsize(self.filename)
        nblocks = (nbytes - self.hdr_size) // self.block_size_bytes # truncated

    def get_block(self, blockid):
        if blockid < 0:
            raise ValueError('Invalid blockid {}'.format(blockid))
        
        if blockid >= self.nblocks:
            raise ValueError('BlockID {} past the end of file'.format(blockid))

        assert 0 <= blockid < self.nblocks

        self.fin.seek(self.hdr_size + blockid*self.block_size_bytes)
        count = np.prod(self.shape)
        d = np.fromfile(self.fine, count=count, dtype=self.dtype)
        d.shape = self.shape
        return d

    def blocks(self):
        '''
        Iterates over blocks in the data
        '''
        blockid = 0
        # check number of blocks on every interation in case someone has written
        # to the file since we last looked.
        while blockid < self.nblocks:
            yield self[blockid]

    def __getitem__(self, index):
        return self.get_block(index)

    def __len__(self):
        return self.nblocks

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

