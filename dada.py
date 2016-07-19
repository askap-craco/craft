#!/usr/bin/env python
"""
Dada file utilities

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import logging
import numpy as np

DEFAULT_HEADER_SIZE = 4096

class DadaFile(object):
    def __init__(self, filename):
        self.filename = filename
        self.fin = open(self.filename, 'r')
        self.header_size = DEFAULT_HEADER_SIZE

        # parse twice, in case the first time the header size isnt' enough 
        self.header_size = self._parse_header(self.header_size)
        self.header_size = self._parse_header(self.header_size)
        
    def _parse_header(self, size):
        self.fin.seek(0)
        hdrbytes = self.fin.read(size)
        self.hdr = {}
        for line in hdrbytes.split('\n'):
            bits = line.split(None, 1)
            if len(bits) == 1:
                continue
            
            name, value = line.split(None, 1)
            self.hdr[name.strip()] = value.strip()

        return int(self.hdr['HDR_SIZE'])

    def __getitem__(self, key):
        if isinstance(key, int):
            self.fin.seek(self.header_size+key);
            return self.fin.read(key)
        elif isinstance(key, slice):
            assert slice.step is None, 'Cant do strided access'
            self.fin.seek(self.header_size + key.start)
            sz = key.start - key.stop
            assert sz > 0, 'Cant do negative directions'
            return self.fin.read(sz)

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

