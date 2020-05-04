#!/usr/bin/env python
"""
Boxcar class - for FDMT pipeline testing

Copyright (C) CSIRO 2020
"""
import numpy as np
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Boxcar(object):
    def __init__(self, nd, nbox, dtype=np.float32):
        '''
        Make a boxcar object

        :nd: number of DM trials - i.e. number of independant values to boxcar
        :nbox: Number of boxcar trials
        '''
        assert nd > 0, 'Invalid number nd'
        assert nbox > 0, 'Invalid nbox'
        self.history = np.zeros((nd, nbox), dtype=dtype)
        self.nd = nd
        self.nbox = nbox

    def do_boxcar(self, din):
        nd, nt = din.shape
        assert nd == self.nd, 'Invalid input shape={}'.format(din.shape)
        assert nt >= self.nbox, 'Invalid input shape. Need more nt than nobox'
    
        assert din.dtype == self.history.dtype, 'Invalid input dtype'
        dout = np.zeros((nd, nt, self.nbox), dtype=din.dtype)

        # deliberately not efficient
        # We'll just concatenate the history and the input
        # First nbox samples need to use history
        
        dincat = np.hstack([self.history, din])
        nbox = self.nbox
        assert dincat.shape == (self.nd, nt+self.nbox), 'Got funny shape'

        # boxcar 0 == the input
        dout[:, :, 0] = dincat[:, nbox:nt+nbox]

        for b in xrange(1, nbox):
            dout[:, :, b] = dincat[:, nbox-b:nt+nbox-b] + dout[:, :, b-1]

        self.history[:,:] = din[:, nt-nbox:nt].copy()

        return dout

    __call__ = do_boxcar
        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
