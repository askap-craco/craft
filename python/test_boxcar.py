#!/usr/bin/env python
"""
Template for testing

Copyright (C) CSIRO 2019
"""
import logging
from unittest import TestCase, main as unittest_main
try:
    from unittest.mock import Mock, MagicMock # -- python3.3 and greater
except ImportError:
    from mock import Mock, MagicMock

import boxcar
import numpy as np
from pylab import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
    

class TestBoxcar(TestCase):

    def setUp(self):
        self.nd = 1
        self.nt = 10
        self.nbox = 3
        self.boxcar = boxcar.Boxcar(self.nd, self.nbox, np.float32)

    def tearDown(self):
        pass
    
    def test_boxcar_ones(self):
        din = np.ones((self.nd, self.nt), dtype=np.float32)
        fout = self.boxcar(din)
        (nd, nt, nbox) = fout.shape
        self.assertEqual(fout.shape, (self.nd, self.nt, self.nbox))

        # For the first block, it'll be the sum over zeros.
        for t in xrange(nbox):
            for b in xrange(nbox):
                if b < t:
                    expect = b + 1
                else:
                    expect = t + 1
                self.assertTrue(np.all(fout[:, t, b] == expect))

        for t in xrange(nbox, nt):
            for b in xrange(nbox):
                self.assertTrue(np.all(fout[:, t, b] == (b+1)))

        # Run another block through  - this tiem the whole thing should be right
        # because it kept the history
        fout = self.boxcar(din)

        # the next block should include the history, so it should be all b+1
        for t in xrange(nt):
            for b in xrange(nbox):
                self.assertTrue(np.all(fout[:, t, b] == (b+1)))

    def test_boxcar_rand(self):
        din = np.random.randn(self.nd, self.nt).astype(np.float32)
        fout = self.boxcar(din)
        (nd, nt, nbox) = fout.shape
        self.assertEqual(fout.shape, (self.nd, self.nt, self.nbox))



        for t in xrange(nbox, nt):
            for b in xrange(nbox):
                myd = din[:, t-b:t+1]
                mysum = myd.sum(axis=1)
                diff = fout[:, t, b] - mysum
                #print 'din', din[0, :]
                #print 'myd', myd[0, :]
                #print 't', t, 'b', b, 'fout', fout[0, t, b], 'mysum',  mysum
                #print 'diff', diff
                self.assertTrue(np.all(abs(diff) < 1e-5))

        
        


def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
