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

from . import boxcar
import numpy as np
from pylab import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
    

class TestBoxcar(TestCase):

    def setUp(self):
        self.nd = 10
        self.nt = 21
        self.nbox = 13
        self.boxcar = boxcar.Boxcar(self.nd, self.nbox, np.float32)

    def tearDown(self):
        pass
    
    def test_boxcar_ones(self):
        din = np.ones((self.nd, self.nt), dtype=np.float32)
        fout = self.boxcar(din)
        (nd, nt, nbox) = fout.shape
        self.assertEqual(fout.shape, (self.nd, self.nt, self.nbox))

        # For the first block, it'll be the sum over zeros.
        for t in range(nbox):
            for b in range(nbox):
                if b < t:
                    expect = b + 1
                else:
                    expect = t + 1
                self.assertTrue(np.all(fout[:, t, b] == expect))

        for t in range(nbox, nt):
            for b in range(nbox):
                self.assertTrue(np.all(fout[:, t, b] == (b+1)))

        # Run another block through  - this tiem the whole thing should be right
        # because it kept the history
        fout = self.boxcar(din)

        # the next block should include the history, so it should be all b+1
        for t in range(nt):
            for b in range(nbox):
                self.assertTrue(np.all(fout[:, t, b] == (b+1)))

    def test_boxcar_rand(self):
        din = np.random.randn(self.nd, self.nt).astype(np.float32)
        fout = self.boxcar(din)
        (nd, nt, nbox) = fout.shape
        self.assertEqual(fout.shape, (self.nd, self.nt, self.nbox))



        for t in range(nbox, nt):
            for b in range(nbox):
                myd = din[:, t-b:t+1]
                mysum = myd.sum(axis=1)
                diff = fout[:, t, b] - mysum
                #print 'din', din[0, :]
                #print 'myd', myd[0, :]
                #print 't', t, 'b', b, 'fout', fout[0, t, b], 'mysum',  mysum
                #print 'diff', diff
                self.assertTrue(np.all(abs(diff) < 1e-5))


def davg(x, b):
    return x.sum(axis=0)/(b+1)

def dsum(x, b):
    return x.sum(axis=0)

def dsqrt(x, b):
    return x.sum(axis=0)/np.sqrt(b+1)
    
class TestBoxcarImage(TestCase):

    def setUp(self):
        self.nd = 10
        self.nt = 21
        self.nbox = 13
        self.npix = 4
        self.nt = self.nbox + 5

    def check_box(self, ib, indata, func):
        for d in range(self.nd):
            for t in range(self.nt):
                dout = ib(d, indata[t, :, :]*d)
                for b in range(self.nbox):
                    start = 0 if t < b else t - b
                    expected = func(indata[start:t+1:,:]*d, b)
                    self.assertTrue(np.allclose(dout[:,:,b], expected))
        
        
    def test_avg_impulse(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'avg')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[0, :, :] = 1
        self.check_box(ib, indata, davg)

    def test_avg_step(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'avg')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[:, :, :] = 1
        self.check_box(ib, indata, davg)


    def test_sum_impulse(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'sum')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[0, :, :] = 1
        self.check_box(ib, indata, dsum)

    def test_sum_step(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'sum')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[:, :, :] = 1
        self.check_box(ib, indata, dsum)

    def test_sqrt_impulse(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'sqrt')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[0, :, :] = 1
        self.check_box(ib, indata, dsqrt)

    def test_sqrt_step(self):
        ib = boxcar.ImageBoxcar(self.nd, self.npix, self.nbox, 'sqrt')
        nt, npix = self.nt, self.npix
        indata = np.zeros((nt, npix, npix))
        indata[:, :, :] = 1
        self.check_box(ib, indata, dsqrt)


def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
