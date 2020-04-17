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

import fdmt
import numpy as np
from pylab import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
    

class TestFdmt(TestCase):

    def setUp(self):
        self.nf = 336 # number of channels - must be a power of 2 currently.
        self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0 # Channel bandwidth in MHz
        self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.nd = 1024 # Number of DM trials to do
        self.nt = 256 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT

    def tearDown(self):
        pass
    
    def test_initialise_ones(self):
        din = np.ones((self.nf, self.nt))
        fout = self.thefdmt.initialise(din)
        (nc, ndt, nt) = fout.shape
        for idt in xrange(ndt):
            validfout = fout[:, idt, idt:self.nt] # only test out to nt for now
            self.assertTrue(np.all(validfout == float(idt+1)))

    def test_initialise_rand(self):
        din = np.random.randn(self.nf, self.nt)
        fout = self.thefdmt.initialise(din)
        (nc, ndt, nt) = fout.shape
        for idt in xrange(ndt):

            mysum = np.zeros((self.nf, self.nt))

            # The sum for t = sum of din(0->t)
            for i in xrange(idt+1):
                mysum[:, 0:self.nt-i] += din[:, i:]

            #print 'mysum', mysum[0, 0:10]
            #print 'fout', fout[0, idt, 0:10]
            nmax = 10
            diff = fout[:, idt, idt:nmax+idt] - mysum[:, 0:nmax]
            #print 'diff', diff[0, 0:10]
            #print 'diffend', diff[0, -10:-1]
            badidxs= np.where(fout[:, idt, idt:nmax+idt] != mysum[:, 0:nmax])
                      
            self.assertTrue(np.all(abs(diff) < 1e-6))





def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
