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

from . import fdmt
from . import sample_fdmt
from . import unit_fdmt
import numpy as np
from pylab import *
from .craco import printstats

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

    

class TestSampleFdmtsMatchBlock(TestCase):

    def setUp(self):
        self.nf = 16 # number of channels 
        #self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0e-3 # Channel bandwidth in GHz
        #self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.fmin = 0.976 # freq in GHz
        #self.nd = 1024 # Number of DM trials to do
        self.nd = 16
        #self.nt = 256 # Number of samples per block
        self.nt = 16
        self.tsamp = 0.864 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT - not with history
        self.nbox = 32
        self.doplot = True

    def tearDown(self):
        pass

    def _test_fdmt(self, sampfdmt_class):
        blockin = np.random.randn(self.nf, self.nt)
        goldout = self.thefdmt(blockin)
        sampfdmt = sampfdmt_class(self.thefdmt)
        sampout = sampfdmt(blockin)
        print(goldout.shape, sampout.shape)
        g = goldout[:,:self.nt]
        ok = np.allclose(g, sampout)

        if self.doplot and not ok:
            fig, axs = subplots(1,3)
            printstats(g, 'golden')
            printstats(sampout, 'calculated')
            printstats(g - sampout, 'golden - calculated')
            axs[0].imshow(g)
            axs[1].imshow(sampout)
            axs[2].imshow(g - sampout)
            show()
            
        self.assertTrue(ok)

    def test_max_fifo(self):
        self._test_fdmt(sample_fdmt.MaxFifo)

    def test_max_fifo_per_iteration(self):
        self._test_fdmt(sample_fdmt.MaxFifoPerIteration)

    def test_individual_fifos(self):
        self._test_fdmt(sample_fdmt.IndividualFifos)

    def test_unit_fdmt(self):
        self._test_fdmt(unit_fdmt.UnitFdmt)

def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
