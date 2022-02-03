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
from .unit_fdmt import *
import numpy as np
from pylab import *
from .craco import printstats

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

    

class TestUnitFdmt(TestCase):

    def setUp(self):
        self.nf = 16 # number of channels 
        #self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0e-3 # Channel bandwidth in GHz
        #self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.fmin = 0.976 # freq in GHz
        #        self.nd = 1024 # Number of DM trials to do
        self.nd = 16
        #self.nt = 256 # Number of samples per block
        self.nt = 16
        self.tsamp = 0.864 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT - not with history
        self.nbox = 32
        self.doplot = True
        self.unitfdmt = UnitFdmt(self.thefdmt)

    def test_state_problem(self):
        d = 1
        iterno = 2
        iunit = 1
        unitidx = 0
        outidx = 2
        value = 0.0
        read_dm = 4
        input_channel = 1
        toffset = 5
        iout = 2

        state = UnitFdmtState(self.unitfdmt)
        metadata = (iterno, read_dm, input_channel, toffset)
        state.set(iterno, iunit, outidx, value, metadata=metadata)

    def test_state2conf_problem(self):
        d = 1
        iterno = 2
        iunit = 1
        unitidx = 0
        outidx = 2
        value = 0.0
        read_dm = 4
        input_channel = 1
        toffset = 5
        iout = 2
        inconfig = IterConfig(self.unitfdmt, iterno)
        outconfig = IterConfig(self.unitfdmt, iterno+1)
        print(inconfig)
        print(outconfig)
        
        read_dm, toffset = inconfig.state2conf(d, iterno, input_channel, unitidx, iout)
        cdex, ciiter, ciunit, ciout = inconfig.conf2state(read_dm, input_channel, toffset)

        self.assertEqual(inconfig.maxout, self.unitfdmt.maxout)
        self.assertEqual(cdex, d)
        self.assertEqual(ciiter, iterno)
        self.assertEqual(ciunit, unitidx)
        self.assertEqual(ciout, iout)



def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
