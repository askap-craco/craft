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
    

class TestFdmtSimple(TestCase):

    def setUp(self):
        self.nf = 336 # number of channels - must be a power of 2 currently.
        self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0 # Channel bandwidth in MHz
        self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.nd = 1024 # Number of DM trials to do
        self.nt = 256 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT
        self.nbox = 32

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

    def test_effective_variance_calcs_equal(self):
        for idt in xrange(self.nd):
            for w in xrange(self.nbox):
                sig1 = self.thefdmt.get_eff_sigma(idt, w+1)
                var1 = self.thefdmt.get_eff_var_recursive(idt, w+1)
                self.assertLess(abs(sig1 -  np.sqrt(var1)), 1e-6,
                                  'Variances not equal idt={} w={} {}!={}'.format(idt, w, sig1, np.sqrt(var1)))


class TestFdmtSingleBlockHits(TestCase):

    def setUp(self):
        self.nf = 336 # number of channels - must be a power of 2 currently.
        self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0 # Channel bandwidth in MHz
        self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.nd = 1024 # Number of DM trials to do
        self.nt = 256 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT
        self.nbox = 32

    def tearDown(self):
        pass

    def test_self_made_frb_nt(self):
        idt = self.nt-1
        # Check we don't get any extra
        d = np.zeros((self.nf, self.nt))+1 
        frb = self.thefdmt.add_frb_track(idt, d)
        frbsum = (frb.sum() - np.prod(frb.shape))*2
        frbout = self.thefdmt(frb)
        maxpos = np.argmax(frbout)
        maxd, maxt = np.unravel_index(maxpos, frbout.shape)
        self.assertEqual(frbsum, frbout.max(), 'Didnt get all the hits')
        self.assertEqual(maxt, idt, 'Peak at Wrong time')
        self.assertEqual(maxd, idt, 'Peak at Wrong DM')


    def test_self_made_frbs_le_nt(self):
        # Only tests up to nt sized frbs. After that we'll need to test overlap and sum
        for idt in xrange(self.nt):
            d = np.zeros((self.nf, self.nt))+1 
            frb = self.thefdmt.add_frb_track(idt, d)
            frbsum = (frb.sum() - np.prod(frb.shape))*2
            frbout = self.thefdmt(frb)
            maxpos = np.argmax(frbout)
            maxd, maxt = np.unravel_index(maxpos, frbout.shape)
            self.assertEqual(frbsum, frbout.max(), 'Didnt get all the hits {} {} at {} {}'
                             .format(frbsum, frbout.max(), maxd, maxt))
            self.assertEqual(maxt, idt, 'Peak at Wrong time')
            self.assertEqual(maxd, idt, 'Peak at Wrong DM')


class TestFdmtWithHistoryHits(TestCase):

    def setUp(self):
        self.nf = 336 # number of channels - must be a power of 2 currently.
        self.fmax = 1448. +0.5#  Freuency of the top of the band in MHz
        self.df = 1.0 # Channel bandwidth in MHz
        self.fmin = self.fmax - self.nf*self.df # Frequency of the bottom of the band in MHz
        self.nd = 1024 # Number of DM trials to do
        self.nt = 256 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt, history_dtype=np.float32) # make FDMT
        self.nbox = 32

    def tearDown(self):
        pass

    def test_init_with_history(self):
        nt = self.nt
        ramp = np.zeros((self.nf, self.nt*2))
        r = np.arange(0, self.nt*2)
        ramp[:, :] = r[np.newaxis, :]
        f1 = self.thefdmt.initialise(ramp[:, 0:nt])
        f2 = self.thefdmt.initialise(ramp[:, nt:])
        initdt = f1.shape[1]
        for idt in xrange(initdt):
            # f2 at t = 0 should = f1 at nt-1 plus idt + 1
            self.assertLess((f1[:, idt, nt-1] + idt + 1 - f2[:, idt, 0]).max(), 1e-6)
                
    def _test_self_made_frb(self, idt, testwithones=False):

        # If testwithones=1 it make sure you've only added where the FRB is
        # and not elsewhere. The FRB is 2 where the FBR exists and 1 elsewhere
        osum = fdmt.OverlapAndSum(self.nd, self.nt)
        d = np.zeros((self.nf, self.nd), dtype=np.float32)
        if testwithones:
            d += 1
            
        frb = self.thefdmt.add_frb_track(idt, d)
        if testwithones:
            frbsum = (frb.sum() - np.prod(frb.shape))*2
        else:
            frbsum = frb.sum()
            
        nblocks = self.nd/self.nt
        expected_blk = idt//self.nt
        expected_t = idt % self.nt
        total_sum = 0
        
        for blk in xrange(nblocks):
            din = d[:, blk*self.nt:(blk+1)*self.nt]
            fdmtout = self.thefdmt(din)
            frbout = osum(fdmtout)

            
            maxpos = np.argmax(frbout)
            maxd, maxt = np.unravel_index(maxpos, frbout.shape)
            if blk == expected_blk:
                self.assertEqual(frbsum, frbout.max(),
                                 'Didnt get all the hits. idt={} blk={} maxd={} maxt={} frbsum={} frboutmax={}'\
                                 .format(idt, blk, maxd, maxt, frb.sum(), frbout.max()))
                self.assertEqual(maxt, expected_t, 'Peak at Wrong time')
                self.assertEqual(maxd, idt, 'Peak at Wrong DM')

    def test_self_made_frbs_eq_2nt(self):
        # test FRBs up to ND.
        # Requires overlap and sum

        idt = self.nt*2-1
        self._test_self_made_frb(idt)


    def test_self_made_frbs_eq_nd(self):
        # test FRBs up to ND.
        # Requires overlap and sum

        idt = self.nd-1
        self._test_self_made_frb(idt)

    def test_self_made_frbs_eq_2nt_ones(self):
        # test FRBs up to ND.
        # Requires overlap and sum

        idt = self.nt*2-1
        self._test_self_made_frb(idt, True)


    def test_self_made_frbs_eq_nd_ones(self):
        # test FRBs up to ND.
        # Requires overlap and sum

        idt = self.nd-1
        self._test_self_made_frb(idt, True)


    def test_self_made_frbs_le_nd(self):
        # test FRBs up to ND.
        # Requires overlap and sum
        for idt in xrange(self.nd):
            self._test_self_made_frb(idt)


def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
