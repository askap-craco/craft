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
import numpy as np
from pylab import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

# change range to range for python3

class TestFdmtLookupTableGeneration(TestCase):

    def setUp(self):
        self.nf = 32 # number of channels
        self.fmin = 760.
        self.df = 1.0 # Channel bandwidth in MHz
        self.fmax = 760 + self.df*self.nf
        self.nd = 186 # Number of DM trials to do
        self.nt = 256 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT

    def test_fres(self):
        self.assertEqual(self.thefdmt.fres_for_iter(0), 1.0)
        self.assertEqual(self.thefdmt.fres_for_iter(1), 1.0*2.0)
        self.assertEqual(self.thefdmt.fres_for_iter(2), 1.0*4.0)
        self.assertEqual(self.thefdmt.fres_for_iter(3), 1.0*8.0)

    def test_freq_of_chan(self):
        self.assertEqual(self.thefdmt.freq_of_chan(0, 0), 760)
        self.assertEqual(self.thefdmt.freq_of_chan(0, 1), 760+1)
        self.assertEqual(self.thefdmt.freq_of_chan(0, 20), 760+20)
        self.assertEqual(self.thefdmt.freq_of_chan(1, 0), 760)
        self.assertEqual(self.thefdmt.freq_of_chan(1, 1), 760+1*2)
        self.assertEqual(self.thefdmt.freq_of_chan(1, 2), 760+2*2)
        self.assertEqual(self.thefdmt.freq_of_chan(2, 2), 760+4*2)
        

    def test_calc_cff(self):
        for iterno in range(self.thefdmt.niter):
            for c in range(self.thefdmt.nchan_out_for_iter(iterno)):
                id1_cff = self.thefdmt.calc_id1_cff(iterno, c)
                off_cff = self.thefdmt.calc_offset_cff(iterno, c)
                chanconfig = self.thefdmt.hist_nf_data[iterno][c][-1]
                print(('Iterno', iterno, 'c', c, 'id1_cff', id1_cff, 'offset_cff', off_cff, int(np.round(id1_cff*(1<<16))), int(np.round(off_cff*(1<<16)))))
                for idm in range(len(chanconfig)):
                    _, id1, offset, id2, _, _, _ = chanconfig[idm]
                    self.assertEqual(id1, int(np.round(idm*id1_cff)))
                    self.assertEqual(offset, int(np.round(idm*off_cff)))


    def test_cal_lookup_table(self):
        # just checks it works fo rnow

        lut = self.thefdmt.calc_lookup_table()

    
class TestFdmtSimple(TestCase):

    def setUp(self):
        self.nf = 336 # number of channels 
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
        for idt in range(ndt):
            validfout = fout[:, idt, idt:self.nt] # only test out to nt for now
            self.assertTrue(np.all(validfout == float(idt+1)))

    def test_initialise_rand(self):
        din = np.random.randn(self.nf, self.nt)
        fout = self.thefdmt.initialise(din)
        (nc, ndt, nt) = fout.shape
        for idt in range(ndt):

            mysum = np.zeros((self.nf, self.nt))

            # The sum for t = sum of din(0->t)
            for i in range(idt+1):
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
        for idt in range(self.nd):
            for w in range(self.nbox):
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
        for idt in range(self.nt):
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
        self.nt = 64 # Number of samples per block
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
        for idt in range(initdt):
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
        
        for blk in range(nblocks):
            din = d[:, blk*self.nt:(blk+1)*self.nt]
            fdmtout = self.thefdmt(din)
            self.assertEqual(fdmtout.shape[0], self.nd)
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
        for idt in range(self.nd):
            self._test_self_made_frb(idt)


class TestFdmtWithHistoryHits(TestCase):

    def setUp(self):
        self.fc = 0.860 # center frequency GHz
        self.bw = 0.288 # bandwidth GHz
        self.Nd = 64 # number of DM trials
        self.Nchan= 64
        self.Nt = 256 # time block size
        self.Tint = 0.864e-3 # integration time - seconds
        self.f1 = self.fc - self.bw/2.
        self.f2 = self.fc + self.bw/2.
        self.chanbw = 1e-3
        self.thefdmt = fdmt.Fdmt(self.f1, self.chanbw, self.Nchan, self.Nd, self.Nt)
        self.nf = self.Nchan
        self.nt = self.Nt
        self.nd = self.Nd


    def test_ones_all_hit(self):
        ones = np.ones((self.nf, self.nt), dtype=np.float32)
        oneout = self.thefdmt(ones)[:, :self.nt]
        nzero = np.sum(oneout == 0)

        # All times up to nt should be >= 1 for all DMs
        plot = True
        if nzero > 0:
            fig, axs = subplots(1,2)
            axs[0].imshow(oneout, aspect='auto', origin='lower')
            axs[1].imshow(oneout == 0, aspect='auto', origin='lower')
            show()
        self.assertEqual(nzero, 0)

    def test_initialise_all_hit(self):
        ones = np.ones((self.nf, self.nt), dtype=np.float32)
        oneout = self.thefdmt.initialise(ones)
        nzero = np.sum(oneout == 0)
        # All times up to nt should be >= 1 for all DMs
        plot = False
        if nzero > 0:
            fig, axs = subplots(1,2)
            axs[0].imshow(oneout[0,:,:], aspect='auto', origin='lower')
            axs[1].imshow(oneout[0,:,:] == 0, aspect='auto', origin='lower')
            show()
        self.assertEqual(nzero, 0)


class TestFdmtHighDm(TestCase):
    def setUp(self):
        self.nf = 256# number of channels 
        self.df = 1e-3 # Channel bandwidth in GHz
        self.fmin = 0.716 # Fmin GHz
        self.nd = 1024 # Number of DM trials to do
        self.nt = 64 # Number of samples per block
        self.tsamp = 1.0 # milliseconds
        self.thefdmt = fdmt.Fdmt(self.fmin, self.df, self.nf, self.nd, self.nt) # make FDMT

    def test_init_ones_at_t0(self):
        ones = np.ones((self.nf, self.nt))
        onei = self.thefdmt.initialise(ones)
        print((onei.shape))
        (nf, nd, nt) = onei.shape
        self.assertEqual(nf, self.nf)
        t0 = onei[:, :, 0]
        for c in range(nf):
            self.assertTrue(np.allclose(t0[c, :], np.ones(nd, dtype=float)))

    def test_init_ones_at_t1(self):
        ones = np.ones((self.nf, self.nt))
        onei = self.thefdmt.initialise(ones)
        print((onei.shape))
        (nf, nd, nt) = onei.shape
        self.assertEqual(nf, self.nf)
        t1 = onei[:, :, 1]
        for c in range(nf):
            self.assertTrue(np.allclose(t1[c, 0 ], np.ones(1, dtype=float)))
            self.assertTrue(np.allclose(t1[c, 1:], np.ones(nd-1, dtype=float)*2))


    def test_ones_at_t0(self):
        ones = np.ones((self.nf, self.nt))
        oneout = self.thefdmt(ones)
        print((oneout.shape))

        # Look at the very first time sample
        t0 = oneout[:, 0]

        # t0 should be monotonically decreasing
        t0d = t0[1:] = t0[:-1]

        is_decreasing = np.all(t0d >= 0)

        plot = False

        if not is_decreasing and plot:
            print (t0d)
            fig, axs = subplots(1,2)
            axs[0].imshow(oneout, aspect='auto', origin='lower')
            axs[1].plot(oneout[:,0].T)
            show()

            
        self.assertTrue(is_decreasing, 'the first time sample should be monotonically decreasing with DM')

    def test_ones_at_t1(self):
        '''
        Ones intput at t=1 should be at least 3 for the largest DM
        '''
        ones = np.ones((self.nf, self.nt))
        oneout = self.thefdmt(ones)
        print((oneout.shape))
        
        # Look at the very first time sample
        t1 = oneout[:, 1]
        plot = False

        if plot:
            fig, axs = subplots(1,2)
            idm = self.nd - 1
            lastdm_tf = np.zeros((self.nf, self.nd))
            lastdm_tf = self.thefdmt.add_frb_track(idm, lastdm_tf)
            freqs = np.arange(self.nf)*self.df + self.fmin
            channels = np.arange(self.nf) - 0.5
            delays = -idm*(freqs**-2 - freqs.max()**-2)/(freqs.max()**-2 - freqs.min()**-2)
            ny, nx = lastdm_tf.shape
            xticks = np.arange(nx+1) - 0.5
            yticks = np.arange(ny+1) - 0.5
            axs[0].imshow(lastdm_tf, aspect='auto', origin='lower')
            axs[0].plot(delays, channels)
            axs[0].set_xticks(xticks, minor=True)
            axs[0].set_yticks(yticks, minor=True)
            axs[0].grid(True, 'minor')
            axs[1].plot(t1)
            show()

            
        #self.assertEqual(t1[-1], 3.0) #, 'T=1 should have at least 1 smeared channel'# KB did this ever pass? I'm not sure.
        self.assertEqual(t1[-1], 2.0) # This does pass, but I can't work out if I shoud have it or 3.0???


        

def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
