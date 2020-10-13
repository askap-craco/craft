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

import craco
import numpy as np
from pylab import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
    

class TestCracoFdmtTranspose(TestCase):

    def test_tranpose_and_xinping_agree(self):
        # Xinping's comments below
    # // Input order is assumed to be [DM-TIME-UV][DM-TIME-UV]
    #// Each [DM-TIME-UV] block has all DM and all UV, but only half TIME
    #// in1 is attached to the first block and in2 is attached to the second block
    #// The first half TIME has [1,2, 5,6, 9,10 ...] timestamps
    #// The second half TIME has [3,4, 7,8, 11,12 ...] timestamps

        ncu = 2
        nd = 2
        nt = 4
        nuv = 1
        d = np.arange(1, nd*nt*nuv + 1, dtype=np.complex64).reshape(nuv, nd, nt)
        dr = craco.fdmt_transpose(d, ncu=ncu)
        print 'Input', d.real.astype(int).flatten()
        print 'Transposed', dr.real.astype(int).flatten()
        expected = np.array([1, 2, 5, 6, 3, 4, 7, 8], dtype=int)
        print 'expected', expected
        print 'test', expected == dr.real.flatten().astype(int)
        self.assertTrue(np.all(expected == dr.flatten().real.astype(int)))

    def test_tranpose_and_inverse_agree_ncu2(self):
        ncu = 2
        nd = 2
        nt = 4
        nuv = 3
        d = np.arange(1, nd*nt*nuv + 1, dtype=np.complex64).reshape(nuv, nd, nt)
        dr = craco.fdmt_transpose(d, ncu=ncu)
        drr = craco.fdmt_transpose_inv(dr, ncu=ncu)

        #print 'drr', drr.real.astype(int).flatten()
        self.assertTrue(np.all(d == drr))

    def test_tranpose_and_inverse_agree_ncu8(self):
        ncu = 8
        nd = 16
        nt = 32
        nuv = 318
        d = np.arange(1, nd*nt*nuv + 1, dtype=np.complex64).reshape(nuv, nd, nt)
        dr = craco.fdmt_transpose(d, ncu=ncu)
        drr = craco.fdmt_transpose_inv(dr, ncu=ncu)
        self.assertTrue(np.all(d == drr))


    def test_tranpose_and_inverse_agree_flattened_ncu8(self):
        ncu = 8
        nd = 16
        nt = 32
        nuv = 318
        d = np.arange(1, nd*nt*nuv + 1, dtype=np.complex64).reshape(nuv, nd, nt)
        dr = craco.fdmt_transpose(d, ncu=ncu)
        drr = craco.fdmt_transpose_inv(dr.flatten(), ncu=ncu, nt=nt, ndm=nd, nuv=nuv)
        self.assertTrue(np.all(d == drr))



def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()
