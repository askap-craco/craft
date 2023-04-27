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

from craft import craco
from craft import craco_plan
import numpy as np
from craft import uvfits
from pylab import *
from astropy.coordinates import SkyCoord, Angle


__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class TestCoordTransformations(TestCase):
    def test_psitheta_to_coord(self):
        psi = Angle('0.6d') # offset degrees - RA direction
        theta = Angle('0.7d') # offset degrees - dec direction
        phase_center = SkyCoord('19h21m11s +63d12m10s')
        coord = craco.psitheta2coord((psi,theta), phase_center)
        print(coord.to_string('hmsdms'), phase_center.to_string('hmsdms'))

        lm = craco.coord2lm(coord, phase_center)
        expected_lm = np.sin([psi.rad, theta.rad])
        print(lm, expected_lm)
        self.assertAlmostEqual(lm[0], expected_lm[0])
        self.assertAlmostEqual(lm[1], expected_lm[1])
        

class TestCracoBaselineCell(TestCase):
    def test_uvcell(self):
        uvpix = (6, 12)
        npix = 256
        c = craco.BaselineCell(float(1+256*2), uvpix, 10, 19, np.arange(10.0), npix)
        #__init__(self, blid, uvpix, chan_start, chan_end, freqs, npix):
        self.assertTrue(c.is_lower)
        self.assertFalse(c.is_upper)
        self.assertEquals(c.uvpix == uvpix)
        self.assertEquals(c.lower_uvpix == uvpix)
        self.assertEquals(c.upper_uvpix == (npix - 6, npix - 12))


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
        print(('Input', d.real.astype(int).flatten()))
        print(('Transposed', dr.real.astype(int).flatten()))
        expected = np.array([1, 2, 5, 6, 3, 4, 7, 8], dtype=int)
        print(('expected', expected))
        print(('test', expected == dr.real.flatten().astype(int)))
        self.assertTrue(np.all(expected == dr.flatten().real.astype(int)))

    def test_tranpose_and_inverse_agree_ncu2(self):
        ncu = 2
        nd = 2
        nt = 4
        nuv = 3
        d = np.arange(1, nd*nt*nuv + 1, dtype=np.complex64).reshape(nuv, nd, nt)
        dr = craco.fdmt_transpose(d, ncu=ncu)
        drr = craco.fdmt_transpose_inv(dr, ncu=ncu)

        #print ('drr', drr.real.astype(int).flatten())
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


class TestBaseline2Uv(TestCase):

    def setUp(self):
        from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
        parser = ArgumentParser(description='Plans a CRACO scan', formatter_class=ArgumentDefaultsHelpFormatter)
        craco_plan.add_arguments(parser)
        values = parser.parse_args()
        s = '/data/craco/ban115/test_data/frb_d0_t0_a1_sninf_lm100200/frb_d0_t0_a1_sninf_lm100200.fits'
        values.uv = s
        logging.info('Loading UV coordinates from file %s ', values.uv)
        f = uvfits.open(values.uv)
        plan = craco_plan.PipelinePlan(f, values)
        baseline_shape = (plan.nt, plan.nbl, plan.nf)
        uv_shape = (plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide)
        input_data = next(f.time_blocks(plan.nt)) # get first block of data
        print(f'uv_shape={uv_shape}')
        self.input_data = input_data
        self.plan = plan
        self.uv_shape = uv_shape
        f.close()

        
    def test_baseline2uv_works(self):
        dout = np.zeros(self.uv_shape, dtype=np.complex64)
        start = time.time()
        craco.baseline2uv(self.plan, self.input_data, dout)
        end = time.time()
        duration = end - start
        print(f'Baseline2uv took {duration} seconds')

    def test_baseline2uv_and_fast_version_agree(self):
        dout = np.zeros(self.uv_shape, dtype=np.complex64)
        craco.baseline2uv(self.plan, self.input_data, dout)
        fast_version = craco.FastBaseline2Uv(self.plan)
        input_flat = craco.bl2array(self.input_data)
        dout2 = np.zeros_like(dout)
        fast_version(input_flat, dout2)

        start = time.time()
        fast_version(input_flat, dout2)
        end = time.time()
        duration = end - start
        print(f'fast version of Baseline2uv took {duration} seconds')
        self.assertTrue(np.all(dout == dout2))

    def test_conjugate_lower(self):
        dout = np.zeros(self.uv_shape, dtype=np.complex64)
        craco.baseline2uv(self.plan, self.input_data, dout)
        fast_version = craco.FastBaseline2Uv(self.plan, conjugate_lower_uvs=True)
        input_flat = craco.bl2array(self.input_data)
        dout2 = np.zeros_like(dout)
        fast_version(input_flat, dout2)
        start = time.time()
        fast_version(input_flat, dout2)
        end = time.time()
        duration = end - start
        plan = self.plan
        print(f'fast version with conjugation Baseline2uv took {duration} seconds')
        for irun, run in enumerate(plan.fdmt_plan.runs):
            for iuv, uv in enumerate(run.cells):
                if uv.is_lower:
                    self.assertTrue(np.all(dout[irun,:,:,iuv] == np.conj(dout2[irun,:,:,iuv])))
                else:
                    self.assertTrue(np.all(dout[irun,:,:,iuv] == dout2[irun,:,:,iuv]))




def _main():
    unittest_main()
    

if __name__ == '__main__':
    _main()

