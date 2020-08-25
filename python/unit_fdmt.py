#!/usr/bin/env python
"""
Unit based FDMT - meant to be more implementable than the sample FDMT in HLS

Boy am I getting bored of this.


Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import fdmt
import sample_fdmt

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

MAX_DM_PER_UNIT = 16

class IterConfig(object):
    def __init__(self, thefdmt, iterno):
        ndm = thefdmt.ndm_out_for_iter(iterno)
        nchan = thefdmt.nchan_out_for_iter(iterno)
        nunit_per_chan = (ndm + MAX_DM_PER_UNIT - 1) // MAX_DM_PER_UNIT
        ndm_per_unit = (ndm + nunit_per_chan - 1) // nunit_per_chan
        assert ndm_per_unit*nunit_per_chan >= ndm
        total_nunit = nchan * nunit_per_chan
        self.ndm = ndm
        self.nchan = nchan
        self.nunit_per_chan = nunit_per_chan
        self.ndm_per_unit = ndm_per_unit
        self.total_nunit = total_nunit
        self.thefdmt = thefdmt
        self.iterno = iterno


    def ndm_out_for_iter_chan(self, iterno, chan):
        return len(self.thefdmt.hist_nf_data[iterno][chan][-1])
        
    def get_cff(self, ochan):
        iterno = self.iterno
        thefdmt = self.thefdmt
        correction = thefdmt.d_f/2.0
        fres = thefdmt.d_f * (2**(iterno + 1))
        f_start = fres * float(ochan) + thefdmt.f_min
        f_end = f_start + fres
        f_middle = f_start + fres/2.0 - correction
        f_middle_larger = f_middle + 2*correction
        cff1 = fdmt.cff(f_middle, f_start, f_end, f_start)
        cff2 = fdmt.cff(f_middle_larger, f_start, f_end, f_start)
        
        return (cff1, cff2)
        
    def get_config(self, ochan, odm):
        '''
        Calculate the DM config given the iteration number, output channel and output dm
        '''
        cff1, cff2 = self.get_cff(ochan)
        id1 = int(round(odm * cff1))
        offset = int(round(odm * cff2))
        id2 = odm - offset
        
        return (id1, id2, offset)
        
        
    def __str__(self):
        s = 'Iteration {iterno} nchan={nchan} ndm={ndm} ndm_per_unit={ndm_per_unit} nunit_per_chan={nunit_per_chan} nunits={total_nunit}'.format(**vars(self))
        return s
    
    __repr__ = __str__

class UnitFdmt(sample_fdmt.MaxFifoPerIteration):
    '''
    Runs an IndividualFifoFdmt but with a blocke-based parallelisable FDMT execute
    function
    First - need to make better execute.
    Next, need to organise how to setup FIFOs in a more statically-knowable way.
    '''
    def __init__(self, *args, **kwargs):
        super(UnitFdmt, self).__init__(*args, **kwargs)

    def fdmt_process(self, din):
        '''
        Takes in a single initialised time sample and returns the FDMT of it.
        
        :din: np.array of shape (NCHAN, ND_IN)
        :dout: np array of shape (ND)
        '''
        thefdmt = self.thefdmt

        assert din.shape[0] == thefdmt.n_f
        assert din.shape[1] == thefdmt.init_delta_t

        nc = thefdmt.n_f
        nunit = nc/2 # Max number of units is half the numberof channels
        maxout = 4
        niter = len(thefdmt.hist_nf_data)
        dout = np.zeros(thefdmt.max_dt)

        # Not all of state is used
        state = np.zeros((niter, nunit, maxout))
        
        # Push the input data in to the iteration 0 FIFOs
        for c in xrange(din.shape[0]):
            for d in xrange(din.shape[1]):
                self.shift(0, d, c, din[c, d])

        for iterno, theiter in enumerate(thefdmt.hist_nf_data):
            iconfig = IterConfig(self.thefdmt, iterno)
            for output_channel in xrange(iconfig.nchan):
                ndm = iconfig.ndm_out_for_iter_chan(iterno, output_channel)
                chanconfig = thefdmt.hist_nf_data[iterno][output_channel][-1]
                assert len(chanconfig) == ndm
                #for odm, config in enumerate(chanconfig):
                for chanu in xrange(iconfig.nunit_per_chan):
                    for d in xrange(iconfig.ndm_per_unit):
                        odm = chanu*iconfig.ndm_per_unit + d
                        if odm >= ndm: # Because there are extra units, we occasionally overstep the mark.
                            break
                        try:
                            id1, id2, offset = self.thefdmt.get_config(iterno, output_channel, odm)
                            in_d1, in_d2, time_offset = iconfig.get_config(output_channel, odm)
                            assert id1 == in_d1
                            assert id2 == in_d2
                            assert offset == time_offset
                            in_chan1 = 2*output_channel
                            in_chan2 = 2*output_channel+1
                            v1 = self.read(iterno, in_d1, in_chan1, 0)
                            v2 = self.read(iterno, in_d2, in_chan2, time_offset)
                            vout = v1 + v2
                            if iterno == niter - 1: # final iteration write to output
                                dout[odm] = vout
                            else:
                                self.shift(iterno+1, odm, output_channel, vout)
                        except Exception, e:
                            print 'Could not calculate', str(e), iterno, output_channel, chanu, d, odm, in_d1, in_d2, time_offset, in_chan1, in_chan2

        return dout

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
