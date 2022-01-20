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
from . import fdmt
from . import sample_fdmt

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

MAX_DM_PER_UNIT = 16

class IterConfig(object):
    def __init__(self, unitfdmt, iterno):
        thefdmt = unitfdmt.thefdmt
        ndm = thefdmt.ndm_in_for_iter(iterno)
        nchan = thefdmt.nchan_in_for_iter(iterno)
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
        self.unitfdmt = unitfdmt
        self.maxout = min(self.unitfdmt.maxout, self.ndm_per_unit)

    def ndm_out_for_chan(self, chan):
        # To minimise wastage, this will be channel dependant.
        try:
            ndm = len(self.thefdmt.hist_nf_data[self.iterno][chan][-1])
        except IndexError:
            raise IndexError('No configuraiton for iterno={} chan={}'.format(iterno, chan))

        return ndm

    def ndm_in_for_chan(self, chan):
        # To minimise wastage, this will be channel dependant.
        try:
            ndm = len(self.thefdmt.hist_nf_data[self.iterno-1][chan][-1])
        except IndexError:
            raise IndexError('No configuraiton for iterno={} chan={}'.format(iterno, chan))

        return ndm
        
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


        cid1, cid2, coffset = self.thefdmt.get_config(self.iterno, ochan, odm)
        assert id1 == cid1
        assert id2 == cid2
        assert offset == coffset
        
        return (id1, id2, offset)

    def state2conf(self, d, iterno, input_channel, unitidx, iout):
        iunit = input_channel*self.nunit_per_chan + unitidx
        maxout = self.maxout
        assert 0<= unitidx < self.nunit_per_chan
        assert 0<= iunit < self.total_nunit
        assert 0 <= iout < maxout
        
        iidx = d*maxout + iout # Index of where we are in the whole state
        dmidx = iidx % self.ndm_per_unit # which FIFO DM we'll use
        ioffset = iidx // self.ndm_per_unit # Which offset within that fifo
        read_dm = dmidx + unitidx*self.ndm_per_unit
        #if read_dm < inconfig.ndm_out_for_iter_chan(iterno, output_channel):
        fifo_length = self.unitfdmt.fifo_length(iterno, read_dm, input_channel)
        toffset = fifo_length - 1 - ioffset # time offset from end of FIFO

        print('state2conf d=', d, 'iterno', iterno, 'input_channel', input_channel, 'unitidx', unitidx, 'iout', iout, 'fifo_length', fifo_length, 'toffset', toffset, 'read_dm', read_dm)

        return read_dm, toffset

    def conf2state(self, ind, inchan, toffset):
        iiter = self.iterno
        maxout = self.maxout
        # Units are in Channel, DM order
        iunit = inchan*self.nunit_per_chan + ind // self.ndm_per_unit
        assert 0 <= iunit < self.total_nunit, 'Unexpected iunit {} config={}'.format(iunit, self)

        # the index in [0, NDM_PER_UNIT) for this unit
        dmidx = ind % self.ndm_per_unit
        assert 0 <= dmidx < self.ndm_per_unit

        fifo_length = self.unitfdmt.fifo_length(self.iterno, ind, inchan)

        ioffset = fifo_length - 1 - toffset
        assert 0 <= ioffset <= 3 # Only ever read up to 3 values back from the end of the FIFO

        iout = dmidx % maxout
        assert 0 <= iout < maxout

        dex = ioffset * maxout + dmidx // maxout
        assert 0 <= dex < MAX_DM_PER_UNIT

        print('conf2state', 'ind',ind, 'inchan',inchan, 'toffset', toffset, 'dex',dex, 'iiter', iiter, 'iunit',iunit, 'iout', iout)
        return dex, iiter, iunit, iout
        
    def __str__(self):
        s = 'Iteration {iterno} nchan={nchan} ndm={ndm} ndm_per_unit={ndm_per_unit} nunit_per_chan={nunit_per_chan} nunits={total_nunit}'.format(**vars(self))
        return s
    
    __repr__ = __str__

class UnitFdmtState(object):
    def __init__(self, unitfdmt):
        self.unitfdmt = unitfdmt
        niter = len(self.unitfdmt.thefdmt.hist_nf_data)
        self.state = np.zeros((niter, self.unitfdmt.nunit, self.unitfdmt.maxout))*np.nan
        self.metastate = np.empty((niter, self.unitfdmt.nunit, self.unitfdmt.maxout), dtype=np.object)
                
    def set(self, iterno, iunit, outidx, value, metadata=None):
        self.state[iterno, iunit, outidx] = value
        self.metastate[iterno, iunit, outidx] = metadata
        print('Set state iterno', iterno, 'iunit', iunit, 'outidx', outidx, value, 'meta=', metadata)
        (iterno, read_dm, input_channel, toffset) = metadata
        conf = IterConfig(self.unitfdmt, iterno)
        cdex, ciiter, ciunit, ciout = conf.conf2state(read_dm, input_channel, toffset)
        assert ciiter == iterno, 'Bad set IITER conf={} iterno={}'.format(ciiter, iterno)
        assert ciunit == iunit, 'Bad set UNIT conf={} iunit={}'.format(ciunit, iunit)
        assert ciout == outidx, 'Bad set OUTIDX conf={} outidx={}'.format(ciout, outidx)

    def get(self, iterno, iunit, outidx, metadata=None):
        v = self.state[iterno, iunit, outidx]

        if metadata is not None:
            m = self.metastate[iterno, iunit, outidx]
            assert m == metadata, 'Metadata not equal for iterno={} iunit={} outidx={} expected={} != actual={}'.format(iterno, iunit, outidx, metadata, m)

        return v

class UnitFdmt(sample_fdmt.IndividualFifos):
    '''
    Runs an IndividualFifoFdmt but with a blocke-based parallelisable FDMT execute
    function
    First - need to make better execute.
    Next, need to organise how to setup FIFOs in a more statically-knowable way.
    '''
    def __init__(self, *args, **kwargs):
        super(UnitFdmt, self).__init__(*args, **kwargs)
        self.maxout = 4
        self.nunit = self.thefdmt.n_f


    def fifo_length(self, iterno, dm, channel):
        '''
        Returns size of FIFO - requires IndividualFifos class

        TODO: Make this statically computable - somehow.
        '''
        # Sometimes this FIFO is never used, in which case get_fifo returns KEYERROR
        # and we should return 0
        try:
            fifo_length = len(self._get_fifo(iterno, dm, channel))
        except KeyError:
            fifo_length = -1;

        return fifo_length

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
        nunit = self.nunit # Max number of units is half the numberof channels
        niter = len(thefdmt.hist_nf_data)
        dout = np.zeros(thefdmt.max_dt)

        
        # Push the input data in to the iteration 0 FIFOs
        for c in range(din.shape[0]):
            for d in range(din.shape[1]):
                self.shift(0, d, c, din[c, d])

        for d in range(MAX_DM_PER_UNIT): # "d" is less about DMs and more about a state machine
            # Not all of state is used
            state = UnitFdmtState(self)

            # not all of inputs is used
            inputs = np.ones((niter, nunit, MAX_DM_PER_UNIT))*np.nan

            # Loop through the state and set the outputs appropriately for this d
            for iterno, theiter in enumerate(thefdmt.hist_nf_data):
                inconfig = IterConfig(self, iterno)
                print('d=', d, 'interno', iterno, 'inconfig', inconfig)
                for input_channel in range(inconfig.nchan):
                    for unitidx in range(inconfig.nunit_per_chan):
                        # For each channel we need to update the state, then calculate the
                        iunit = input_channel*inconfig.nunit_per_chan + unitidx
                        # Copy data from FIFOs to state
                        for iout in range(inconfig.maxout):
                            read_dm, toffset = inconfig.state2conf(d, iterno, input_channel, unitidx, iout)
                            if toffset >= 0:
                                v = self.read(iterno, read_dm, input_channel, toffset)
                                state.set(iterno, iunit, iout, v, metadata=(iterno, read_dm, input_channel, toffset))


            for iterno, theiter in enumerate(thefdmt.hist_nf_data):
                inconfig = IterConfig(self, iterno)
                outconfig = IterConfig(self, iterno+1)
                print('Finished writing state. Now calculating inputs to FIFOs', d, iterno, outconfig)
                for output_channel in range(outconfig.nchan):
                    ndm = outconfig.ndm_in_for_chan(output_channel)
                    print('Channel', output_channel, 'ndm', ndm)
                    for unitidx in range(outconfig.nunit_per_chan):
                        # For each channel we need to update the state, then calculate the
                        iunit = output_channel*outconfig.nunit_per_chan + unitidx
                        assert iunit < nunit
                                    
                        for idm in range(outconfig.ndm_per_unit): # For each DM this unit is responsible for
                            odm = unitidx*outconfig.ndm_per_unit + d
                            if odm >= ndm:
                                print('ODM', odm, 'ndm', ndm, 'Breaking out of loop')
                                break

                            # OK so we need to convert odm, output_channel into id1, id2, offset, ichan1 and ichan2
                            in_d1, in_d2, time_offset = inconfig.get_config(output_channel, odm)
                            in_chan1 = 2*output_channel
                            in_chan2 = 2*output_channel+1

                            # then we need to convert that into iterno, iunit, iout, and d
                            dex1, iter1, iunit1, iout1 = inconfig.conf2state(in_d1, in_chan1, 0)
                            dex2, iter2, iunit2, iout2 = inconfig.conf2state(in_d2, in_chan2, time_offset)

                            print('Calc', d, iterno, output_channel, unitidx, odm, end=' ')
                            print(dex1, iter1, iunit1, iout1, in_d1, in_chan1)
                            print(dex2, iter2, iunit2, iout2, in_d2, in_chan2, time_offset)
                            
                            # If current d == expected d, then read value from state
                            if d == dex1:
                                s = state.get(iter1, iunit1, iout1, metadata=(iter1, in_d1, in_chan1, 0))
                                assert not np.isnan(s)
                                inputs[iterno, iunit, d] = s
                                
                            if d == dex2:
                                s = state.get(iter2, iunit2, iout2, metadata=(iter2, in_d2, in_chan2, time_offset))
                                assert not np.isnan(s)
                                inputs[iterno, iunit, d] += s

                            # should have finished everything - do a big old shift
                            if d == MAX_DM_PER_UNIT - 1:
                                vout = inputs[iterno, iunit, d]
                                assert not np.isnan(vout)
                                if iterno == niter - 1: # final iteration write to output
                                    dout[odm] = vout
                                else:
                                    self.shift(iterno+1, odm, output_channel, vout)
                                

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
