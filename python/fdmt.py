#!/usr/bin/env python
"""
FDMT class - borrowed from Barak Zackay

Copyright (C) CSIRO 2017
"""
import numpy as np
import logging
from numba import njit, jit, prange
import math

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

@jit(nopython=True)
def isquare(f):
    return 1./(f*f)

@jit(nopython=True)
def cff(f1_start, f1_end, f2_start, f2_end):
    return (isquare(f1_start) - isquare(f1_end))/(isquare(f2_start) - isquare(f2_end))

@jit(nopython=True)
def calc_delta_t(f_min, f_max, f_start, f_end, max_dt):
    rf = cff(f_start, f_end, f_min, f_max)
    delta_tf = (float(max_dt) - 1.0)*rf
    delta_t = int(math.ceil(delta_tf)) + 1
    
    return delta_t


def copy_kernel(v1, vout):
    '''
    Pretty simple really - just copies v1 to the output

    v1 and vout have differnet shapes - so we write the final part with zero
    '''
    assert len(vout) >= len(v1)
    n = len(v1)
    vout[:n] = v1[:]
    vout[n:] = 0

    assert not np.any(np.isnan(vout))

def add_offset_kernel2(v1, v2, vout, toff):
    '''
    Add vector v1 to v2 applying an offset of toff samples, save result in vout.
    If shift is larger enough to go past the end of vout, then ignore the stuff past the end of vout.
    In the FMDT, toff is always small enough that there's at least 1 sample of overlap
    
    Raise error if: len(v1) != len(v2) or len(vout) < len(v1) or toff < 0 or toff > len(v1)

    Implementation notes: There's a few max/min and conditions in here that could be done in the plan
    Saving the odd branch or two in the kernel, if you care.

    This is the key kernel of the FDMT.
    :v1: input vector. Always the same size as v2
    :v2: Input vector. Always the same size as v1
    :vout: Output vector. Always same or larger size than v1/v2.
    :toff: Number of samples to shift v2 by. Always >= 0.

    >>> add_offset_kernel2(np.arange(4), np.arange(4)+1, np.zeros(6, dtype=int), 0)
    array([1, 3, 5, 7, 0, 0])

    >>> add_offset_kernel2(np.arange(4), np.arange(4)+1, np.zeros(6, dtype=int), 1)
    array([0, 2, 4, 6, 4, 0])

    >>> add_offset_kernel2(np.arange(4), np.arange(4)+1, np.zeros(6, dtype=int), 2)
    array([0, 1, 3, 5, 3, 4])

    >>> add_offset_kernel2(np.arange(4), np.arange(4)+1, np.zeros(6, dtype=int), 3)
    array([0, 1, 2, 4, 2, 3])

    >>>  # Can't do larger offsets - add_offset_kernel2(np.arange(4), np.arange(4)+1, np.zeros(6, dtype=int), 4)
    #array([0, 1, 3, 5, 3, 4])
    
    '''
    assert len(v1) == len(v2)
    assert toff >= 0
    assert len(vout) >= len(v1)
    assert not np.any(np.isnan(v1)), 'V1 contains nans'
    assert not np.any(np.isnan(v2)), 'V2 contains nans'

    nt_in = len(v1)
    nt_out = len(vout)
    nsum = nt_in - toff
    if toff >= nt_in:
        raise NotImplementedError('Cant do disjoint offsets yet - dont need to anyway')

    t = 0
    vout[0:toff] = v1[0:toff]
    t += toff
    vout[t:t+nsum] = v1[t:t+nsum] + v2[0:nsum]
    t+= nsum
    nrest = min(nt_in - nsum, nt_out - t)

    if nrest > 0:
        vout[t:t+nrest] = v2[nsum:nsum+nrest]
        t += nrest

    # this part where we write zeros isn't really necessary if we didn't want to look at a pretty square plot at the end. 
    # Ideally you wouldn't waste memory, or memory bandwidth on this. But I really like looking at pretty plots.
    if t < nt_out:
        vout[t:nt_out+1] = 0

    assert not np.any(np.isnan(vout)), 'Data not written nt_in={} nt_out={} nsum={} toff={} {}'.format(nt_in, nt_out, nsum, toff, vout)
                 
    return vout

class Fdmt(object):
    def __init__(self, f_min, f_off, n_f, max_dt, n_t):
        '''
        Make an object that does the FDMT.
        This contsructor makes the FDMT plan which can then be called with .execute() or just __call__

        :f_min: minimum frequency in MHz
        :f_off: Maximum frequency in MHz
        :n_f: Numberof channels
        :max_dt: Number of DM trials
        :n_t: Number of samples in input block
        '''
        self.f_min = float(f_min)
        self.d_f = f_off
        self.n_f = int(n_f)
        assert n_f > 0
        self.bw = self.n_f * f_off
        self.f_max = self.f_min + self.bw
        assert(self.f_min < self.f_max)
        self.niter = int(np.ceil(np.log2(n_f)))
        self.max_dt = int(max_dt)
        self.n_t = int(n_t)
        freqs = np.arange(n_f)*self.d_f + self.f_min
        self.init_delta_t = self._calc_delta_t(self.f_min, self.f_min + self.d_f)
        self._state_shape = np.array([self.n_f, self.init_delta_t, self.n_t + self.init_delta_t])
        self.hist_delta_t = [self.init_delta_t]
        self.hist_state_shape = [self._state_shape]

        # This hist_nf_data is the guts of the plan. I need to organise this to be more than just a list of tuples.
        # 
        self.hist_nf_data = []

        # channel width of top (copied) and bottom (i.e. all remaining) channels
        self._df_top = self.d_f
        self._df_bot = self.d_f
        self._ndt_top = self.init_delta_t

        # make the plan by pretend running through the iterations
        for i in xrange(1, self.niter+1):
            self._save_iteration(i)

    def _calc_delta_t(self, f_start, f_end):
        return calc_delta_t(self.f_min, self.f_max, f_start, f_end, self.max_dt)

    @jit(nopython=False, parallel=True)
    def _save_iteration(self, intnum):
        '''
        Appends the iteration description to self.hist_nf_data

        If there are some quirks in here it's partly becuase I'm trying to keep closely to Barak's original code
        and I'm a little wary about changing everything otherwise it's tricky to debug.

        Currently can only handle non power of 2

        :intnum: iteration number. starts at 1 according to Barak's convention.
        '''
        assert intnum > 0
        
        n_t = self.n_t
        df = self.d_f
        s = self.hist_state_shape[intnum-1] # input state shape
        nf_in = s[0]
        nf = nf_in//2 + nf_in % 2 # output number of channels - Includes copied channel, if required
        do_copy = nf_in % 2 == 1 # True if we have an odd number of input channels and the top one will be copied
        fjumps = float(nf) # output number of subbands

        if do_copy:
            pass # top channel width unchanged
        else:
            self._df_top += self._df_bot # Top channel will be wider by the new channel

        self._df_bot *= 2.0 # Bottom channels will be added together

        if nf == 1: # if this is the last iteration
            delta_f = self._df_top
        else:
            delta_f = self._df_bot

        fres = self._df_bot
        # delta_f = 2**(intnum)*self.d_f # channel width in MHz - of the normal channels
        delta_t = self._calc_delta_t(self.f_min, self.f_min + delta_f) # Max IDT for this iteration
        ndt = delta_t

        state_shape = np.array([nf, delta_t, n_t + ndt])
        
        correction = 0.0
        if intnum > 0: # this is always invoked - it's a leftover from Barak's code
            correction = self.d_f/2.0*0

        # shift input and shift output are never used - they're leftovers from barak's code
        shift_input = 0
        shift_output = 0

        # keep delta_t, state_shape, and per_subband parameters for posterity
        self.hist_delta_t.append(delta_t)
        self.hist_state_shape.append(state_shape)
        nf_data = []
        self.hist_nf_data.append(nf_data)

        #print 'Iteration = {} output bottom channel bandwidth = {} bottom channel output bandwidth={} top channel input bandwidth{} copy? {} inshape={} outshape={}'.format(intnum, fres, self._df_bot, self._df_top, do_copy,
        #s, state_shape)

        # for each ouput subband
        for iif in xrange(nf):
            is_top_subband = iif == nf - 1 # True if it's the final subband
            f_start = fres * float(iif) + self.f_min # frequency at the bottom of the subband
            copy_subband = False
            if not is_top_subband: # if it's one of the bottom channels
                f_end = f_start + fres
                f_middle = f_start + fres/2.0 - correction # Middle freq of subband less 0.5x resolution
                delta_t_local = self._calc_delta_t(f_start, f_end) 
            else: # if this is the top output subband
                if do_copy:
                    f_end = f_start + self._df_top*2.0
                    f_middle = f_start + self._df_top - correction # middle freq of subband less 0.5 resolution
                    copy_subband = True
                    delta_t_local = self._ndt_top
                else: # there are 2 subbands available in the input data. The width of the output subband is the sum fo the input width (which is fres/2.0 plus the subband width)
                    f_end = f_start + self._df_top # frequency of the top of the subband
                    f_middle = f_start + fres/2.0 - correction
                    delta_t_local = self._calc_delta_t(f_start, f_end) 
                    self._ndt_top = delta_t_local
                    copy_subband = False

            #f_middle_larger = (f_end - f_start)/2.0 + f_start + correction # Frequency of the middle - with a bit extra - for rounding calculation
            f_middle_larger = f_middle + 2*correction # Middle freq of subband + 0.5x resolution

            #print 'Iterno={} iif={} f_start={} f_end={} f_middle={} delta_t_local={} copy?={}'.format(intnum, iif, f_start, f_end, f_middle, delta_t_local, copy_subband)

            # Fix detla_t_local if we've made a mistake. Gross but it happens for some parameters
            if delta_t_local > self.max_dt:
                delta_t_local = self.max_dt

            # save per-subband info for posterity
            idt_data = []
            nf_data.append((f_start, f_end, f_middle, f_middle_larger, delta_t_local, idt_data))
            
            # for each DM in this subband
            for idt in xrange(delta_t_local):
                dt_middle = int(round(idt * cff(f_middle, f_start, f_end, f_start))) # DM of the middle
                dt_middle_index = dt_middle + shift_input # same as dt_middle
                dt_middle_larger = int(round(idt*cff(f_middle_larger, f_start, f_end, f_start))) # dt at slightly larger freq
                dt_rest = idt - dt_middle_larger # remaining dt
                dt_rest_index = dt_rest + shift_input # same as dt_rest

                # The sum_* values are the whole point of all this stuff. They're 3 tuples containing
                # (subband, dm, time offset) for the 
                sum_dst_start = (iif, idt+shift_output, dt_middle_larger)  # output
                sum_src1_start = (2*iif, dt_middle_index, dt_middle_larger) # lower channel of input
                sum_src2_start = (2*iif + 1, dt_rest_index, 0) # upper channel of input

                # I like this nomenclature better - keeping the previous nomenclature to make easy reference to Barak's code
                id1 = dt_middle_index
                id2 = dt_rest_index
                offset = dt_middle_larger

                # we'll use id2=-1 as a flag to say copy
                if is_top_subband and do_copy:
                    id1 = idt
                    id2 = -1
                    offset = 0
                    

                # dt_middle_larger is also known as mint
                idt_data.append((dt_middle, id1, offset, id2, sum_dst_start, sum_src1_start, sum_src2_start))

    def initialise(self, din):
        '''
        Returns the initial state of the FDMT given the supplied input.

        The state[c, idt, t] is equal to state[c,idt-1, t] + the input

        :TODO: Save last idt of the current block for future block and apply
        :TODO: This will also need to be a kernel that runs on the GPU.

        :din: input array must have shape (nf, nt)
        :returns: array with shape (nf, nd, nt) and same dtype as input.
        '''
        assert din.shape == (self.n_f, self.n_t)
        outshape = (self.n_f, self.init_delta_t, self.n_t+self.init_delta_t)
        state = np.zeros(outshape, dtype=din.dtype)
        idt = 0
        state[:, 0, 0:self.n_t] = din
        for idt in xrange(1, self.init_delta_t):
            state[:, idt, idt:idt+self.n_t] = state[:, idt-1, idt:idt+self.n_t] + din


        return state
        

    def _execute_iteration(self, iterno, din):
        '''
        Executes a single iteration of the FDMT with the supplied input state
        
        :iterno: iteration number to apply  > 0
        :din: input state - shape should be equal to self.hist_state_shape[iterno]
        :returns: array with output state shape self.hist_state_shape[iterno+1] with same dtype as input
        '''
        assert iterno >=0
        nfd = self.hist_nf_data[iterno]
        assert np.all(din.shape == self.hist_state_shape[iterno]), 'Invalid input shape. Was: {} expected {}'.format(din.shape, self.hist_state_shape[iterno])
        out_shape = self.hist_state_shape[iterno+1]
        dout = np.ones(out_shape, dtype=din.dtype)*np.nan # Set to NaN so we can check aftwards that we've filled everything in correctly - but this can be empty
        nchan, ndt, nt_out = out_shape
        # Not all of the out_shape is filled - ndt is filled at the lower end of the band, but we don't do all of it at the higher
        # end - it wastes operations.
        assert len(nfd) == nchan, 'Expected {} channels from iteration. Got {} chanels'.format(len(nfd), nchan)

        # Size of input vector in samples
        nt_in = din.shape[2]

        for ichan in xrange(nchan):
            chanconfig = self.hist_nf_data[iterno][ichan][-1]
            assert ndt >= len(chanconfig)
            # Only do idt up to len(chanconfig) - even though this is less than ndt - this saves on flops
            
            for idt in xrange(len(chanconfig)):
                config = chanconfig[idt]
                # TODO: Make this config a little easier to grok than just a tuple
                _, id1, offset, id2, _, _, _ = config

                # TODO: for those interested in caching, id1 and id2 are

                #print 'ichan={} idt={} id1={} id2={} offset={} din.shape={} dout.shape={}'.format(ichan,  idt, id1, id2, offset, din.shape, dout.shape)
                in1 = din[2*ichan, id1, :] # first channel 
                out = dout[ichan, idt, :]
                if id2 == -1: # copy - all idt in this channel. id2 should be -1 for all idt for this channel
                    # This should be the highest input subband
                    assert ichan == nchan - 1, 'Help, invalid subband'
                    copy_kernel(in1, out)
                else:
                    in2 = din[2*ichan+1, id2, :] # second channel is the next channel up from the first channel
                    add_offset_kernel2(in1, in2, out, offset)

        #assert not np.any(np.isnan(dout)), 'Help! Some data not written'

        return dout

    def execute(self, din):
        '''
        Executes this FDMT on the supplied input
        
        :din: input data. Size=(nf, nt)
        :returns: output data. Size=(nd, nt)
        '''
        state = self.initialise(din)
        niter = len(self.hist_nf_data)

        # This naieve loop does a new malloc every iteration. Ideally you'd malloc 2 buffers and pingpong
        # between them, but numpy doesn't really like doing that, as the size of the state changes between iterations.
        for i in xrange(niter):
            state = self._execute_iteration(i, state)

        # final state has a single 'channel' - remove this axis
        return state[0, :, :]

    @property
    def max_state_size(self):
        '''
        Return the largest state in elements required to run this FDMT
        '''
        return max([s.prod() for s in self.hist_state_shape])
    
    def __call__(self, din):
        return self.execute(din)


class OverlapAndSum(object):
    '''
    Implements an overlap and sum operation so you can get full S/N
    on FRBs with DMS that are larger than the block size

    It keeps an (nd, nd+nt) sized history buffer which it uses to maintain the state. About half of the buffer is unused. Future versions could fix this.

    
    '''
    def __init__(self, nd, nt, dtype=None):
        '''
        Creates a new overlap and sum buffer - the history size is
        (nd, nd+nt). The output is stored in the first nt samples
        
        :nd: number of dispersion trials must be >= nt
        :nt: block size in samples > 0
        '''
        self.nd = nd
        self.nt = nt
        assert nd > 0
        assert nt > 0
        assert nd >= nt
        self.history = np.zeros((nd, nd + nt), dtype)
        
    def process(self, block):
        '''
        Processes the given block of data and returns the most recent block
        
        :block: input data - shape=(nd, nt)
        '''
        
        nd, nt = self.nd, self.nt
        assert block.shape == (nd, nd+nt), 'Invalid block shape {}'.format(block.shape)
        
        # Update history - left most (lowest time index values) get updated to the previous history
        # shifted by nt plus the input block
        self.history[:, 0:nd] = self.history[:, nt:nd+nt] + block[:, 0:nd]
        
        # for times > nd, we just copy the input block in - we have nothing to add to it
        # We'll explicit with the slice boundaries here, for clarity
        self.history[:, nd:nd+nt] = block[:, nd:nd+nt]
        
        # THe output block is the first nt samples of the history
        output = self.history[:, 0:nt]
        
        return output
    

    def __call__(self, block):
        return self.process(block)
        


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
