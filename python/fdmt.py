#!/usr/bin/env python
"""
FDMT class - borrowed from Barak Zackay

Copyright (C) CSIRO 2017
"""
import numpy as np
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def isquare(f):
    return 1./(f*f)

def cff(f1_start, f1_end, f2_start, f2_end):
    return (isquare(f1_start) - isquare(f1_end))/(isquare(f2_start) - isquare(f2_end))

class Fdmt(object):
    def __init__(self, f_min, f_max, n_f, max_dt, n_t):
        '''
        :f_min: minimum frequency in MHz
        :f_max: Maximum frequency in MHz
        :n_f: Numberof channels
        :max_dt: Number of DM trials
        :n_t: Number of samples in input block
        '''
        self.f_min = float(f_min)
        self.f_max = float(f_max)
        assert(self.f_min < self.f_max)
        self.bw = self.f_max - self.f_min
        self.n_f = int(n_f)
        self.d_f = self.bw / float(self.n_f)
        self.niter = int(np.ceil(np.log2(n_f)))
        self.max_dt = int(max_dt)
        self.n_t = int(n_t)
        freqs = np.arange(n_f)*self.d_f + self.f_min
        self.init_delta_t = self.calc_delta_t(self.f_min, self.f_min + self.d_f)
        self._state_shape = np.array([self.n_f, self.init_delta_t, self.max_dt])
        self.hist_delta_t = []
        self.hist_state_shape = [self._state_shape]
        self.hist_nf_data = []

        # frequencies of subbands, indexed by iteration - useful for tracking non log2 setups, or stuff with gaps (one day).
        self.hist_cfreqs = [freqs]
        self.hist_bws = np.ones(len(freqs))*self.d_f

        # channel width of top (copied) and bottom (i.e. all remaining) channels
        self._df_top = self.d_f
        self._df_bot = self.d_f
        self._ndt_top = self.init_delta_t
        
        for i in xrange(1, self.niter+1):
            self.save_iteration(i)

    def calc_delta_t(self, f_start, f_end):
        rf = cff(f_start, f_end, self.f_min, self.f_max)
        delta_tf = (float(self.max_dt) - 1.0)*rf
        delta_t = int(np.ceil(delta_tf))
        return delta_t
                                                 
    def save_iteration(self, intnum):
        n_t = self.n_t
        df = self.d_f
        s = self.hist_state_shape[intnum-1] # input state shape
        nf_in = s[0]
        nf = nf_in//2 + nf_in % 2 # output number of channels - Includes copied channel, if required
        do_copy = nf_in % 2 == 1 # True if we have an odd number of input channels and the top one will be copied
        fjumps = float(nf) # output number of subbands
        state_shape = np.array([nf, delta_t + 1, n_t])
        if do_copy:
            pass # top channel width unchanged
        else:
            self._df_bot += self._df_bot # Top channel will be wider by the new channel

        self._df_bot *= 2.0 # Bottom channels will be added together

        if nf == 1: # if this is the last iteration
            delta_f = self._df_top
        else:
            delta_f = self._df_bot

        fres = self._df_bot
        # delta_f = 2**(intnum)*self.d_f # channel width in MHz - of the normal channels
        delta_t = self.calc_delta_t(self.f_min, self.f_min + delta_f) # Max IDT for this iteration
        
        correction = 0.0
        if intnum > 0: # this is always invoked - it's a leftover from Barak's code
            correction = self.d_f/2.0

        # shift input and shift output are never used - they're leftovers from barak's code
        shift_input = 0
        shift_output = 0

        # keep delta_t, state_shape, and per_subband parameters for posterity
        self.hist_delta_t.append(delta_t)
        self.hist_state_shape.append(state_shape)
        nf_data = []
        self.hist_nf_data.append(nf_data)

        # for each ouput subband
        for iif in xrange(nf):
            is_top_subband = iif == nf - 1 # True if it's the final subband
            f_start = frange/fjumps * float(iif) + self.f_min # frequency at the bottom of the subband
            if not is_top_subband: # if it's one of the bottom channels
                f_end = f_start + fres
                f_middle = f_start + fres/2.0 - correction # Middle freq of subband less 0.5x resolution
                delta_t_local = self.calc_delta_t(f_start, f_end) + 1
            else: # if this is the top output subband
                if do_copy:
                    f_end = f_start + self._df_top*2.0
                    f_middle = f_start + self._df_top - correction # middle freq of subband less 0.5 resolution
                    copy_subband = True
                    delta_t_local = self._ndt_top
                else: # there are 2 subbands available in the input data. The width of the output subband is the sum fo the input width (which is fres/2.0 plus the subband width)
                    f_end = self.f_min + self._df_top # frequency of the top of the subband
                    f_middle = f_start + fres/2.0 - correciton
                    delta_t_local = self.calc_delta_t(f_start, f_end) + 1
                    self._ndt_top = delta_t_local

            #f_middle_larger = (f_end - f_start)/2.0 + f_start + correction # Frequency of the middle - with a bit extra - for rounding calculation
            f_middle_larger = f_middle + 2*correction # Middle freq of subband + 0.5x resolution

            # Fix detla_t_local if we've made a mistake. Gross but it happens for some parameters
            if delta_t_local > self.ndt:
                delta_t_local = self.ndt
            

            # save per-subband info for posterity
            idt_data = []
            nf_data.append((f_start, f_end, f_middle, f_middle_larger, delta_t_local, idt_data))
            
            # for each DM in this subband
            for idt in xrange(delta_t_local):
                dt_middle = np.round(idt * cff(f_middle, f_start, f_end, f_start)) # DM of the middle
                dt_middle_index = dt_middle + shift_input # same as dt_middle
                dt_middle_larger = np.round(idt*cff(f_middle_larger, f_start, f_end, f_start)) # dt at slightly larger freq
                dt_rest = idt - dt_middle_larger # remaining dt
                dt_rest_index = dt_rest + shift_input # same as dt_rest

                # The sum_* values are the whole point of all this stuff. They're 3 tuples containing
                # (subband, dm, time offset) for the 
                sum_dst_start = (iif, idt+shift_output, dt_middle_larger)  # output
                sum_src1_start = (2*iif, dt_middle_index, dt_middle_larger) # lower channel of input
                sum_src2_start = (2*iif + 1, dt_rest_index, 0) # upper channel of input

                idt_data.append((dt_middle, dt_middle_index, dt_middle_larger, dt_rest_index, sum_dst_start, sum_src1_start, sum_src2_start))


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
