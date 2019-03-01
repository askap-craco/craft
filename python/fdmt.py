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
        self.f_min = float(f_min)
        self.f_max = float(f_max)
        assert(self.f_min < self.f_max)
        self.bw = self.f_max - self.f_min
        self.n_f = int(n_f)
        self.d_f = self.bw / float(self.n_f)
        self.niter = int(np.ceil(np.log2(n_f)))
        self.max_dt = int(max_dt)
        self.n_t = int(n_t)
        freqs = np.arange(nf)*self.df + self.f_min

        self.init_delta_t = self.calc_delta_t(self.f_min, self.f_min + self.d_f)
        self._state_shape = np.array([self.n_f, self.init_delta_t, self.max_dt])
        self.hist_delta_t = []
        self.hist_state_shape = [self._state_shape]
        self.hist_nf_data = []

        # frequencies of subbands, indexed by iteration - useful for tracking non log2 setups, or stuff with gaps (one day).
        self.hist_freqs = [freqs] 
        
        for i in xrange(1, self.niter+1):
            self.iteration(i, n_t)


    def calc_delta_t(self, f_start, f_end):
        rf = cff(f_start, f_end, self.f_min, self.f_max)
        delta_tf = (float(self.max_dt) - 1.0)*rf
        delta_t = int(np.ceil(delta_tf))
        return delta_t
                                                 
    def iteration(self, intnum, n_t):
        delta_f = 2**(intnum)*self.d_f # channel width in MHz
        delta_t = self.calc_delta_t(self.f_min, self.f_min + delta_f) # Max IDT for this iteration
        frange = self.bw # bandwidth
        s = self.hist_state_shape[intnum-1] # input state shape
        nf = s[0]//2 + s[0]% 2 # output number of channels - accounts for non power of 2.
        fjumps = float(nf) # output number of subbands
        state_shape = np.array([nf, delta_t + 1, n_t])
        
        correction = 0.0
        if intnum > 0: # this is never invoked - it's a leftover from Barak's code
            correction = self.d_f/2.0

        # shift input and shift output are never used - they're leftovers from baraks code
        shift_input = 0
        shift_output = 0

        # keep delta_t, state_shape, and per_subband parameters for posterity
        self.hist_delta_t.append(delta_t)
        self.hist_state_shape.append(state_shape)
        nf_data = []
        self.hist_nf_data.append(nf_data)

        # for each ouput subband
        for iif in xrange(nf):
            f_start = frange/fjumps * float(iif) + self.f_min # frequency at the bottom of the subband
            f_end = frange/fjumps*float(iif + 1) + self.f_min # frequency of the top of the subband
            f_middle = (f_end - f_start)/2.0 + f_start - correction # Frequency of the middle of th subband
            f_middle_larger = (f_end - f_start)/2.0 + f_start + correction # Frequency of the middle - with a bit extra - for rounding calculation
            # max DM that we'll comute for this subband
            delta_t_local = self.calc_delta_t(f_start, f_end) + 1
            
            idt_data = []

            # save per-subband info for posterity
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
