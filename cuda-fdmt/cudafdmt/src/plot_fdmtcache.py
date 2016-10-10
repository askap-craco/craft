#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def isquare(f):
    return 1./(f*f)

def cff(f1_start, f1_end, f2_start, f2_end):
    return (isquare(f1_start) - isquare(f1_end))/(isquare(f2_start) - isquare(f2_end))

class Fdmt(object):
    def __init__(self, f_min, f_max, n_f, max_dt, n_t):
        assert(f_min < f_max)
        self.f_min = float(f_min)
        self.f_max = float(f_max)
        self.bw = self.f_max - self.f_min
        self.n_f = int(n_f)
        self.d_f = self.bw / float(self.n_f)
        self.niter = int(np.ceil(np.log2(n_f)))
        self.max_dt = int(max_dt)
        self.n_t = int(n_t)

        self.init_delta_t = self.calc_delta_t(self.f_min, self.f_min + self.d_f)
        self._state_shape = np.array([self.n_f, self.init_delta_t, self.max_dt])
        self.hist_delta_t = []
        self.hist_state_shape = [self._state_shape]
        self.hist_nf_data = []
        
        
        for i in xrange(1, self.niter+1):
            self.iteration(i, n_t)


    def calc_delta_t(self, f_start, f_end):
        rf = cff(f_start, f_end, self.f_min, self.f_max)
        delta_tf = (float(self.max_dt) - 1.0)*rf
        delta_t = int(np.ceil(delta_tf))
        return delta_t
                                                 
    def iteration(self, intnum, n_t):
        delta_f = 2**(intnum)*self.d_f
        delta_t = self.calc_delta_t(self.f_min, self.f_min + delta_f)
        frange = self.bw
        s = self.hist_state_shape[intnum-1]
        nf = s[0]/2 + s[0]% 2
        fjumps = float(nf)
        state_shape = np.array([nf, delta_t + 1, n_t])
        
        correction = 0.0
        if intnum > 0:
            correction = self.d_f/2.0

        shift_input = 0
        shift_output = 0

        self.hist_delta_t.append(delta_t)
        self.hist_state_shape.append(state_shape)
        nf_data = []
        self.hist_nf_data.append(nf_data)
        
        for iif in xrange(nf):
            f_start = frange/fjumps * float(iif) + self.f_min
            f_end = frange/fjumps*float(iif + 1) + self.f_min
            f_middle = (f_end - f_start)/2.0 + f_start - correction
            f_middle_larger = (f_end - f_start)/2.0 + f_start + correction
            delta_t_local = self.calc_delta_t(f_start, f_end) + 1
            
            idt_data = []
            nf_data.append((f_start, f_end, f_middle, f_middle_larger, delta_t_local, idt_data))
            

            for idt in xrange(delta_t_local):
                dt_middle = np.round(idt * cff(f_middle, f_start, f_end, f_start))
                dt_middle_index = dt_middle + shift_input
                dt_middle_larger = np.round(idt*cff(f_middle_larger, f_start, f_end, f_start))
                dt_rest = idt - dt_middle_larger
                dt_rest_index = dt_rest + shift_input

                sum_dst_start = (iif, idt+shift_output, dt_middle_larger)
                sum_src1_start = (2*iif, dt_middle_index, dt_middle_larger)
                sum_src2_start = (2*iif + 1, dt_rest_index, 0)

                idt_data.append((dt_middle, dt_middle_index, dt_middle_larger, dt_rest_index, sum_dst_start, sum_src1_start, sum_src2_start))
                

def plot_cache(fdmt):
    print(fdmt.hist_state_shape)

    pylab.plot(fdmt.hist_delta_t)
    pylab.xlabel('Iteration')
    pylab.ylabel('Delta_t')

    pylab.figure()
    niter = len(fdmt.hist_nf_data)

    n_cache_hits = np.zeros(niter)
    n_mem_hits = np.zeros(niter)

    for iiter, nf_data in enumerate(fdmt.hist_nf_data):
        print 'IITER', iiter
        fig, axes = pylab.subplots(2,2)
        axes = axes.flatten()
        
        for iif, subband_data in enumerate(nf_data):
            idt_data = subband_data[5]
            dt_range = np.arange(len(idt_data))
            dt_middle_larger = [idd[2] for idd in idt_data]
            tshift = dt_middle_larger
            dt1  = np.array([idd[1] for idd in idt_data])
            dt2  = np.array([idd[3] for idd in idt_data])
            axes[0].plot(dt_range, dt1)
            axes[0].set_ylabel('src1 dt')
            axes[1].plot(dt_range, dt2)
            axes[1].set_ylabel('src2 dt')
            axes[2].plot(dt_range, dt1 - dt2)
            axes[2].set_ylabel('src1 dt - src2 dt')
            axes[3].plot(dt_range, tshift)
            axes[3].set_ylabel('Tshift')
            #cache_hit_ratio = float(sum(dt1==dt2) - 1)/float(len(dt1))
            #print iiter, iif, cache_hit_ratio
            n_cache_hits[iiter] += (sum(dt1 == dt2) - 1)
            n_mem_hits[iiter] += len(dt1)
            
            #pylab.plot(dt_range, dt2, label='Subband {}'.format(iif))

        fig.text(0.5, 0.98, 'Iteration {}'.format(iiter), ha='center', va='top')
        fig.text(0.5, 0.02, 'Delta t', ha='center', va='bottom')


        pylab.title('Iteration {}'.format(iiter))

    pylab.figure()
    pylab.plot(n_cache_hits/n_mem_hits)
    pylab.xlabel('Iteration')
    pylab.ylabel('Cache hit ratio')
    print 'Fraction of cache hits', n_cache_hits.sum()/n_mem_hits.sum(), 'Nhits', n_cache_hits.sum(), 'Nmem', n_mem_hits.sum()
    pylab.show()

    

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    nchan = 512
    fdmt = Fdmt(1440.-nchan, 1440., nchan, 512, 512)
    plot_cache(fdmt)
    

if __name__ == '__main__':
    _main()
