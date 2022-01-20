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
from . import fdmt as fdmt_class

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

                

def plot_cache(fdmt):
    print((fdmt.hist_state_shape))

    pylab.plot(fdmt.hist_delta_t)
    pylab.xlabel('Iteration')
    pylab.ylabel('Delta_t')

    pylab.figure()
    niter = len(fdmt.hist_nf_data)

    n_cache_hits = np.zeros(niter)
    n_mem_hits = np.zeros(niter)

    for iiter, nf_data in enumerate(fdmt.hist_nf_data):
        print('IITER', iiter)
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

        fig.text(0.5, 0.98, 'Iteration {} nf={} ndt={}'.format(iiter, len(nf_data), len(dt_range)), ha='center', va='top')
        fig.text(0.5, 0.02, 'Delta t', ha='center', va='bottom')


        pylab.title('Iteration {}'.format(iiter))

    pylab.figure()
    pylab.plot(n_cache_hits/n_mem_hits)
    pylab.xlabel('Iteration')
    pylab.ylabel('Cache hit ratio')
    print('Fraction of cache hits', n_cache_hits.sum()/n_mem_hits.sum(), 'Nhits', n_cache_hits.sum(), 'Nmem', n_mem_hits.sum())
    pylab.show()

    

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    nchan = 512
    fdmt = fdmt_class.Fdmt(1440.-nchan, 1440., nchan, 512, 128)
    plot_cache(fdmt)
    

if __name__ == '__main__':
    _main()
