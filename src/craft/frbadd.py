#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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


    nchan = 336 # channels
    fch1 = 1.448 # GHz
    foff =  -1/1e3 # GHz
    tsamp = 1.2 # milliseconds
    
    while True:
        t = np.random.rand()*5*0 + 5 # random seconds between 5 and 10
        dm = np.random.rand()*400
        dm = 3000
        widthms = np.random.rand()*5 + 0.064
        widthsamp = widthms/tsamp
        widthsamp = 0.1
        
        #d = file.readthismeanyseconds(t)
        nsamp = int(t*1000/tsamp)
        d = (np.random.randn(nchan, nsamp)*18 + 128).astype(np.uint8)
        print(t)

        frb = np.zeros((nchan, nsamp))
        times = np.arange(nsamp)
        for c in range(nchan):
            freq = fch1 + c*foff
            dispersion_delay_ms = 4.15*dm*(fch1**-2 - freq**-2)
            dispersion_delay_samp = abs(dispersion_delay_ms/tsamp) + 50 # offset by a bit so the first few samples aren't off the end of the thing

            smearing_width_ms = 4.15*dm*((freq-0.5*foff)**-2 - (freq+0.5*foff)**-2)
            smearing_width_samp = abs(smearing_width_ms/tsamp)

            total_width2 = smearing_width_samp**2 + widthsamp**2
            total_width = np.sqrt(total_width2)
            amplitude = 1.

            print(amplitude, dm, widthms, widthsamp, freq, dispersion_delay_ms, dispersion_delay_samp, total_width2)
            x = amplitude*np.exp(-(times-dispersion_delay_samp)**2/total_width2/2.)
            tstart = int(dispersion_delay_samp - total_width)
            tend = int(dispersion_delay_samp + total_width)
            frb[c, :] = x
            frb[c, tstart+200:tend+200+1] += amplitude
            #$pylab.plot(times,x )
            #pylab.show()

        pylab.figure()
        pylab.imshow(frb, aspect='auto', interpolation='nearest')
        pylab.figure()
        pylab.plot(frb.T)
        


        pylab.show()
        
    
    

if __name__ == '__main__':
    _main()
