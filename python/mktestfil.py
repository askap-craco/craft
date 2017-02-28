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
import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', type=float, default=128.)
    parser.add_argument('-n','--noiserms', type=float, default=0.)
    parser.add_argument('-d','--dm', type=float, help='DM of a single pulse')
    parser.add_argument('-p','--pulse-amplitude', type=float, help='Amplitude of single pulse', default=1.)
    parser.add_argument('-s','--show', action='store_true')
    parser.add_argument('-c','--nchan', type=int, help='Number of channels', default=336)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    hdr = {'source_name':'test',
           'tsamp': 0.001265625,
           'pulsarcentric':0,
           'az_start':1.,
           'nbits':8,
           'foff':-1.0,
           'fch1':1448.0,
           'nchans': values.nchan,
           'telescope_id':1,
           'src_dej': -30724.3618819,
           'src_raj': 164015.762746,
           'tstart':57681.5338261,
           'nifs':1
       }
    
    nchans = hdr['nchans']
    nifs = hdr['nifs']
    ntimes = 16384
    shape = (ntimes, nifs, nchans)
    f1 = hdr['fch1']
    foff = hdr['foff']
    f2 = hdr['fch1'] + float(nchans) * foff

    
    d = np.random.randn(np.prod(shape))*values.noiserms + values.offset
    d.shape = shape
    print 'Shape', d.shape

    offset = 617
    if values.dm is not None:
        t = 0
        while True:
            tint = t * hdr['tsamp']*1000.
            currfreq = ((f1/1e3)**-2 + tint/4.15/values.dm)**(-0.5)*1e3
            currchan = int(np.round(-(f1 - currfreq)/foff))
            #print 'dm', values.dm, 't',t, 'tint', tint, 'currfreq', currfreq,'currchan', currchan
            t += 1
            if currchan < 0 or currchan >= nchans or (t + offset) >= ntimes:
                break
            else:
                d[t+offset, 0, currchan] += values.pulse_amplitude

            
    if values.show:
        pylab.imshow(d)
        pylab.show()
                            


    d = d.astype(np.uint8)
    sfile = sigproc.SigprocFile(values.files[0], 'w', hdr)
    sfile.seek_data(0)
    print 'Writing ', d.shape, d.dtype, 'to', values.files[0]
    d.tofile(sfile.fin)
    sfile.fin.flush()
    sfile.fin.close()
    
    

if __name__ == '__main__':
    _main()
