#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2019
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def dmchanof(f1, foff, tint, dm):
    currfreq = ((f1)**-2 + tint/4.15/dm)**(-0.5)
    currchan = int(np.round(-(f1 - currfreq)/foff))
    #print 'dm', dm, 't',t, 'tint', tint, 'currfreq', currfreq,'currchan', currchan

    return currchan, currfreq


def dmdelay(dm, f1, f2):
    return 4.15*dm*(f1**-2 - f2**-2)


def mkfrb(f1, foff, nchans, tsamp, dm, amp=1, offset=0, noiserms=0, ntimes=4096):
    '''
    Make a simple time-frequency waterfall plot containing a width=1 FRB.

    :f1: frequency of the first channel in GHz
    :foff: Ofset between channels in GHz
    :nchans: Number of channels
    :tsamp: sampleing time, in milliseconds
    :dm: dispersion measure - in PC/CM3
    :amp: Amplitude
    :offset: Time offset in samples
    :noiserms: RMS of noise to add
    :ntimes: block size
    :returns: Numpy array with shape (ntimes, nchans)
    '''
    assert f1 > 0
    assert f1 < 10, 'You probably have the wrong units for frequency. f1 is in GHz. f1={}'.format(f1)
    assert abs(foff) < 1, 'You probably have the wrong units for foff. It is in GHz. foff={}'.format(foff)
    assert tsamp > 0
    assert tsamp < 100, 'You probably have the wrong units for tsamp. Its in milliseconds. tsamp={}'.format(tsamp)
    assert amp > 0
    assert offset >= 0

    shape = (ntimes, nchans)
    f2 = f1 + float(nchans) * foff
    fstart = f1 + foff*0.5
    d = np.random.randn(np.prod(shape))*noiserms 
    d.shape = shape
    t = 0    

    while True:
        t1 = t * tsamp
        t2 = (t + 1) * tsamp
        
        c1,frqt1 = dmchanof(f1,foff, t1, dm)
        c2,frqt2 = dmchanof(f1,foff, t2, dm)
        
        fc1 = f1 + c1*foff
        fc2 = f1 + c2*foff
        
        # This logic might not be right
        if c2 >= nchans:
            c2 = nchans-1
            currchan = c1

        #print 'dm', dm, 't', t, c1, c2, t1, t2, frqt1, frqt2
        if c1 < 0 or c1 >= nchans or (t + offset) >= ntimes:
 #           print 'breaking', c1, nchans, (t+offset), ntimes
            break

        #d[t+offset, 0, c1:c2+1] += values.pulse_amplitude
        tin1 = dmdelay(dm, frqt1, fc1-foff/2.) # time spent between the pulse start freq, and bottom of the channel
        tin11 = dmdelay(dm, fc1+foff/2., fc1-foff/2.) # Time we could have spent, had the pulse started in the top, rather than the side

        tin2 = dmdelay(dm, fc2+foff/2., frqt2)
        tin22 = dmdelay(dm, fc2+foff/2., fc2-foff/2.)

#        print 'c1',c1,'f1',fc1,'tin1',tin1, 'tin11', tin11, tin1/tin11
#        print 'c2',c2,'f2',fc2,'tin2',tin2, 'tin22', tin22, tin2/tin22
        assert tin1 > 0
        assert tin2 > 0

        assert c2 >= c1
        if c1 == c2:
            d[t+offset, c1] += amp
        else:
            d[t+offset, c1] += amp*(tin1/tin11)**1
            d[t+offset, c2] += amp*(tin2/tin22)**1
            #print 'Setting %d=%f and %d=%f'% (c1, amp*tin1/tin11, c2, amp*tin2/tin22), d[t+offset, c1], d[t+offset, c2]

        if c2 - c1 > 1:
            d[t + offset, c1+1:c2] += amp
            #print 'Setting', c1+1, c2-1, 'to', amp, d[t+offset, c1+1]

        t += 1

    return d


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

    frb = mkfrb(1.320, -0.001, 256, 1, 100, ntimes=512).T
    pylab.imshow(frb)
    pylab.show()

    
if __name__ == '__main__':
    _main()
