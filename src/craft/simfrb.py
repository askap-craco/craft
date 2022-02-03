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
import warnings

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def dmchanof(f1, foff, tint, dm):
    if dm == 0:
        currfreq = f1
    else:
        currfreq = ((f1)**-2 + tint/4.15/dm)**(-0.5)
        
    currchan = int(np.round(-(f1 - currfreq)/foff))
    #print 'dm', dm, 't',t, 'tint', tint, 'currfreq', currfreq,'currchan', currchan

    return currchan, currfreq


def dmdelay(dm, f1, f2):
    '''
    Calculates DM delay in milliseconds 
    
    Assumes DM constant of 4.15
    
    :dm: DM (milliseconds)
    :f1: Frequency (GHz)
    :f2: Frequency (GHz)
    :returns: DM delay (Millsecondss)
    '''

    return 4.15*dm*(f1**-2 - f2**-2)


def calc_delta_dm(fch1, foff, nchans, tsamp):
    '''
    Calculates the DM resolution in PC/CM3 given the system parameters
    
    Assumes DM constant of 4.15

    :fch1: First channel center frequency (GHz)
    :foff: Channel offset (GHz)
    :nchans: Number of channels
    :tsamp: samplign time (ms)'
    :returns: Delta DM in PC/CM3

    >>> calc_delta_dm(0.736, 0.001, 256, 1.7)
    0.4948473447941156
    '''
    fend = fch1 + foff*float(nchans - 1)
    ddm = np.abs(tsamp / dmdelay(1.0, fch1, fend))
    return ddm

def dm2idm(fch1, foff, nchans, tsamp, dm):
    '''
    Calculate the DM index in samples of a given DM in PC/CM3

    Assumes DM constant of 4.15

    :fch1: First channel center frequency (GHz)
    :foff: Channel offset (GHz)
    :nchans: Number of channels
    :tsamp: samplign time (ms)
    :dm: Dipserison measure in PC/CM3
    :returns: Delta DM index in samples as a float. you mayneed to round/cast if you want an in

    >>> int(np.round(dm2idm(0.736, 0.001, 256, 1.7, 100)))
    202
    '''
    ddm = calc_delta_dm(fch1, foff, nchans, tsamp)
    idm = dm/ddm
    return idm

def idm2dm(fch1, foff, nchans, tsamp, idm):
    '''
    Calculates DM in pc/cm3 given a dm index in samples

    Assumes DM constant of 4.15

    :fch1: First channel center frequency (GHz)
    :foff: Channel offset (GHz)
    :nchans: Number of channels
    :tsamp: samplign time (ms)
    :idm: DM index in samples

    :returns: DM in pc/cm3

    >>> int(np.round(idm2dm(0.736, 0.001, 256, 1.7, 202)))
    100
    
    '''
    ddm = calc_delta_dm(fch1, foff, nchans, tsamp)
    dm = idm*ddm
    return dm



def mkfrb(f1, foff, nchans, tsamp, dm, amp=1, offset=0, noiserms=0, ntimes=4096, dclevel=0):
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
    :dclevel: DC level
    :returns: Numpy array with shape (ntimes, nchans)
    '''
    assert f1 > 0
    assert f1 < 10, 'You probably have the wrong units for frequency. f1 is in GHz. f1={}'.format(f1)
    assert abs(foff) < 1, 'You probably have the wrong units for foff. It is in GHz. foff={}'.format(foff)
    assert tsamp > 0
    if tsamp > 100:
        warnings.warn('You probably have the wrong units for tsamp. Its in milliseconds. tsamp={}'.format(tsamp))
        
    assert amp > 0
    assert offset >= 0

    f2 = f1 + float(nchans) * foff
    fstart = f1 + foff*0.5
    shape = (ntimes, nchans)
    d = np.random.randn(np.prod(shape))*noiserms + dclevel
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
        #assert tin1 > 0, 'Invalid tin1={}'.format(tin1)
        #assert tin2 > 0

        if c2 < c1: # Swap
            (c1, c2) = (c2, c1)
            

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

def mkfrb2(f1, foff, nchans, tsamp, dm, amp=1, toffset=0, noiserms=0, ntimes=4096, dclevel=0):
    '''
    Make a simple time-frequency waterfall plot containing a width=1 FRB.
    (mkfrb logic sux - trying again)

    :f1: frequency of the first channel in GHz
    :foff: Ofset between channels in GHz
    :nchans: Number of channels
    :tsamp: sampleing time, in milliseconds
    :dm: dispersion measure - in PC/CM3
    :amp: Amplitude
    :offset: Time offset in ms
    :noiserms: RMS of noise to add
    :ntimes: block size
    :dclevel: DC level
    :returns: Numpy array with shape (ntimes, nchans)
    '''
    assert f1 > 0
    assert f1 < 10, 'You probably have the wrong units for frequency. f1 is in GHz. f1={}'.format(f1)
    assert tsamp > 0
    if tsamp > 100:
        warnings.warn('You probably have the wrong units for tsamp. Its in milliseconds. tsamp={}'.format(tsamp))

    ntimes = int(ntimes)
    nchans = int(nchans)

    freqs = f1 + np.arange(nchans)*foff
    ftop = freqs.max()
    shape = (ntimes, nchans)
    d = np.random.randn(np.prod(shape))*noiserms + dclevel
    d.shape = shape

    for t in range(ntimes):
        tstart_ms = t*tsamp # Beginnignof this integration
        tend_ms = (t+1)*tsamp # end of this integration
        for ic, f in enumerate(freqs):
            ttop_ms = dmdelay(dm, f+0.5*foff, ftop) + toffset # time FRB enters the top of this channel
            tbot_ms = dmdelay(dm, f-0.5*foff, ftop) + toffset # time FRB exits the bottom of the channel
            tsmear = abs(ttop_ms - tbot_ms) # Time from top to bottom

            # Time FRB arrived in this time/frequency cell
            tentry = max(ttop_ms, tstart_ms)

            # Time FRB exited this time/frequency cell
            texit = min(tbot_ms, tend_ms)
            
            # time spent inside this time-frequency cell
            tinside = texit - tentry

            # fraction of the smearing time spent in this channel
            if tsmear == 0:
                tfrac = 1
            else:
                tfrac = tinside/tsmear

            if tinside > 0 or (tinside == 0 and tstart_ms <= tentry < tend_ms):
                #print t, ic, f, tstart_ms, tend_ms, ttop_ms, tbot_ms, tsmear, tinside, tfrac, tentry, texit
                assert tfrac > 0
                d[t,ic] += tfrac*amp

    return d


def mkfrb_fdmt(f1, foff, nchans, tsamp, dm, amp=1, toffset=0, noiserms=0, ntimes=4096, dclevel=0):
    '''

    '''


    assert f1 > 0
    assert f1 < 10, 'You probably have the wrong units for frequency. f1 is in GHz. f1={}'.format(f1)
    assert tsamp > 0
    if tsamp > 100:
        warnings.warn('You probably have the wrong units for tsamp. Its in milliseconds. tsamp={}'.format(tsamp))


    idm = int(np.round(dm2idm(f1, foff, nchans, tsamp, dm)))
    ntimes = int(ntimes)
    nchans = int(nchans)
    assert idm < ntimes
    freqs = f1 + np.arange(nchans)*foff
    shape = (ntimes, nchans)
    d = np.random.randn(np.prod(shape))*noiserms + dclevel
    d.shape = shape
    from . import fdmt
    thefdmt = fdmt.Fdmt(f1, foff, nchans, max_dt=idm+1, n_t=ntimes)
    toffset_samp = int(np.round(float(toffset)/float(tsamp)))
    d = thefdmt.add_frb_track(idm, d.T, amp, toffset_samp)
    return d.T


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
