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
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def dmchanof(f1, foff, tint, dm):
    currfreq = ((f1)**-2 + tint/4.15/dm)**(-0.5)
    currchan = int(np.round(-(f1 - currfreq)/foff))
    #print 'dm', dm, 't',t, 'tint', tint, 'currfreq', currfreq,'currchan', currchan

    return currchan, currfreq


def dmdelay(dm, f1, f2):
    return 4.15*dm*(f1**-2 - f2**-2)

def tinside(f1, foff, c, dm):
    fmid = f1  + c * foff
    ftop = fmid + foff*0.5
    fbot = fmid - foff*0.5
    t = abs(dmdelay(dm, fbot, ftop))
    return t


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', type=float, default=128.5)
    parser.add_argument('-n','--noiserms', type=float, default=0.)
    parser.add_argument('-d','--dm', type=float, help='DM of a single pulse')
    parser.add_argument('--idt', type=int, help='DM of single pulse in FDMT dt units')
    parser.add_argument('-p','--pulse-amplitude', type=float, help='Amplitude of single pulse', default=1.)
    parser.add_argument('-s','--show', action='store_true')
    parser.add_argument('-c','--nchan', type=int, help='Number of channels', default=336)
    parser.add_argument('--toffset', type=int, help='Sampel offset to start dpulse', default=1500)
    parser.add_argument('--ntimes', type=int,help='Duration of filterbank in samples', default=4096)
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
    
    assert 0 <= values.offset <= 255
    if values.idt is not None:
        t_total = float(values.idt) * tsamp
        values.dm = -t_total/(f1**-2 - f2**-2)/4.15
        assert values.dm >= 0
        print('idt=', values.idt, 'dm=', values.dm)

    sfile = sigproc.SigprocFile(values.files[0], 'w', hdr)
    sfile.seek_data(0)

    if values.dm is None:
        for dm in range(1, 2000):
            print(dm, hdr, values)
            d = mkfrb(float(dm), hdr, values)
            print('Writing ', d.shape, d.dtype, 'to', values.files[0])
            d.tofile(sfile.fin)

    else:
        d = mkfrb(values.dm, hdr, values)
        print('Writing ', d.shape, d.dtype, 'to', values.files[0])
        d.tofile(sfile.fin)
        
    sfile.fin.flush()
    sfile.fin.close()



def mkfrb(dm,  hdr, values):
    nchans = hdr['nchans']
    nifs = hdr['nifs']
    ntimes = values.ntimes
    shape = (ntimes, nifs, nchans)
    f1 = hdr['fch1']/1e3
    foff = hdr['foff']/1e3
    f2 = hdr['fch1']/1e3 + float(nchans) * foff
    tsamp = hdr['tsamp']*1000.
    fstart = f1 + foff*0.5

    offset = values.toffset

    d = np.random.randn(np.prod(shape))*values.noiserms + values.offset
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

        print('dm', dm, 't', t, c1, c2, t1, t2, frqt1, frqt2)
        if c1 < 0 or c1 >= nchans or (t + offset) >= ntimes:
            print('breaking', c1, nchans, (t+offset), ntimes)
            break

        #d[t+offset, 0, c1:c2+1] += values.pulse_amplitude
        tin1 = dmdelay(dm, frqt1, fc1-foff/2.) # time spent between the pulse start freq, and bottom of the channel
        tin11 = dmdelay(dm, fc1+foff/2., fc1-foff/2.) # Time we could have spent, had the pulse started in the top, rather than the side

        tin2 = dmdelay(dm, fc2+foff/2., frqt2)
        tin22 = dmdelay(dm, fc2+foff/2., fc2-foff/2.)

        print('c1',c1,'f1',fc1,'tin1',tin1, 'tin11', tin11, tin1/tin11)
        print('c2',c2,'f2',fc2,'tin2',tin2, 'tin22', tin22, tin2/tin22)
        assert tin1 > 0
        assert tin2 > 0

        assert c2 >= c1
        amp = values.pulse_amplitude
        if c1 == c2:
            d[t+offset, 0, c1] += amp
        else:
            d[t+offset, 0, c1] += amp*(tin1/tin11)**1
            d[t+offset, 0, c2] += amp*(tin2/tin22)**1
            print('Setting %d=%f and %d=%f'% (c1, amp*tin1/tin11, c2, amp*tin2/tin22), d[t+offset, 0, c1], d[t+offset, 0, c2])

        if c2 - c1 > 1:
            d[t + offset, 0, c1+1:c2] += amp
            print('Setting', c1+1, c2-1, 'to', amp, d[t+offset, 0, c1+1])

                

            #d[t+offset, 0, c1:c2+1] += values.pulse_amplitude
        t += 1
            
    if values.show:
        pylab.imshow(d)
        pylab.show()
                            


    d = d.astype(np.uint8)
    return d
    

if __name__ == '__main__':
    _main()
