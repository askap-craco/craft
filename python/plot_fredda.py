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
import subprocess
from craftsim import dispersed_voltage, dispersed_stft
from FDMT import *
import matplotlib.gridspec as gridspec

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}'.format(x, y, z)

def myimshow(ax, *args, **kwargs):
    kwargs['aspect'] = 'auto'
    kwargs['interpolation'] = 'nearest'
    im = ax.imshow(*args, **kwargs)
    ax.format_coord = Formatter(im)
    return im

def make_signal2():
    freq = np.linspace(fmin, fmax, nf)
    dm = 150
    assert(len(freq) == nf)
    #d = np.ones((nf, nt), dtype=np.float32)
    d = np.zeros((nf, nt), dtype=np.float32)
    d += np.random.randn(nf, nt)
    return d


    N_total = N_f*N_t*N_bins
    PulseLength = N_f*N_bins

def load4d(fname, dtype=np.float32):
    fin = open(fname, 'r')
    theshape =  np.fromfile(fin, dtype=np.int32, count=4)
    d = np.fromfile(fin, dtype=dtype, count=theshape.prod())
    d.shape = theshape
    fin.close()
    return d

def file_series(prefix, start=0):
    i = start
    while True:
        fname = prefix % i
        if os.path.exists(fname):
            yield fname
        else:
            break
        i += 1

def show_inbuf_series(prefix, theslice, start=0, maxn=10):
    for ifname, fname in enumerate(file_series(prefix, start)):
        ostate = load4d(fname)
        fig = pylab.figure()
        v = ostate[theslice]
        print fname, ostate.shape, len(v), 'zeros?', np.all(ostate == 0), 'max', v.max(), np.unravel_index(v.argmax(), v.shape), 'NaNs?', np.sum(np.isnan(v))
        
        v = np.ma.masked_invalid(v)
        nfreq, ntime = v.shape
        vmid = np.ma.median(v)
        voff = np.std((v - vmid))

        

        myimshow(axes[0], (v), aspect='auto', origin='lower')

        #bins = np.linspace(vmid-3,vmid+3,50)
        bins = 50

        for f in xrange(nfreq):
            axes[1].hist(v[f,:], bins=bins, histtype='step')
            #axes[1].loglog(abs(np.fft.rfft(v[120:140, :].sum(axis=0))))

        pylab.title(fname)

        if maxn is not None and ifname >= maxn:
            print 'Quitting as maxn exceeded'
            break

def show_fdmt_series(prefix, theslice, values, start=0, maxn=10, ibeam=0):
    for ifname, fname in enumerate(file_series(prefix, start)):
        ostate = load4d(fname)
        gs = gridspec.GridSpec(2,5)
        p = plt.subplot
        
        rawax = p(gs[0, 0:2])

        rawname = fname.replace('fdmt','inbuf')
        rawdat = load4d(rawname)
        rawd = rawdat[ibeam, :, 0, :]
        print 'Raw inbuf', rawdat.shape
        nchans, ntimes = rawd.shape
        myimshow(rawax, rawd, aspect='auto', origin='lower')

        fdmtax = p(gs[1, 0:2])
        v = ostate[theslice]
        v = np.ma.masked_invalid(v)
        nfreq, ntime = v.shape

        maxpos = np.unravel_index(v.argmax(), v.shape)
        print fname, ostate.shape, len(v), 'zeros?', np.all(ostate == 0), 'max', v.max(), maxpos , 'NaNs?', np.sum(np.isnan(v))
        vmid = np.ma.median(v)
        voff = np.std((v - vmid))
        myimshow(fdmtax, (v), aspect='auto', origin='lower')
        fdmtax.set_title(fname)
        fdmtax.set_ylabel('Delta_t')
        fdmtax.set_xlabel('Sample')
        rawax.set_title(rawname)
        rawax.set_ylabel('Channel')

        specax = p(gs[0,2])
        chans = np.arange(nchans)
        specax.plot(rawd.std(axis=1), chans)
        specax.plot(rawd.mean(axis=1), chans)

        dmax2 = p(gs[1,2])
        ndt, _ = v.shape
        dts = np.arange(ndt)
        dmax2.plot(v.std(axis=1), dts)
        dmax2.plot(v.mean(axis=1), dts)
        dmax2.set_ylabel('Delta_t')
        dmax2.set_xlabel('StdDev/Mean')
        
        dmax = p(gs[:, 3:6])
        maxdm, maxt = maxpos
        if values.dmrange is None:
            dmrange = slice(max(maxdm-10, 0), min(maxdm+10, ndt-1) )
        else:
            dmrange = slice(*values.dmrange)

        dmax.plot(v[dmrange, :].T)
        dmax.set_xlabel('Sample')
        dmax.set_ylabel('S/N')
        

        pylab.show()

        if maxn is not None and ifname >= maxn - 1:
            print 'Quitting as maxn exceeded'
            break

def comma_list(s):
    return map(int, s.split(','))


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--beam', type=float, help='beam number')
    parser.add_argument('-d','--dmrange', type=comma_list, help='Dm range to show in time series plot')
    parser.add_argument('-s','--start', type=int, help='Start block')
    parser.add_argument('-n','--maxn', type=int, help='max number of blocks ot plot')
    parser.set_defaults(verbose=False, beam=0, start=4, maxn=10)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    #show_inbuf_series('inbuf_e%d.dat', [0, slice(None), 0, slice(None)], start=20, maxn=1)
    show_fdmt_series('fdmt_e%d.dat', [values.beam, 0, slice(None), slice(None)], values, start=values.start, maxn=values.maxn, ibeam=values.beam)

    #pylab.show()
    



if __name__ == '__main__':
    _main()
