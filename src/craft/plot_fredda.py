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
    v = np.fromfile(fin, dtype=dtype, count=theshape.prod())
    v.shape = theshape
    fin.close()
    print(fname, v.shape, len(v), 'Nzeros=', np.sum(v == 0), \
            'max/min/mean/sum {}/{}/{}/{}'.format(v.max(), v.min(), v.mean(), v.sum()), \
            'max at', \
            np.unravel_index(v.argmax(), v.shape), 'NaNs?', np.sum(np.isnan(v)), 'Ninfs?=', np.sum(np.isinf(v)))

    return v

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
    fig = pylab.figure()
    
    for ifname, fname in enumerate(file_series(prefix, start)):
        fig.clear()
        gs = gridspec.GridSpec(1,1)
        p = plt.subplot
        ax = p(gs[0])

        ostate = load4d(fname)
        v = ostate[theslice]
        v = np.ma.masked_invalid(v)
        v = np.ma.masked_equal(v, 0)
        nfreq, ntime = v.shape
        vmid = np.ma.median(v)
        voff = np.std((v - vmid))
        myimshow(ax, (v), aspect='auto', origin='lower')
        ax.set_xlabel('t')
        ax.set_ylabel('dt')
        ax.set_title(fname)

        pylab.show()


        if maxn is not None and ifname >= maxn:
            print('Quitting as maxn exceeded')
            break

def statplot(ax, fname, name):
    d = None
    try :
        d = load4d(fname.replace('fdmt',name))
        sl = [0,0,slice(None),slice(None)]
        ax.plot(d[sl].T)
        ax.set_ylabel(name)
    except IOError:
        print('No data for ', fname, name)
    return d

def show_fdmt_series(prefix, theslice, values, start=0, maxn=10, ibeam=0):
    fig = pylab.figure()
    
    for ifname, fname in enumerate(file_series(prefix, start)):
        ostate = load4d(fname)
        gs = gridspec.GridSpec(2,10)
        p = plt.subplot
        rawax = p(gs[0, 0:2])

        rawname = fname.replace('fdmt','inbuf')
        rawdat = load4d(rawname)
        rawd = rawdat[ibeam, :, 0, :]
        print('Raw inbuf', rawdat.shape)
        nchans, ntimes = rawd.shape
        rawd = np.ma.masked_equal(rawd, 0.0)
        myimshow(rawax, rawd, aspect='auto', origin='lower')

        fdmtax = p(gs[1, 0:2])
        v = ostate[theslice]
        v = np.ma.masked_invalid(v)
        nfreq, ntime = v.shape

        maxpos = np.unravel_index(v.argmax(), v.shape)
        print(fname, ostate.shape, len(v), 'zeros?', np.all(ostate == 0), \
            'max/min/mean/sum {}/{}/{}/{}'.format(v.max(), v.min(), v.mean(), v.sum()), \
            'max at', \
            np.unravel_index(v.argmax(), v.shape), 'NaNs?', np.sum(np.isnan(v)))

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
        specax.set_ylim(min(chans), max(chans))

        dmax2 = p(gs[1,2])
        ndt, _ = v.shape
        dts = np.arange(ndt)
        stdline, = dmax2.plot(v.std(axis=1), dts)
        meanline, = dmax2.plot(v.mean(axis=1), dts)
        dmax2.set_ylabel('Delta_t')
        dmax2.set_xlabel('StdDev/Mean')
        dmax2.set_ylim(min(dts), max(dts))
        dmax2.axvline(0, c=meanline.get_color(), ls=':')
        dmax2.axvline(1, c=stdline.get_color(), ls=':')
        
        dmax = p(gs[:, 3:6])
        maxdm, maxt = maxpos
        if values.dmrange is None:
            dmrange = slice(max(maxdm-10, 0), min(maxdm+10, ndt-1) )
        else:
            dmrange = slice(*values.dmrange)

        dmax.plot(v[dmrange, :].T)
        dmax.set_xlabel('Sample')
        dmax.set_ylabel('S/N')

        meanax = p(gs[0, 6:8])
        stdax = p(gs[0, 8:10])
        kurtax = p(gs[1, 6:8])
        dm0ax = p(gs[1, 8:10])

        statplot(meanax, fname, 'mean')
        statplot(stdax, fname, 'std')
        statplot(kurtax, fname, 'kurt')
        statplot(dm0ax, fname, 'dm0')
        rawdm0 = rawd.sum(axis=1).T
        dm0ax.plot(rawd.sum(axis=0).T)
        try:
            dm0count = load4d(fname.replace('fdmt','dm0count'))
            dm0countax = dm0ax.twinx()
            dm0countax.plot(dm0count[0,0,ibeam,:])
            dm0countax.set_ylim(0, None)
        except:
            print('COunlt plot count')

        pylab.show()

        if maxn is not None and ifname >= maxn - 1:
            print('Quitting as maxn exceeded')
            break

def comma_list(s):
    return list(map(int, s.split(',')))


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--beam', type=int, help='beam number')
    parser.add_argument('-d','--dmrange', type=comma_list, help='Dm range to show in time series plot')
    parser.add_argument('-s','--start', type=int, help='Start block')
    parser.add_argument('-n','--maxn', type=int, help='max number of blocks ot plot')
    parser.set_defaults(verbose=False, beam=0, start=4, maxn=10)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    #show_inbuf_series('ostate_e%d.dat', [0, 0, slice(None),slice(None)], start=values.start, maxn=10)
    #pylab.show()
    show_fdmt_series('fdmt_e%d.dat', [values.beam, 0, slice(None), slice(None)], values, start=values.start, maxn=values.maxn, ibeam=values.beam)

    pylab.show()
    



if __name__ == '__main__':
    _main()
