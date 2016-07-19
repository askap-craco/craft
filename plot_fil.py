#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sigproc
from plot_dada import Formatter
import numpy as np
import pylab


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-t', '--times', help='Integration range to plot')
    parser.add_argument('--dtype', help="override dtype", nargs="?")
    parser.add_argument('-n','--normalise', action='store_true', help='Normalise each channel before plotting', default=False)
    parser.add_argument('-p','--period',  help='Folding period', type=float, default=0.089328385024)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = sigproc.SigprocFile(values.files[0])

    if values.dtype is None:
        if f.nbits == 32:
            dtype = np.float32
        elif f.nbits == 8:
            dtype = np.uint8
        else:
            raise ValueError('Cant parse it: %s' % f.nbits)
    else:
        dtype = np.dtype(values.dtype)

    if values.times:
        bits = map(int, values.times.split(','))
        if len(bits) == 1:
            tstart = bits[0]
            ntimes = 128*8
        elif len(bits) == 2:
            tstart, ntimes = bits
        else:
            raise ValueError('Invalid times: %s', values.times)
    else:
        tstart = 0
        ntimes = 128*8

    assert ntimes > 0
    if f.data_type == 1:
        plot_spectra(f, tstart, ntimes, dtype, values)
    elif f.data_type == 2:
        plot_time_series(f, tstart, ntimes, dtype, values)
    else:
        raise ValueError('Unknown data type %s' % f.data_type)


def plot_time_series(f, tstart, ntimes, dtype, values):
    tend = tstart + ntimes
    nelements = ntimes*f.nifs*f.nchans

    f.seek_data(f.bytes_per_element*tstart)
    print nelements, ntimes, f.nchans, f.nifs, dtype
    v = np.fromfile(f.fin, dtype=dtype, count=nelements )
    t = np.arange(len(v))*f.tsamp
    pylab.plot(t, v)
    pylab.title(f.filename)
    pylab.show()


def fold(data, tsamp, period):
    nbins = int(period/tsamp)
    print 'Folding with tsamp {} period {} nbins {} data.shape {}'.format(tsamp, period, nbins, data.shape)
    ntimes, nchans = data.shape
    folded = np.zeros((nbins, nchans), dtype=np.float32)
    
    for i in xrange(ntimes):
        t = i*tsamp;
        phase = t/period
        thebin = i % nbins
        folded[thebin, :] += data[i, :]

    return folded

def plot_spectra(f, tstart, ntimes, dtype, values):
    tend = tstart + ntimes
    nelements = ntimes*f.nifs*f.nchans

    f.seek_data(f.bytes_per_element*tstart)
    v = np.fromfile(f.fin, dtype=dtype, count=nelements )
    print 'Nelements', nelements, 'Ntimes', ntimes, 'nchans', f.nchans, 'nifs', f.nifs, dtype, 'Actual length', len(v)
    
    v.shape = (ntimes, f.nifs, f.nchans)
    v = np.ma.masked_array(v, mask=np.isnan(v))
    nrows = f.nifs/6
    ncols = 6
    nrows = 1
    ncols = 2
    
    fig, axes = pylab.subplots(nrows, ncols, sharex=True, sharey=True)

    fft_fig, fft_axes = pylab.subplots(nrows, ncols, sharex=True, sharey=True)

    hist_fix, hist_axes = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
    
    bandpass_fig, bandpass_axes = pylab.subplots(nrows, ncols, sharex=True, sharey=True)

    if values.period > 0:
        foldfig, foldaxes = pylab.subplots(nrows, ncols, sharex=True, sharey=True)
        if f.nifs == 1:
            foldaxes = [foldaxes]
    else:
        foldfix, foldaxes = None, None


    for ifn in xrange(nrows*ncols):
        ax = axes.flat[ifn]
        hax = hist_axes.flat[ifn]
        bpax = bandpass_axes.flat[ifn]
        fftax = fft_axes.flat[ifn]

        if foldaxes is not None:
            fax = foldaxes[ifn]

        data = v[:,  ifn, :]

        if values.normalise:
            dmean = data.mean(axis=0)
            dstd = data.std(axis=0)

            print 'data', data.shape, dmean.shape, dstd.shape
            data = (data - dmean + 100)/dstd
            #vmin = data.min()
            #vmax = data.max()
            vmin = None
            vmax = None

        else:
            vrange = v.max() - v.min()
            vmid = data.mean()
            vstd = v.std()
            vmax = vmid+2*vstd
            vmin = vmid-2*vstd

        print 'VMIN', vmin, 'VMAX', vmax, data.shape
        im = ax.imshow(data.T, aspect='auto', vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
        ax.text(0, 0, 'ifn %d' % ifn, va='top', ha='left') 
        ax.format_coord = Formatter(im)
        ax.set_xlabel('Integration')
        ax.set_ylabel('Channel')

        #fig.title(f.filename)

        for c in xrange(f.nchans):
            hax.hist(v[:,ifn, c].flatten(), label='IF %d'%ifn, bins=100, histtype='step')

        pylab.title(f.filename)
        pylab.savefig('{}_hist.png'.format(f.filename))

        if values.period > 0:
            vfolded = fold(data, f.tsamp, values.period)
            print 'Plotting folded data'
            fax.imshow(vfolded.T, label='IF %s'%ifn, aspect='auto')
            

        T = ntimes*f.tsamp
        times = np.arange(ntimes)
        periods = 1.0/times


        for chan in xrange(f.nchans):
            fftax.plot(periods, abs(np.fft.fft(v[:, ifn, chan])))

        pylab.title(f.filename)
        pylab.xlabel('Period (s)')
        pylab.savefig('{}_fft.png'.format(f.filename))

        bpax.plot(v[:, ifn, :].mean(axis=0))
        bpax.set_xlabel('Channel')
        bpax.set_ylabel('Mean bandpass')
        bpax2 = bpax.twinx()
        bpax2.plot(v[:, ifn, :].std(axis=0), 'r:')
        bpax2.set_ylabel('Bandpass std')

        pylab.figure()
        dm0 = v[:,ifn,:].mean(axis=1)
        fig, dm0ax = pylab.subplots(2,1)
        dm0ax[0].plot(times, dm0);
        dm0ax[0].set_xlabel('Time (ms)')
        dm0ax[0].set_ylabel('DM=0 time series')
        
        dm0ax[1].semilogy(np.abs(np.fft.rfft(dm0))**2)
        dm0ax[1].set_xlabel('Digital frequency / sample')
        dm0ax[1].set_ylabel('FFT(DM0)**2')

        
        pylab.show()


    
    fig.savefig('{}_dynspec.png'.format(f.filename))

    for i in xrange(f.nifs):
        ifv = v[:, i, 0]
        print 'Channel 0 mean, if=', i, ifv.mean(), ifv.std()
        ifv = v[:, i, 72]
        print 'Channel 72 mean, if=', i, ifv.mean(), ifv.std()


    #pylab.subplots_adjust(wspace=0, hspace=0)
    pylab.show()

    
    
    

if __name__ == '__main__':
    _main()
