#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
from . import sigproc
from .plot_dada import Formatter
import numpy as np
import pylab
import sys


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-t', '--times', help='Integration range to plot')
    parser.add_argument('--dtype', help="override dtype", nargs="?")
    parser.add_argument('-n','--normalise', action='store_true', help='Normalise each channel before plotting', default=False)
    parser.add_argument('-p','--period',  help='Folding period vela=0.089328385024)', type=float,default=0)
    parser.add_argument('--plot-fft', help='plot fft', action='store_true', default=False)
    parser.add_argument('--plot-image',help='plot image', action='store_true', default=False)
    parser.add_argument('--plot-hist',help='plot hist', action='store_true', default=False)


    parser.add_argument('--nxy', help='rows & columns', default='1,1')
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
        bits = list(map(int, values.times.split(',')))
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
    print(nelements, ntimes, f.nchans, f.nifs, dtype)
    v = np.fromfile(f.fin, dtype=dtype, count=nelements )
    t = np.arange(len(v))*f.tsamp
    pylab.plot(t, v)
    pylab.title(f.filename)
    pylab.show()

def fold(data, tsamp, period):
    nbins = int(period/tsamp)
    print('Folding with tsamp {} period {} nbins {} data.shape {}'.format(tsamp, period, nbins, data.shape))
    ntimes, nchans = data.shape
    folded = np.zeros((nbins, nchans), dtype=np.float32)
    
    for i in range(ntimes):
        t = i*tsamp;
        phase = t/period
        thebin = i % nbins
        folded[thebin, :] += data[i, :]

    return folded

def mysubplots(nrows, ncols, *args, **kwargs):
    fig, axes = pylab.subplots(nrows, ncols, *args, **kwargs)
    if not hasattr(axes, 'flat'):
        axes = [axes]

    return fig, np.array(axes).flatten()
    

def plot_spectra(f, tstart, ntimes, dtype, values):
    tend = tstart + ntimes
    nelements = ntimes*f.nifs*f.nchans

    f.seek_data(f.bytes_per_element*tstart)
    v = np.fromfile(f.fin, dtype=dtype, count=nelements )
    v = v.astype(np.float)
    print('Nelements', nelements, 'Ntimes', ntimes, 'nchans', f.nchans, 'nifs', f.nifs, dtype, 'Actual length', len(v))
    
    v.shape = (ntimes, f.nifs, f.nchans)
#    v.shape = (ntimes, f.nchans, f.nifs)
#    v = np.swapaxes(v, 2, 1)
    v = np.ma.masked_array(v, mask=np.isnan(v))
    #nrows = f.nifs/6

    nrows, ncols = list(map(int, values.nxy.split(',')))

    #plot_cov(v[:, :, 150:151], f)
    
    fig, axes = mysubplots(nrows, ncols, sharex=True, sharey=True)
    
    if values.plot_fft:
        fft_fig, fft_axes = mysubplots(nrows, ncols, sharex=True, sharey=True)

    hist_fix, hist_axes = mysubplots(nrows, ncols, sharex=True, sharey=True)
    
    bandpass_fig, bandpass_axes = mysubplots(nrows, ncols, sharex=True, sharey=True)
    dm0_fig, dm0_ax = mysubplots(1,1)
    dm0ax= dm0_ax[0]

    dm0_fig, dm0_axs = mysubplots(1,1)
    dm0axs = dm0_axs[0]

    if values.period > 0:
        foldfig, foldaxes = mysubplots(nrows, ncols, sharex=True, sharey=True)
    else:
        foldfig, foldaxes = None, None


    ntimes, nbeams, nchans = v.shape

    beam_mean = v.mean(axis=1)

    for ifn in range(nrows*ncols):
        ax = axes.flat[ifn]
        hax = hist_axes.flat[ifn]
        bpax = bandpass_axes.flat[ifn]
        #dm0ax = dm0_ax.flat[ifn]

        if values.plot_fft:
            fftax = fft_axes.flat[ifn]

        if foldaxes is not None:
            fax = foldaxes[ifn]

        data = v[:,  ifn, :]
        dm0 = data.mean(axis=1)
        dm0s = data.std(axis=1)

        if values.normalise:
            dmean = data.mean(axis=0)
            dstd = data.std(axis=0)

            print('data', data.shape, dmean.shape, dstd.shape)
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

        print('VMIN', vmin, 'VMAX', vmax, data.shape)
        if values.plot_image:
            im = ax.imshow(data.T, aspect='auto', vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
            ax.text(0, 0, 'ifn %d' % ifn, va='top', ha='left') 
            ax.format_coord = Formatter(im)
            ax.set_xlabel('Integration')
            ax.set_ylabel('Channel')

        #fig.title(f.filename)

        if values.plot_hist:
            '''
            for c in xrange(nchans):
                hax.hist(data[:, c].flatten(), label='IF %d'%ifn, bins=100, histtype='step')
            '''

            #bins=np.arange(0, 25r5, 1)
            #bins=256
            #hax.hist(data[0::2, :].flatten(), bins=bins, histtype='step', color='k')
            #hax.hist(data[1::2, :].flatten(), bins=bins, histtype='step', color='r')
            hax.scatter(data.flatten(), beam_mean.flatten(), marker='.')

        if values.period > 0:
            vfolded = fold(data, f.tsamp, values.period)
            print('Plotting folded data')
            fax.imshow(vfolded.T, label='IF %s'%ifn, aspect='auto')
            

        T = ntimes*f.tsamp*1e3
        #T = float(ntimes)
        times = np.arange(ntimes)*f.tsamp
        periods = T/np.arange(ntimes, dtype=np.float)
        funit='ms'
        pltharmonics = np.arange(2, 30)*f.tsamp*1e3

        if values.plot_fft:
            for chan in range(0, nchans,16):
                afixed = (data[:, chan])
                af = abs(np.fft.rfft(afixed))
                if np.any(af[1:] > 0):
                    fftax.loglog(periods[1:len(af)], af[1:])
                    fftax.set_xlabel('period (%s)' % funit)

                #fftax.set_xlim(2,10)
                #fftax.set_ylim(1,100)

        pylab.title(f.filename)
        bpax.plot(data.mean(axis=0))
        bpax.set_xlabel('Channel')
        bpax.set_ylabel('Mean bandpass')
        bpax2 = bpax.twinx()
        bpax2.plot(data.std(axis=0), 'r:')
        bpax2.set_ylabel('Bandpass std')

        #dm0ax.plot(times, dm0);
        #dm0ax.set_xlabel('Time (ms)')
        #dm0ax.set_ylabel('DM=0 time series')
        
        dm0f = np.abs(np.fft.rfft(dm0))**2

        dm0ax.semilogy(periods[1:len(dm0f)], dm0f[1:])
        dm0ax.set_xlabel('Period (%s)' % funit)
        dm0ax.set_ylabel('FFT(DM0)**2')
        
        for harmonic in pltharmonics:
            dm0ax.vlines(harmonic, min(dm0f[1:]), max(dm0f), color='k')

        #dm0ax.set_xlim(0, 30)
        #dm0ax.set_ylim(1e2, 1e9)
        
        dm0sf = np.abs(np.fft.rfft(dm0))**2
        dm0axs.semilogy(periods[1:len(dm0sf)], dm0sf[1:])
        dm0axs.set_xlabel('Period (%s)' % funit)
        dm0axs.set_ylabel('FFT(std(DM0))**2')
        #dm0axs.set_xlim(0, 30)

            
    fig.savefig('{}_dynspec.png'.format(f.filename))

    #for i in xrange(nchans):
        #ifv = v[:, i, 0]
   #     print 'Channel 0 mean, if=', i, ifv.mean(), ifv.std()
        #ifv = v[:, i, 72]
   #     print 'Channel 72 mean, if=', i, ifv.mean(), ifv.std()


    #mysubplots_adjust(wspace=0, hspace=0)
    new_tick_locations = np.array([.2, .5, .9])

    def tick_function(X):
        V = 1/(1+X)
        return ["%.3f" % z for z in V]

    #dm0axy = dm0ax.twiny()
    #ax2.set_xlim(ax1.get_xlim())
    #ax2.set_xticks(new_tick_locations)
    ##ax2.set_xticklabels(tick_function(new_tick_locations))
    #ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")


    pylab.show()

    
def plot_cov(v, f, save_all=False):
    vf = np.zeros_like(v)
    vmax = 0
    vmin = -3

    ntimes, nbeams, nchans = v.shape

    for c in range(nchans):
        R = np.cov(v[:, :, c].T)
        lam, evec = np.linalg.eig(R)
        evec = np.matrix(evec)
        for t in range(v.shape[0]):
            vf[t, :, c] -= np.dot(v[t, :, c], evec)*lam

        pylab.clf()
        pylab.imshow(np.log10(R/np.diag(R)), interpolation='none', vmax=vmax, vmin=vmin)
        pylab.xlabel('Beam number')
        pylab.ylabel('Beam number')
        pylab.title('{} channel {}'.format(f.filename, c))
        cb = pylab.colorbar()
        cb.set_label('Correlation coefficient (dB)')

        if save_all or c == 0:
            pylab.savefig('{}_beamr_c{:03d}.png'.format(f.filename, c))
        #pylab.show()



def plot_beam_cov(v, f):
    nchans = f.nchans
    nchans = 48
    nbeams = nrows*ncols
    v = v[:, 0:nbeams, 32:nchans+32]
    #v = np.swapaxes(v, 1, 2)
    vr = v.reshape(v.shape[0], v.shape[1]*v.shape[2])
    print(vr.shape)
    vrr = np.cov(vr.T)
    pylab.imshow(np.log10(vrr/np.diag(vrr)), vmax=0, vmin=-2, aspect='auto', origin='upper', interpolation='none')
    c = pylab.colorbar()
    c.set_label('Correlation coefficient (dB)')
    pylab.xlabel('beam x channel')
    pylab.ylabel('beam x channel')

    

if __name__ == '__main__':
    _main()
