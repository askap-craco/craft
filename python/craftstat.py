#!/usr/bin/env python
"""
Calculates statistics on lots of beams

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import itertools
from craftobs import load_beams
from plotutil import subplots
import matplotlib.gridspec as gridspec
from influxdb import InfluxDBClient

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def annotate(fig, title, xlabel, ylabel):
    fig.text( 0.5, 0.98,title, ha='center', va='top')
    fig.text(0.5, 0.02, xlabel, ha='center', va='bottom')
    fig.text(0.02, 0.5, ylabel, rotation=90, ha='center', va='top')

def commasep(s):
    return map(int, s.split(','))

def floatcommasep(s):
    return map(float, s.split(','))

def next_divisor(x, n):
    for t in xrange(x, n/2):
        if n % t == 0:
            return t

def divisors(n):
    d = [t for t in xrange(1, n/2 + 1) if n % t == 0]
    return d

MJD_UNIX_EPOCH = 40587.0

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('--step', type=int, help='time step', default=10)
    parser.add_argument('--ntimes', type=int, help='Numerb of samples per block', default=1024)
    parser.set_defaults(verbose=False, nxy="1,1")
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    st = CraftStatMaker(values.files, values)
    #influxout = open('/tmp/craftstat.influxdb', 'w')
    fullpath  = os.path.abspath(values.files[0])
    pathbits = fullpath.split('/') # hacky way of getting sbid, scanid and ant
    print len(pathbits), pathbits
    sbid = pathbits[-5]
    scanid = pathbits[-4]
    ant = pathbits[-3]
    startid = pathbits[-2]
    client = InfluxDBClient(host='akingest01', database='craft')
    while True:
        s, unix_start = st.next_stat()
        nbeams, nstat = s.shape
        # nanoseconds since unix epoch
        unix_nanosec = int(np.round(unix_start * 1e9))
        for b in xrange(nbeams):
            outs = 'craftstat,sbid={},scanid={},ant={},beam={} '.format(sbid, scanid, ant, b)
            statbits = ['{}={}'.format(statname, stat) for (statname, stat) in zip(st.stat_names, s[b,: ])]
            outs += ','.join(statbits)
            outs += ' {}\n'.format(unix_nanosec)
            field_dict = {}
            for sname, v in zip(st.stat_names, s[b, :]):
                field_dict[sname] = v
            #influxout.write(outs)
            body = {'measurement':'craftstat',
                    'tags':{'sbid':sbid,'scanid':scanid,'ant':ant,'beam':b},
                    'time':unix_nanosec,
                    'fields':field_dict
                    }
            client.write_points([body])
            #print outs

class CraftStatMaker(object):
    def __init__(self, files, values):
        self.files = files
        self.values = values
        self.ntimes = values.ntimes
        self.tstep = values.step
        self.tstart = 0
        self.fft_freqs = np.array([38.5, 50, 100, 200, 300])
        beams, self.spfiles = load_beams(self.files, self.tstart, self.ntimes, return_files=True)
        self.stat_names = ['bmean','bstd']
        for f in self.fft_freqs:
            self.stat_names.append('f{}'.format(f))

    def next_data(self):

        for ifin, f in enumerate(self.spfiles):
            f.seek_sample(self.tstart)
            nelements = f.nchans*self.ntimes
            dtype = np.uint8
            v = np.fromfile(f.fin, dtype=dtype, count=nelements )

            

    def next_stat(self):
        beams, files = load_beams(self.files, self.tstart, self.ntimes, return_files=True)
        mjdstart = files[0].tstart
        tsamp = files[0].tsamp

        self.mjdstart = mjdstart
        self.tsamp = tsamp
        mjdtime = self.mjdstart + self.tstart * self.tsamp/3600./24.
        unix_start = (self.mjdstart - MJD_UNIX_EPOCH)*86400.
        unix_time = unix_start + self.tstart * self.tsamp 
        self.tstart += self.tstep * self.ntimes

        # Normalise to nominal 0 mean, unit stdev - HACK!
        beams -= 128
        beams /= 18

        ntimes, nbeams, nfreq = beams.shape
        stat = np.zeros((nbeams, 2 + len(self.fft_freqs)))

        for i in xrange(nbeams):
            bi = beams[:, i, :]
            ntimes, nfreq = bi.shape
            # spectra
            beam_mean = bi.mean(axis=0)
            beam_std = bi.std(axis=0)
            #bmm = np.tile(beam_mean, (ntimes, 1))
            #bsm = np.tile(beam_std, (ntimes, 1))
            #bi_znorm = (bi - bmm)/bsm
            #beam_kurtosis = np.mean((bi_znorm)**4, axis=0)/np.mean((bi_znorm)**2, axis=0)**2
            bstd = beam_std.std()
            stat[i, 0] = beam_mean.mean()
            stat[i, 1] = bstd

            # dm0
            dm0 = bi.mean(axis=1)
            dm0f = abs(np.fft.rfft(dm0, axis=0))**2
            ntimes, nchans = bi.shape

            idxs = np.rint(self.fft_freqs * float(ntimes) * self.tsamp).astype(int)
            stat[i, 2:] = dm0f[idxs]/bstd

            
            fftfreqs  = np.arange(len(dm0f))/float(ntimes)/self.tsamp

            #pylab.plot(fftfreqs, dm0f)
            #for f in self.fft_freqs:
            #    pylab.axvline(f)
            #pylab.show()

        return stat, unix_time


def getstats():
    tstart = self.tstart
    ntimes = self.ntimes
    beams, files = load_beams(self.files, tstart, ntimes, return_files=True)


class Plotter(object):
    @staticmethod

    def from_values(values, tstart, ntimes):
        p = Plotter(values.files, values.nxy, fft=values.fft, raw_units=values.raw_units)
        if values.seconds:
            p.set_position_seconds(*values.seconds)
        else:
            p.set_position_sample(tstart, ntimes)

        p.imzrange = values.imzrange
        p.dm = values.dm
        p.draw()

        return p
    
    def __init__(self, filenames, nxy, tstart=0, ntimes=1024, fft=False, raw_units=False):
        self.nrows, self.ncols = nxy
        self.figs = {}
        self.fig_labels = {}
        # Sniff data
        self.files = filenames
        beams, files = load_beams(filenames, tstart, ntimes=1, return_files=True)
        ntimes, self.nbeams, self.nfreq = beams.shape

        self.bnames = [f.filename.split('.')[-2] for f in files]
        nbeams = len(self.files)
        mjdstart = files[0].tstart
        tsamp = files[0].tsamp
        self.raw_units = raw_units
        self.tscrunch_factor = 1
        self.fscrunch_factor = 1
        self.dm = 0.
        if self.raw_units:
            xunit = 'sample'
            yunit = 'channel'
        else:
            xunit= 'sec'
            yunit = 'MHz'
            
        self.mjdstart = mjdstart
        self.tsamp = tsamp
        self.set_position_sample(tstart, ntimes)
        self.imzrange = None
        if nbeams == 1:
            self.mk_single_fig('dynspec', '', 'Time (%s) after %f' % (xunit, mjdstart), 'Frequency (%s)' % yunit)
        else:
            self.mkfig('mean', 'Mean bandpass', 'Frequency (%s)' % yunit,'Mean bandpass')
            self.mkfig('std', 'Bandpass stdDev', 'Frequency (%s)' % yunit ,'Bandpass StdDev')
            self.mkfig('kurt','Kurtosis', 'Frequency (%s)' % yunit, 'Kurtosis')
            self.mkfig('dynspec', 'Dynamic spectrum', 'Time (%s) after %f' % (xunit, mjdstart), 'Frequency (%s)' % yunit)

        if nbeams > 1:
            self.mkfig('cov', 'Beam covariance', 'beam no', 'beamno', 1,1)

        if fft:
            self.mkfig('fftim', 'FFT of all channels', 'Frequency (Hz)', 'Channel')
            self.mkfig('fftplt', 'FFT of DM0', 'Frequency (Hz)', 'FFT (dB)')
            self.mkfig('dm0plt', 'DM0', 'Sample', 'Intensity')

        self.fft = fft

    def set_position_sample(self, tstart, ntimes):
        self.tstart = tstart
        self.ntimes = ntimes

    def set_position_seconds(self, sstart, secs):
        self.tstart = int(sstart/self.tsamp)
        self.ntimes = int(secs/self.tsamp)

    def mk_single_fig(self, name, title, xlab, ylab):
        p = plt.subplot
        fig = plt.figure()
        gs = gridspec.GridSpec(3,3)
        rawax = fig.add_subplot(gs[0:2, 0:2])
        tax = fig.add_subplot(gs[2,0:2], sharex=rawax)
        fax = fig.add_subplot(gs[0:2,2], sharey=rawax)
        fig.canvas.mpl_connect('key_press_event', self.press)
        #annotate(fig, title, xlab, ylab)
        self.figs[name] = (fig, [rawax, tax, fax])
        self.fig_labels[name] = (title, xlab, ylab)

    def mkfig(self, name, title, xlab, ylab, nrows=None, ncols=None):

        if nrows is None:
            nrows = self.nrows
        if ncols is None:
            ncols = self.ncols

        fig, axes = subplots(nrows, ncols, sharex=True, sharey=True)
        axes = axes.flatten()
        fig.canvas.mpl_connect('key_press_event', self.press)
        annotate(fig, title, xlab, ylab)

        self.figs[name] = (fig, axes)

    def getfig(self, name, i):
        fig, axes = self.figs[name]
        
        ax = axes[i]
        print i, len(axes), len(self.bnames), self.bnames
        ax.text(0.98, 0.98, 'B{}'.format(self.bnames[i]), va='top', ha='right', transform=ax.transAxes)
        return fig, ax

    def clearfigs(self):
        for name, (fig, axes) in self.figs.iteritems():
            for ax in axes:
                print 'clearing', fig, ax
                ax.cla()
                

    def saveall(self, prefix):
        for name, (fig, axes) in self.figs.iteritems():
            fout='{}_{}.png'.format(prefix, name)
            print 'Saving', fout
            fig.savefig(fout)

    def closeall(self):
        for name, (fig, axes) in self.figs.iteritems():
            plt.close(fig)

        self.figs = {}
        
    def drawall(self):
        for name, (fig, axes) in self.figs.iteritems():
            fig.draw()


    def draw(self):
        f0 = files[0]
        self.beams = beams
        print 'Loaded beams', beams.shape
        print 'scrunching t=', self.tscrunch_factor, 'f=', self.fscrunch_factor, 'dm', self.dm

        if self.dm != 0:
            beams = dmroll(beams, self.dm, f0.fch1/1e3, f0.foff/1e3, f0.tsamp*1e3)

        if self.tscrunch_factor != 1:
            beams = tscrunch(beams, self.tscrunch_factor)

        if self.fscrunch_factor != 1:
            beams = fscrunch(beams, self.fscrunch_factor)

        ntimes, nbeams, nfreq = beams.shape
        mjdstart = f0.tstart
        
        if self.raw_units:
            tsamp = 1.
            fch1 = 0.
            foff = 1.
        else:
            tsamp = f0.tsamp
            fch1 = f0.fch1
            foff = f0.foff

        src_raj = f0.src_raj
        src_dej = f0.src_dej
        freqs = np.linspace(fch1, fch1 + nfreq*foff, nfreq, endpoint=True)
        times=  np.linspace(tstart*tsamp, (ntimes+tstart)*tsamp, ntimes, endpoint=True)

        assert(len(freqs) == nfreq)
        
        if foff < 0:
            origin = 'upper'
            im_extent = (times[0], times[-1], freqs[-1], freqs[0])
        else:
            origin = 'lower'
            im_extent = (times[0], times[-1], freqs[0], freqs[-1])


        if self.imzrange is None:
            imzmin = None
            imzmax = None
        else:
            imzmin, imzmax = self.imzrange

        if nbeams > 1:
            self.draw_many(beams, im_extent, origin, imzmin, imzmax, freqs)
        else:
            self.draw_single(beams, im_extent, origin, imzmin, imzmax, freqs, times)

        pylab.draw()

    def draw_single(self, beams, im_extent, origin, imzmin, imzmax, freqs, times):
        bi = beams[:, 0, :]
        ntimes, nfreq = bi.shape
        bi = bi.T
        print 'BISHAPE', bi.shape, 'ZRAGE', imzmin, imzmax
        fig, [rawax, tax, fax] = self.figs['dynspec']
        rawax.imshow(bi, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent, interpolation='none')
        if imzmin is None and imzmax is None:
            self.imzrange = (bi.min(), bi.max())
            
        fax.plot(bi.mean(axis=1)*np.sqrt(ntimes), freqs, label='mean')
        #fax.plot(bi.std(axis=1)*np.sqrt(ntimes),freqs, label='std')
        fax.set_ylim(freqs.min(), freqs.max())
        tax.plot(times, bi.mean(axis=0)*np.sqrt(nfreq), label='mean')
        #tax.plot(times, bi.std(axis=0)*np.sqrt(nfreq), label='std')
        tax.set_xlim(times.min(), times.max())

        title, xlab, ylab = self.fig_labels['dynspec']
        rawax.set_ylabel(ylab)
        plt.setp(rawax.get_xticklabels(), visible=False)
        plt.setp(fax.get_yticklabels(), visible=False)
        fig.subplots_adjust(hspace=0, wspace=0)
        tax.set_xlabel(xlab)
        tax.set_ylabel('S/N')
        fax.set_xlabel('Bandpass')
        #tax.legend(frameon=False)
        #fax.legend(frameon=False)


    def draw_many(self, beams, im_extent, origin, imzmin, imzmax, freqs):
        ntimes, nbeams, nfreq = beams.shape
        chan = 150
        dm0 = beams.mean(axis=2)
        #dm0 = beams[:, :, chan]
        dm0 -= dm0.mean(axis=0)

        if nbeams > 1:
            cov = np.cov(dm0.T)
            #b0 = beams[:, 0, :]
            #cov = np.cov(b0.T)
            cov /= np.diag(cov)
            fig, ax = self.figs['cov']
            ax[0].imshow(np.log10(cov), interpolation='none')
            ax[0].set_xlabel('Beam number')
            ax[0].set_ylabel('Beam number')
            
        nplots = min(self.nrows*self.ncols, nbeams)
        for i in xrange(nplots):

            fig1, ax1 = self.getfig('dynspec', i)
            fig2, ax2 = self.getfig('mean', i)
            fig3, ax3 = self.getfig('std', i)
            fig6, ax6 = self.getfig('kurt', i)
            bi = beams[:, i, :]
            print 'bi', bi.shape
            ntimes, nfreq = bi.shape
        
            ax1.imshow(bi.T, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent, interpolation='none')

            beam_mean = bi.mean(axis=0)
            beam_std = bi.std(axis=0)
            bmm = np.tile(beam_mean, (ntimes, 1))
            bsm = np.tile(beam_std, (ntimes, 1))
            bi_znorm = (bi - bmm)/bsm
            print 'Znorm', bi_znorm.shape
            beam_kurtosis = np.mean((bi_znorm)**4, axis=0)/np.mean((bi_znorm)**2, axis=0)**2
            print 'kurt shape', beam_kurtosis.shape

            ax2.plot(freqs, beam_mean)
            ax3.plot(freqs, beam_std)
            ax6.plot(freqs, beam_kurtosis)
            
            if self.fft:
                dm0 = bi.mean(axis=1)
                fig4, ax4 = self.getfig('fftim', i)
                fig5, ax5 = self.getfig('fftplt', i)
                dm0f = abs(np.fft.rfft(dm0, axis=0))**2
                ntimes, nchans = bi.shape
                bf = abs(np.fft.rfft(bi, axis=0).T)**2
                fftfreqs  = np.arange(len(dm0f))/float(ntimes)/self.tsamp
                fft_ext = (fftfreqs.min(), fftfreqs.max(), 0, nchans)
                ax4.imshow(np.log10(bf)[:, 1:], aspect='auto', extent=fft_ext, origin='lower')
                print 'fft', dm0f.shape
                ax5.loglog(fftfreqs[1:], (dm0f[1:]))
                fig6, ax6 = self.getfig('dm0plt', i)
                ax6.plot(dm0)

        #self.drawall()

if __name__ == '__main__':
    _main()
