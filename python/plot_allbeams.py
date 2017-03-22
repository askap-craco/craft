#!/usr/bin/env python
"""
Plots a lot of beams

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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def onpick(event):
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print thisline.get_label(), xdata[ind], ydata[ind]

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

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-t', '--times', help='Integration range to plot (samples)', type=commasep)
    parser.add_argument('-s', '--seconds', help='Integration range to plot (seconds)', type=commasep)
    parser.add_argument('--nxy', help='number of rows,columns in plots', type=commasep)
    parser.add_argument('--imzrange', help='Z range for dynamic spectrum', type=floatcommasep)
    parser.add_argument('--fft', help='plot fft', action='store_true',default=False)
    parser.add_argument('--save', help='Save plots as png', action='store_true', default=False)
    parser.add_argument('--raw-units', help='Use raw unts, rather than physical units on axis', action='store_true', default=False)
    parser.add_argument('--dm', help='Dispersion measure (pc/cm3)', default=0., type=float)
    parser.set_defaults(verbose=False, nxy="1,1")
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.times:
        bits = values.times
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

    plt = Plotter.from_values(values, tstart, ntimes)
    pylab.show()

def tscrunch(beams, factor):
    ntimes, nbeams, nfreq = beams.shape
    newbeams = np.zeros((ntimes/factor, nbeams, nfreq))
    for t in xrange(newbeams.shape[0]):
        tstart = t*factor
        tend = (t+1)*factor
        newbeams[t, :, :] = beams[tstart:tend, :, :].sum(axis=0)/np.sqrt(factor)

    return newbeams


def fscrunch(beams, factor):
    ntimes, nbeams, nfreq = beams.shape
    newbeams = np.zeros((ntimes, nbeams, nfreq/factor))
    for f in xrange(newbeams.shape[2]):
        fstart = f*factor
        fend = (f+1)*factor
        newbeams[:, :, f] = beams[:, :, fstart:fend].sum(axis=2)/np.sqrt(factor)

    return newbeams


def dmroll(beams, dm, fch1, foff, tint):
    newbeams = np.empty_like(beams)
    ntimes, nbeams, nfreq = newbeams.shape
    for f in xrange(nfreq):
        freq = fch1 + f*foff
        tdelay = 4.15*dm*(fch1**-2 - freq**-2)
        shift = int(np.round(tdelay/tint))
        newbeams[:, :, f] = np.roll(beams[:,:,f], shift)

    return newbeams
        
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
        print self.files[0]
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


    def __del__(self):
        self.closeall()
    
    def press(self, event):
        print 'press', event.key
        draw = True
        if event.key == 'right' or event.key == 'n':
            self.tstart += self.ntimes/2
        elif event.key == 'left' or event.key == 'p':
            self.tstart -= self.ntimes/2
            self.tstart = max(self.tstart, 0)
        elif event.key == 'w':
            self.ntimes *= 2
        elif event.key == 'a':
            if self.ntimes > 2:
                self.ntimes /= 2
        elif event.key == 't':
            self.tscrunch_factor += 1
        elif event.key == 'T':
            self.tscrunch_factor = max(1, self.tscrunch_factor - 1)
        elif event.key == 'f':
            fdiv = divisors(self.nfreq)
            self.fscrunch_factor = fdiv[fdiv.index(self.fscrunch_factor) + 1]
        elif event.key == 'F':
            fdiv = divisors(self.nfreq)
            self.fscrunch_factor = fdiv[fdiv.index(self.fscrunch_factor) -1]
        elif event.key == 'd':
            self.dm = float(raw_input('Input DM(pc/cm3)'))
        elif event.key == 'ctrl+c':
            sys.exit(0)
        else:
            draw = False

        if draw:
            self.clearfigs()
            self.draw()


    def draw(self):
        tstart = self.tstart
        ntimes = self.ntimes
        beams, files = load_beams(self.files, tstart, ntimes, return_files=True)
        beams -= 128
        beams /= 18
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
        print 'BISHAPE', bi.shape
        fig, [rawax, tax, fax] = self.figs['dynspec']
        rawax.imshow(bi, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent, interpolation='none')
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
            
            dm0 = bi.mean(axis=1)
            
            if self.fft:
                fig4, ax4 = self.getfig('fftim', i)
                fig5, ax5 = self.getfig('fftplt', i)
                dm0f = abs(np.fft.rfft(dm0, axis=0))**2
                ntimes, nchans = bi.shape
                bf = abs(np.fft.rfft(bi, axis=0).T)**2
                fftfreqs  = np.arange(len(dm0f))/float(ntimes)/tsamp
                fft_ext = (fftfreqs.min(), fftfreqs.max(), 0, nchans)
                ax4.imshow(np.log10(bf)[:, 1:], aspect='auto', extent=fft_ext, origin='lower')
                ax5.loglog(fftfreqs[1:], (dm0f[1:]))
                fig6, ax6 = self.getfig('dm0plt', i)
                ax6.plot(dm0)

        #self.drawall()

if __name__ == '__main__':
    _main()
