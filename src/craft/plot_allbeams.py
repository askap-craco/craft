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
from .craftobs import load_beams
from .plotutil import subplots
import matplotlib.gridspec as gridspec

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def onpick(event):
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print(thisline.get_label(), xdata[ind], ydata[ind])

def annotate(fig, title, xlabel, ylabel):
    fig.text( 0.5, 0.98,title, ha='center', va='top')
    fig.text(0.5, 0.02, xlabel, ha='center', va='bottom')
    fig.text(0.02, 0.5, ylabel, rotation=90, ha='center', va='top')

def commasep(s):
    return list(map(int, s.split(',')))

def floatcommasep(s):
    return list(map(float, s.split(',')))

def next_divisor(x, n):
    for t in range(x, n/2):
        if n % t == 0:
            return t

def divisors(n):
    d = [t for t in range(1, n/2 + 1) if n % t == 0]
    return d

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-t', '--times', help='Integration range to plot (samples)', type=commasep)
    parser.add_argument('-s', '--seconds', help='Integration range to plot (seconds)', type=commasep)
    parser.add_argument('-m','--mjd', help='Integration range to plot(mjd,sec)', type=floatcommasep)
    parser.add_argument('--nxy', help='number of rows,columns in plots', type=commasep)
    parser.add_argument('--imzrange', help='Z range for dynamic spectrum', type=floatcommasep)
    parser.add_argument('--fft', help='plot fft', action='store_true',default=False)
    parser.add_argument('--save', help='Save plots as png', action='store_true', default=False)
    parser.add_argument('--raw-units', help='Use raw unts, rather than physical units on axis', action='store_true', default=False)
    parser.add_argument('-d', '--dm', help='Dispersion measure (pc/cm3)', default=0., type=float)
    parser.add_argument('-T', '--tscrunch', help='TScrunch by this factor', default=1, type=int)
    parser.add_argument('-F','--fscrunch', help='Fscrunch by this factor', default=1, type=int)
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
    for t in range(newbeams.shape[0]):
        tstart = t*factor
        tend = (t+1)*factor
        newbeams[t, :, :] = beams[tstart:tend, :, :].sum(axis=0)/np.sqrt(factor)

    return newbeams


def fscrunch(beams, factor):
    ntimes, nbeams, nfreq = beams.shape
    newbeams = np.zeros((ntimes, nbeams, nfreq/factor))
    for f in range(newbeams.shape[2]):
        fstart = f*factor
        fend = (f+1)*factor
        newbeams[:, :, f] = beams[:, :, fstart:fend].sum(axis=2)/np.sqrt(factor)

    return newbeams


def dmroll(beams, dm, fch1, foff, tint):
    newbeams = np.empty_like(beams)
    ntimes, nbeams, nfreq = newbeams.shape
    reffreq = min(fch1, fch1 + nfreq*foff)
    for f in range(nfreq):
        freq = fch1 + f*foff
        tdelay = 4.15*dm*(reffreq**-2 - freq**-2)
        shift = int(np.round(tdelay/tint))
        newbeams[:, :, f] = np.roll(beams[:,:,f], shift, axis=0)

    return newbeams
        
class Plotter(object):
    @staticmethod

    def from_values(values, tstart, ntimes):
        p = Plotter(values.files, values.nxy, fft=values.fft, raw_units=values.raw_units, tscrunch=values.tscrunch, fscrunch=values.fscrunch)
        if values.seconds:
            p.set_position_seconds(*values.seconds)
        elif values.mjd:
            p.set_position_mjd(*values.mjd)
        else:
            p.set_position_sample(tstart, ntimes)

        p.imzrange = values.imzrange
        p.dm = values.dm
        p.draw()

        return p
    
    def __init__(self, filenames, nxy, tstart=0, ntimes=1024, fft=False, raw_units=False, tscrunch=1, fscrunch=1):
        self.nrows, self.ncols = nxy
        self.figs = {}
        self.fig_labels = {}
        # Sniff data
        self.files = filenames
        print(self.files[0])
        beams, files = load_beams(filenames, tstart, ntimes=1, return_files=True)
        self.total_samples = min([f.nsamples for f in files])
        ntimes, self.nbeams, self.nfreq = beams.shape
        self.freq_flags = np.ones(self.nfreq, dtype=np.bool)
        f0 = files[0]
        self.freqs = np.linspace(f0.fch1, f0.fch1 + self.nfreq*f0.foff, self.nfreq, endpoint=True)
        self.bnames = [f.filename.split('.')[-2] for f in files]
        nbeams = len(self.files)
        mjdstart = files[0].tstart
        tsamp = files[0].tsamp
        self.raw_units = raw_units
        self.tscrunch_factor = tscrunch
        self.fscrunch_factor = fscrunch
        self.dm = 0.
        if self.raw_units:
            xunit = 'sample'
            yunit = 'channel'
        else:
            xunit= 'sec'
            yunit = 'MHz'
            
        self.mjdstart = mjdstart
        self.tsamp = tsamp
        self.rescale = False
        self.set_position_sample(tstart, ntimes)
        self.imzrange = None
        self.histfig = None
        if nbeams == 1:
            self.mk_single_fig('dynspec', '', 'Time (%s) after %f' % (xunit, mjdstart), 'Frequency (%s)' % yunit)
        else:
            self.mkfig('mean', 'Mean bandpass', 'Frequency (%s)' % yunit,'Mean bandpass')
            self.mkfig('std', 'Bandpass stdDev', 'Frequency (%s)' % yunit ,'Bandpass StdDev')
            self.mkfig('kurt','Kurtosis', 'Frequency (%s)' % yunit, 'Kurtosis')
            self.mkfig('dynspec', 'Dynamic spectrum', 'Time (%s) after %f' % (xunit, mjdstart), 'Frequency (%s)' % yunit)
            self.mkfig('time','Time series', 'Time (%s) after %f' % (xunit, mjdstart), 'Amplitude')

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

    def set_position_mjd(self, mjdstart, ntimes):
        self.ntimes = int(ntimes)
        self.tstart = int(((mjdstart - self.mjdstart)*86400)/self.tsamp) - self.ntimes # Is my tstart really what I think it is?
        print('MJD', mjdstart, ntimes, self.tstart)

    def mk_single_fig(self, name, title, xlab, ylab):
        p = plt.subplot
        fig = plt.figure()
        gs = gridspec.GridSpec(3,3)
        rawax = fig.add_subplot(gs[0:2, 0:2])
        tax = fig.add_subplot(gs[2,0:2], sharex=rawax)
        fax = fig.add_subplot(gs[0:2,2], sharey=rawax)
        fig.canvas.mpl_connect('key_press_event', self.press)
        fig.canvas.mpl_connect('pick_event', self.pick)
        #fig.canvas.mpl_connect('button_press_event', self.button_press)
        #fig.canvas.mpl_connect('button_release_event', self.button_release)
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
        print(i, len(axes), len(self.bnames), self.bnames)
        ax.text(0.98, 0.98, 'B{}'.format(self.bnames[i]), va='top', ha='right', transform=ax.transAxes)
        return fig, ax

    def clearfigs(self):
        for name, (fig, axes) in self.figs.items():
            for ax in axes:
                print('clearing', fig, ax)
                ax.cla()
                

    def saveall(self, prefix):
        for name, (fig, axes) in self.figs.items():
            fout='{}_{}.png'.format(prefix, name)
            print('Saving', fout)
            fig.savefig(fout)

    def closeall(self):
        for name, (fig, axes) in self.figs.items():
            plt.close(fig)

        self.figs = {}
        
    def drawall(self):
        for name, (fig, axes) in self.figs.items():
            fig.draw()

    def __del__(self):
        self.closeall()

    def pick(self, event):
        ''' Called when data picked'''
        #        print event, dir(event)
        #print event.artist, dir(event.artist)
        #print event.mouseevent, dir(event.mouseevent)
        e = event.mouseevent
        time, freq = e.xdata, e.ydata
        self.toggle_freq_flags(freq)
        self.clearfigs()
        self.draw()

    def toggle_freq_flags(self, f1, f2=None):
        f1idx = np.argmin(np.round(abs(self.freqs - f1)))
        if f2 is None:
            f2idx = f1idx + 1
        else:
            f2idx = np.argmin(np.round(abs(self.freqs - f2)))

        if f1idx > f2idx:
            f1idx, f2idx = (f2idx, f1idx)

        assert f1idx < f2idx
        print('Toggling frequencies', f1, f2, f1idx, f2idx)
        self.freq_flags[f1idx:f2idx] = ~ self.freq_flags[f1idx:f2idx]
    
    def reset_flagged_frequencies(self):
        self.freq_flags[:] = True
        self.draw()

    def flag_frequencies(self):
        try:
            l = input('f1 f2 (with space):')
            f1, f2 = list(map(float, l.split()))
            self.toggle_freq_flags(f1, f2)
        except Exception as e:
            print('Couldnt parse flags', l, e)


    
    def press(self, event):
        print('press', event.key)
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
        elif event.key == 'z':
            fdiv = divisors(self.nfreq)
            self.fscrunch_factor = fdiv[fdiv.index(self.fscrunch_factor) + 1]
        elif event.key == 'Z':
            fdiv = divisors(self.nfreq)
            self.fscrunch_factor = fdiv[fdiv.index(self.fscrunch_factor) -1]
        elif event.key == 'd':
            self.dm = float(input('Input DM(pc/cm3)'))
        elif event.key == 'c':
            self.squeeze_zrange(2.)
        elif event.key == 'C':
            self.squeeze_zrange(0.5)
        elif event.key == 'r':
            self.rescale = not self.rescale
            self.imzrange = None
        elif event.key == 'g':
            self.flag_frequencies()
        elif event.key == 'G':
            self.reset_flagged_frequencies()
        elif event.key == 'e':
            self.goto_end()
        elif event.key == 'b':
            self.goto_beginning()
        elif event.key == 'ctrl+c':
            sys.exit(0)
        elif event.key == 'h':
            self.show_histogram()
        elif event.key == '?':
            self.print_help()
            draw = False
        else:
            draw = False

        if draw:
            self.clearfigs()
            self.draw()

    def goto_start(self):
        self.tstart = 0

    def goto_end(self):
        self.tstart = self.total_samples - self.ntimes

    def show_histogram(self):
        if self.histfig is None:
            self.histfig = self.mkfig('histfig', 'Histograms', 'S/N','Number',1,1)

        #fig, ax = self.getfig('histfig', 0)
        fig,ax = pylab.subplots(1,1)
        beams = self.beams
        bi = beams[:, 0, :]
        bif = bi.flatten()
        ntimes, nfreq = bi.shape
        ax.hist(bif, bins=100)
        fig.show()

    def squeeze_zrange(self, mul):
        zmin, zmax = self.imzrange
        zmid = (zmin + zmax)/2.
        zrange = zmax - zmin
        new_zrange = zrange/mul
        self.imzrange = (zmid - new_zrange/2., zmid + new_zrange/2.)
        return self.imzrange


    def print_help(self):
        s = '''
        Key Mapping
        n or right arrow - Move right by half a window
        p or left arrow - Move left by half a window
        w - zoom out by 2
        a - zoom in by 2
        t - increase tscrunch by 1 bin
        T - decrease tscrunch by 1 bin
        z - increase fscrunch by 1 bin
        Z - decrease fscrunch by 1 bin
        d - Dedisperse (I'll ask for the DM on the cmdline
        c - Increase colormap zoom
        C - Decrease colormap zoom
        g - Set frequency flags from prompt
        G - clear frequency flags
        r - Toggle rescaling
        k - Show histogram'
        h or ? - Print this help
        Ctrl-C - quit'''
        print(s)

        return s


    def draw(self):
        tstart = self.tstart
        ntimes = self.ntimes
        beams, files = load_beams(self.files, tstart, ntimes, return_files=True)
        beams = np.ma.MaskedArray(beams, (beams == 0) | (~self.freq_flags[np.newaxis, np.newaxis, :]))
            
        f0 = files[0]
        self.beams = beams
        print('Loaded beams', beams.shape)
        print('scrunching t=', self.tscrunch_factor, 'f=', self.fscrunch_factor, 'dm', self.dm)

        orig_ntimes, orig_nbeams, orig_nfreq = beams.shape
        
        if self.dm != 0:
            beams = dmroll(beams, self.dm, f0.fch1/1e3, f0.foff/1e3, f0.tsamp*1e3)

        if self.tscrunch_factor != 1:
            beams = tscrunch(beams, self.tscrunch_factor)

        if self.fscrunch_factor != 1:
            beams = fscrunch(beams, self.fscrunch_factor)
            
        if self.rescale:
            print('Doing rescale')
            beams -= beams.mean(axis=0)
            beams /= beams.std(axis=0) 

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
        freqs = np.linspace(fch1, fch1 + nfreq*foff*self.fscrunch_factor, nfreq, endpoint=True)
        times=  np.linspace(tstart*tsamp, (ntimes+tstart)*tsamp*self.tscrunch_factor, ntimes, endpoint=True)

        assert(len(freqs) ==nfreq)
        
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
            self.draw_single(beams, im_extent, origin, freqs, times)

        pylab.draw()

    def draw_single(self, beams, im_extent, origin, freqs, times):
        bi = beams[:, 0, :]
        ntimes, nfreq = bi.shape
        bi = bi.T
        fig, [rawax, tax, fax] = self.figs['dynspec']
        if self.imzrange is None:
            self.imzrange = (bi.min(), bi.max())

        imzmin, imzmax = self.imzrange
        print('BISHAPE', bi.shape, 'ZRAGE', imzmin, imzmax)
        
        rawax.imshow(bi, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent, interpolation='none', picker=3)
        #if imzmin is None and imzmax is None:
        #self.imzrange = (bi.min(), bi.max())
            
        fax.plot(bi.mean(axis=1), freqs, label='mean')
        fax.plot(bi.max(axis=1), freqs, label='max')
        fax.plot(bi.min(axis=1), freqs, label='min')
        fax.plot(bi.std(axis=1), freqs, label='std')
        fax.plot(np.ones(nfreq), freqs, ls=':')
        #fax2 = fax.twiny()
        #fax2.plot(bi.std(axis=1),freqs, 'r', label='std')
        fax.set_ylim(freqs.min(), freqs.max())
        tax.plot(times, bi.mean(axis=0)*np.sqrt(nfreq), label='mean')
        tax.plot(times, bi.std(axis=0), label='std')
        tax.plot(times, np.ones(ntimes), ls=':')
        #tax.plot(times, bi.max(axis=0), label='max')
        #tax.plot(times, bi.min(axis=0), label='min')
        tax.set_xlim(times.min(), times.max())
        #tax2 = tax.twinx()
        #tax2.plot(times, bi.std(axis=0), 'r', label='std')


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
        for i in range(nplots):

            fig1, ax1 = self.getfig('dynspec', i)
            fig2, ax2 = self.getfig('mean', i)
            fig3, ax3 = self.getfig('std', i)
            fig6, ax6 = self.getfig('kurt', i)
            bi = beams[:, i, :]
            ntimes, nfreq = bi.shape
        
            ax1.imshow(bi.T, aspect='auto', origin=origin, vmin=imzmin, vmax=imzmax, extent=im_extent, interpolation='none')

            beam_mean = bi.mean(axis=0)
            beam_std = bi.std(axis=0)
            bmm = np.tile(beam_mean, (ntimes, 1))
            bsm = np.tile(beam_std, (ntimes, 1))
            bi_znorm = (bi - bmm)/bsm
            print('Znorm', bi_znorm.shape)
            beam_kurtosis = np.mean((bi_znorm)**4, axis=0)/np.mean((bi_znorm)**2, axis=0)**2
            print('kurt shape', beam_kurtosis.shape)

            ax2.plot(freqs, beam_mean)
            ax3.plot(freqs, beam_std)
            ax6.plot(freqs, beam_kurtosis)

            fig7,ax7 = self.getfig('time',i)
            ax7.plot(bi.mean(axis=1))
            
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
                print('fft', dm0f.shape)
                ax5.loglog(fftfreqs[1:], (dm0f[1:]))
                fig6, ax6 = self.getfig('dm0plt', i)
                ax6.plot(dm0)

        #self.drawall()

if __name__ == '__main__':
    _main()
