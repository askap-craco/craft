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

        
class Plotter(object):
    @staticmethod

    def from_values(values, tstart, ntimes):
        p = Plotter(values.files, values.nxy, fft=values.fft)
        if values.seconds:
            p.set_position_seconds(*values.seconds)
        else:
            p.set_position_sample(tstart, ntimes)

        p.imzrange = values.imzrange
        p.draw()

        return p
    
    def __init__(self, filenames, nxy, tstart=0, ntimes=1024, fft=False):
        self.nrows, self.ncols = nxy
        self.figs = {}
        # Sniff data
        self.files = filenames
        beams, files = load_beams(filenames, tstart, ntimes=1, return_files=True)
        self.bnames = [f.filename.split('.')[-2] for f in files]
        nbeams = len(self.files)
        mjdstart = files[0].tstart
        tsamp = files[0].tsamp

        self.mjdstart = mjdstart
        self.tsamp = tsamp
        self.set_position_sample(tstart, ntimes)
        self.imzrange = None

        self.mkfig('mean', 'Mean bandpass', 'Frequency (MHz)','Mean bandpass')
        self.mkfig('std', 'Bandpass stdDev', 'Frequency (MHz)','Bandpass StdDev')
        self.mkfig('kurt','Kurtosis', 'Frequency (MHz)', 'Kurtosis')
        self.mkfig('dynspec', 'Dynamic spectrum', 'Time (s) after %f' % mjdstart, 'Frequency (MHz)')
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
        self.beams = beams
        print 'Loaded beams', beams.shape
        ntimes, nbeams, nfreq = beams.shape

        f0 = files[0]
        mjdstart = f0.tstart
        tsamp = f0.tsamp
        fch1 = f0.fch1
        foff = f0.foff
        src_raj = f0.src_raj
        src_dej = f0.src_dej
        freqs = np.linspace(fch1, fch1 + nfreq*foff, nfreq, endpoint=True)
        assert(len(freqs) == nfreq)
        
        if foff < 0:
            origin = 'upper'
            im_extent = (tstart*tsamp, (ntimes+tstart)*tsamp, fch1 + foff*nfreq, fch1)
        else:
            origin = 'lower'
            im_extent = (tstart*tsamp , (ntimes+tstart) *tsamp, fch1, fch1 + foff*nfreq)

        chan = 150
        dm0 = beams.mean(axis=2)
        dm0 = beams[:, :, chan]
        dm0 -= dm0.mean(axis=0)

        if self.imzrange is None:
            imzmin = None
            imzmax = None
        else:
            imzmin, imzmax = self.imzrange


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

        pylab.draw()
        #self.drawall()

if __name__ == '__main__':
    _main()
