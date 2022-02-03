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
from .cmdline import strrange
from .rtdata import FreddaRescaleData

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Formatter(object):
    def __init__(self, ax, xdata, d, values, name, label):
        self.ax = ax
        self.xdata = xdata
        self.d = d
        self.values = values
        self.name = name
        self.label = label

    def __call__(self, x, y):
        x = min(int(np.round(x)), self.d.shape[0])
        y = min(int(np.round(y)), self.d.shape[1])
        
        z = self.d[x,y] # actually d.T is plotted
        lbl = ''
        if self.label is not None:
            if y >= 0 and y < len(self.label):
                lbl = self.label[int(y)]

        s = 'x={:0.1f} y={:0.1f} z={:0.1e} {}'.format(x, y, z, lbl)
        
        return s

def showplot(ax, x, d, values, name, label=None, **kwargs):
    im = np.ma.masked_equal(d, 0)
    
    if values.image:
        r = ax.imshow(im.T, aspect='auto', interpolation='nearest')
        ax.format_coord = Formatter(ax, x, d, values, name, label)
    else:
        (nrow, ncol) = d.shape
        # need to plot one at a time fo rlabelling. Grrr.
        r = []
        lbl = None
        for n in range(ncol):
            if label is not None:
                lbl = label[n]
                
            l = ax.plot(x, im[:, n], picker=3, label=lbl, **kwargs)
            r.append(l)
            
            
    ax.text(0.05, 1, name, ha='left', va='top', transform=ax.transAxes)


    return r

class RescalePlot(object):
    def __init__(self, f, values):
        self.thefile = f
        self.d = FreddaRescaleData(f, values.exclude_ants)
        self.values = values
        self.figs = []
        self.axs = []
        for i in range(5):
            fig, ax = pylab.subplots(3,3)
            self.figs.append(fig)
            self.axs.append(ax.flatten())
            fig.canvas.mpl_connect('pick_event', self.onpick)

        self.plot()

    def onpick(self, event):
        line = event.artist
        xdata, ydata = line.get_data()
        for ind in event.ind:
            print(('Data x[{}]={} y[{}]={} label={}'.format(ind, xdata[ind], ind, ydata[ind], line.get_label())))


    def plot(self):
        values = self.values
        d = self.d
        print(d.hdr)
        print(d.dada_files)
        print(d.dada_files[0].nblocks)
        print(d.nblocks)
        print(d.antennas)
        print(d.nbeams)
        fig, fig2, fig3,fig4,fig5 = self.figs
        ax, ax2, ax3,ax4,ax5 = self.axs

        if values.blkidx >= 0:
            blkidx = values.blkidx
        else:
            blkidx = d.nblocks + values.blkidx

        if values.freq:
            ichan = np.argmin(abs(d.freqs - values.freq))
        else:
            ichan = values.ichan
        ibeam = values.ibeam
        bd = d[blkidx]
        prevbd = None
        if values.ratio:
            prevbd = d[blkidx-1]
            
        iax = 0

        chan_labels = ['chan={}={} MHz'.format(c, freq) for (c, freq) in enumerate(d.freqs)]
        nbeams = d.nbeams_per_antenna

        ant_beam_labels = []
        ant_labels = d.antennas
        antnumbers = list(map(int, [a[2:] for a in d.antennas]))
        if values.ant is None:
            iant = 0
        else:
            iant = antnumbers.index(values.ant)

        # if merging antennas we get way more labels
        if values.antmerge:
            beam_labels = []
            iant = 0
            for a in d.antennas:
                beam_labels.extend(['{} b{}'.format(a, b) for b in range(nbeams)])
            ant_labels = beam_labels
        else:
            beam_labels = ['beam={}'.format(b) for b in range(nbeams)]
                

        xbeams = np.arange(nbeams)
        xant = np.arange(len(d.antennas))
        #plotnames = ['mean','std','kurt','scale','offset', 'decay_offset', 'nsamps', 'gtest']
        plotnames = ['mean','std','kurt','scale','offset', 'decay_offset', 'gtest']
        
        for iname, name in enumerate(plotnames):

            if name == 'gtest':
                tsamp = float(values.nsamps_per_int)
                bdn = bd['mean']**2/bd['std']**2 / tsamp - 1.0
                if prevbd:
                    prevbd['gtest'] = prevbd['mean']**2/prevbd['std']**2 / tsamp - 1.0
            else:
                bdn = bd[name]
                
            if values.log:
                bdn= 10*np.log10(bdn)
            elif values.lognz:
                bdn = 10*np.log10(bdn - bdn.min() + 2) # scalethe log so you dont log0
            if prevbd:
                bdn /= prevbd[name]
                bdn -= 1.

            if values.normant is not None:
                normant = antnumbers.index(values.normant)
                bdn /= bdn[normant, :, :]

            if values.normchan is not None:
                bdn /= bdn[:, :, values.normchan][:,:,np.newaxis]
                
            nant, nbeam, nchan = bdn.shape
            if values.antmerge:
                bdn.shape = (1, nant*nbeam, -1)

            print(name, bdn.shape) # shape is (ant, beam, chan)

            lines1 = showplot(ax[iax],d.freqs, bdn[iant, :, :].T, values, name, beam_labels)
            lines2 = showplot(ax2[iax], xbeams, bdn[:,:,ichan].T, values, name, ant_labels)
            lines3 = showplot(ax3[iax], d.freqs, bdn[:, ibeam, :].T, values, name, ant_labels)
            lines4 = showplot(ax4[iax], xbeams, bdn.mean(axis=2).T, values, name, ant_labels)
            lines5 = showplot(ax5[iax], xbeams, bdn.std(axis=2).T, values, name, ant_labels)
            
            iax += 1

        for iname, name in enumerate(['dm0','dm0count']):
            bdn = bd[name]
            (nant, nbeams, nsamp) = bdn.shape
            xsamp = np.arange(nsamp)
            showplot(ax[iax], xsamp, bdn[iant, :, :].T, values, name, beam_labels)
            showplot(ax2[iax], xbeams, bdn[:,:,values.sample].T , values, name, ant_labels)
            showplot(ax3[iax], xsamp, bdn[:,ibeam,:].T, values, name, ant_labels)
            iax += 1

        iax -= 1
        fig.suptitle('Antenna = {}={}'.format(iant, d.antennas[iant]))
        fig2.suptitle('Channel = {}={} MHz'.format(ichan, d.freqs[ichan]))
        fig3.suptitle('Beam = {}'.format(ibeam))
        fig4.suptitle('Mean over spectrum vs ant/beam')
        fig5.suptitle('Std over spectrum vs ant/beam')
        #fig2.legend(lines2, ant_labels)

        ax[6].set_xlabel('Channel')
        ax2[6].set_xlabel('Beam')
        ax3[6].set_xlabel('Channel')
        ax4[6].set_xlabel('Beam')
        ax5[6].set_xlabel('Beam')
        

        if values.image:
            ax[6].set_ylabel('beam')
            ax2[6].set_ylabel('Antenna')
            ax3[6].set_ylabel('Antenna')

    def show(self):
        pylab.show()

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    #parser.add_argument('-b','--beam', type=int, antenna='bea
    parser.add_argument('-a','--ant', type=int, help='Antenna number', default=None)
    parser.add_argument('-c','--ichan', type=int, help='Channel number', default=0)
    parser.add_argument('-f','--freq', type=float, help='Frequency', default=None)
    parser.add_argument('-t','--sample', type=int, help='Sample number', default=0)
    parser.add_argument('-b','--ibeam', type=int, help='Beam number', default=0)
    parser.add_argument('-i','--blkidx', type=int, help='Block index', default=0)
#    parser.add_argument('-t','--time', type=Stime, help='time to show (UT, MJD, or sample number')
    parser.add_argument('-l','--log', action='store_true', default=False, help='do log on imshow')
    parser.add_argument('-z','--lognz', action='store_true', default=False, help='do non-zero log before plotting')
    parser.add_argument('--image', action='store_true', default=False, help='Show images rathe rthan lines')
    parser.add_argument('-r','--ratio', action='store_true', default=False, help='Plot ratio between this integration and the previous one')
    parser.add_argument('--antmerge', action='store_true', help='Merge antennas with beams so you see everything')
    parser.add_argument('--normant', type=int, help='Normalise by this antenna number')
    parser.add_argument('--normchan', type=int, help='Normalise by this channel')
    parser.add_argument('-I','--nsamps-per-int', type=int, help='Number of samples per integration - for GTEST', default=1)
    parser.add_argument('--exclude-ants', type=strrange, help='Names of antennas to exclude')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        p = RescalePlot(f, values)
        p.show()
    

if __name__ == '__main__':
    _main()
