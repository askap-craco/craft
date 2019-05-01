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
from rtdata import FreddaRescaleData


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def showplot(ax, x, d, values, name, label=None, **kwargs):
    if values.image:
        r = ax.imshow(d.T, aspect='auto')
    else:
        (nrow, ncol) = d.shape
        # need to plot one at a time fo rlabelling. Grrr.
        r = []
        lbl = None
        for n in xrange(ncol):
            if label is not None:
                lbl = label[n]
                
            l = ax.plot(x, d[:, n], picker=3, label=lbl, **kwargs)
            r.append(l)
            
            
    ax.text(0.05, 1, name, ha='left', va='top', transform=ax.transAxes)
    return r

class RescalePlot(object):
    def __init__(self, f, values):
        self.thefile = f
        self.d = FreddaRescaleData(f)
        self.values = values
        self.figs = []
        self.axs = []
        for i in xrange(3):
            fig, ax = pylab.subplots(3,3)
            self.figs.append(fig)
            self.axs.append(ax.flatten())
            fig.canvas.mpl_connect('pick_event', self.onpick)

        self.plot()

    def onpick(self, event):
        line = event.artist
        xdata, ydata = line.get_data()
        for ind in event.ind:
            print('Data x[{}]={} y[{}]={} label={}'.format(ind, xdata[ind], ind, ydata[ind], line.get_label()))


    def plot(self):
        values = self.values
        d = self.d
        print d.hdr
        print d.dada_files
        print d.dada_files[0].nblocks
        print d.nblocks
        print d.antennas
        print d.nbeams
        fig, fig2, fig3 = self.figs
        ax, ax2, ax3 = self.axs
    
        blkidx = values.blkidx
        iant = values.iant
        ichan = values.ichan
        ibeam = values.ibeam
        bd = d[blkidx]
        prevbd = None
        if values.ratio:
            prevbd = d[blkidx-1]
            
        iax = 0
        ant_labels = d.antennas
        chan_labels = ['chan={}={} MHz'.format(c, freq) for (c, freq) in enumerate(d.freqs)]
        nbeams = d.nbeams_per_antenna
        beam_labels = ['beam={}'.format(b) for b in xrange(nbeams)]
        xbeams = np.arange(nbeams)
        xant = np.arange(len(d.antennas))

        for iname, name in enumerate(['mean','std','kurt','scale','offset', 'decay_offset', 'nsamps']):

            bdn = bd[name]
            if values.log:
                bdn= 10*np.log10(bdn)
            elif values.lognz:
                bdn = 10*np.log10(bdn - bdn.min() + 1) # scalethe log so you dont log0
            if prevbd:
                bdn /= prevbd[name]
            
            print name, bdn.shape # shape is (ant, beam, chan)
            lines1 = showplot(ax[iax],d.freqs, bdn[iant, :, :].T, values, name,label=beam_labels)
            lines2 = showplot(ax2[iax], xbeams, bdn[:,:,ichan].T, values, name,label=ant_labels)
            lines3 = showplot(ax3[iax], d.freqs, bdn[:, ibeam, :].T, values, name,label=ant_labels)
            
            iax += 1

        for iname, name in enumerate(['dm0','dm0count']):
            bdn = bd[name]
            (nant, nbeams, nsamp) = bdn.shape
            xsamp = np.arange(nsamp)
            showplot(ax[iax], xsamp, bdn[iant, :, :].T, values, name)
            showplot(ax2[iax], xbeams, bdn[:,:,values.sample].T , values, name)
            showplot(ax3[iax], xsamp, bdn[:,ibeam,:].T, values, name)
            iax += 1

        iax -= 1
        fig.suptitle('Antenna = {}={}'.format(iant, d.antennas[iant]))
        fig2.suptitle('Channel = {}={} MHz'.format(ichan, d.freqs[ichan]))
        fig3.suptitle('Beam = {}'.format(ibeam))
        #fig2.legend(lines2, ant_labels)

        ax[6].set_xlabel('Channel')
        ax2[6].set_xlabel('Beam')
        ax3[6].set_xlabel('Channel')

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
    parser.add_argument('-a','--iant', type=int, help='Antenna number', default=0)
    parser.add_argument('-c','--ichan', type=int, help='Channel number', default=0)
    parser.add_argument('-t','--sample', type=int, help='Sample number', default=0)
    parser.add_argument('-b','--ibeam', type=int, help='Beam number', default=0)
    parser.add_argument('-i','--blkidx', type=int, help='Block index', default=0)
    parser.add_argument('-l','--log', action='store_true', default=False, help='do log on imshow')
    parser.add_argument('-z','--lognz', action='store_true', default=False, help='do non-zero log before plotting')
    parser.add_argument('--image', action='store_true', default=False, help='Show images rathe rthan lines')
    parser.add_argument('-r','--ratio', action='store_true', default=False, help='Plot ratio between this integration and the previous one')
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
