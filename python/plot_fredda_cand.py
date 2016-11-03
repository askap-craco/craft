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
import plotutil
import fnmatch

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def find_files(rootd, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(rootd):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False, show=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for fin in values.files:
        if os.path.isdir(fin):
            plot_dir(fin, values)
        else:
            ax = pylab.gca()
            plot_file(fin, values, ax, title=fin)

        fout = fin
        if fout.endswith('/'):
            fout = fout[0:-1]
            
        fout = fout + '_fredda_cand.png'
        print 'Saving', fout
        pylab.savefig(fout, dpi=200)
        
        if values.show:
            pylab.show()


def plot_dir(din, values):
    candfiles = find_files(din, 'fredda.cand')
    print 'len candfiles', len(candfiles)
    if len(candfiles) == 0:
        print din, 'is empty'
        return

    ncols = 4
    nrows = len(candfiles)/ncols
    if nrows*ncols < len(candfiles):
        nrows += 1

    assert nrows * ncols >= len(candfiles)
    print 'nxy', nrows, ncols
    fig, axes = plotutil.subplots(nrows, ncols, sharex=True, sharey=True)
    fig.set_size_inches([8,6])
    fig.set_dpi(300)
    axes = axes.flatten()
    
    for ic, f in enumerate(candfiles):
        ax=axes[ic]
        subtitle = f.replace(din, '').replace('C000','').replace('fredda.cand','').replace('/','')
        v = plot_file(f, values, ax, labels=False, subtitle=subtitle)

    fig.suptitle(din)
    fig.text(0.5, 0.05, 'Time(s)', ha='center',va='bottom')
    fig.text(0.05, 0.5, 'Delta_t', rotation=90, ha='center', va='top')
    

def plot_file(fin, values, ax, title=None, labels=True, subtitle=None):
    vin = np.loadtxt(fin)

    if subtitle:
        ax.text(0.05, 0.95, subtitle, ha='left', va='top', transform=ax.transAxes)


    if len(vin) == 0:
        print fin, 'is empty'
        if subtitle:
            s = 'Empty' % subtitle
        else:
            s = '%s is empty' % fin.replace('/', '/\n')
        ax.text(0.5, 0.5, s, ha='center', va='center', fontsize=8, transform=ax.transAxes)
        return vin

    print vin.shape
    sn = vin[:, 0]
    sampno = vin[:, 1]
    time = vin[:, 2]
    boxcar = vin[:, 3]
    dm = vin[:, 4]
    beamno = vin[:, 5]

    ubeams = set(beamno)
    for b in sorted(ubeams):
        bmask = beamno == b
        ax.plot(time[bmask], dm[bmask]+1, marker='x', ls='None', label='Beam %d' %b, picker=3)
               
    plotutil.addpick()
    
    if labels:
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('DM (delta_t)')

    if title:
        ax.set_title(fin)


    return vin

if __name__ == '__main__':
    _main()
