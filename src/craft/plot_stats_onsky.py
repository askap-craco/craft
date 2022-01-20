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
from astropy.io import fits
from astropy.coordinates import SkyCoord
from matplotlib.patches import *

import aplpy

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    underlay = values.files[0]
    stats = np.genfromtxt(values.files[1], dtype=None)
    print(stats)
    cmap = plt.cm.get_cmap('viridis')
    nhrs = np.array([s[6]/3600. for s in stats])

    norm = mpl.colors.Normalize(vmin=0, vmax=nhrs.max())

    #fig = aplpy.FITSFigure(underlay, figsize=(5,2))
    #fig.show_grayscale(vmin=3400,vmax=6400, invert=True)
    pos = SkyCoord('11:23:12 -04:58:49', unit='hour,deg')

    hdulist = fits.open(underlay)
    hdr = hdulist[0].header
    print(str(hdr))
    print(hdulist.info())
    d = hdulist[0].data
    fig, ax = plt.subplots()
    rapix = hdr['CRPIX1']
    radelt = hdr['CDELT1']
    raval = hdr['CRVAL1']
    assert hdr['CTYPE1'] == 'RA---CAR'
    decpix,decdelt,decval = hdr['CRPIX2'], hdr['CDELT2'], hdr['CRVAL2']
    assert hdr['CTYPE2'] == 'DEC--CAR'
    print(d.shape)
    ny, nx = d.shape
    left = (0 - rapix)*radelt + raval
    right = (nx - rapix)*radelt + raval
    
    bottom = (0 - decpix)*decdelt + decval
    top = (ny - decpix)*decdelt + decval
    print('top', top, d.shape[0], decpix, decdelt, decval)
    extents = (left, right, bottom, top)

    print(extents)
    #$aspect = 1./abs((right-left)*15/(top - bottom))

    ax.imshow(d, vmin=3400, vmax=6600, cmap='gray_r', origin='bottom', extent=extents)
    xticks = np.linspace(left, right, 9, endpoint=True)
    xlbls = [r'%d$^h$' % (d/15.) for d in xticks]
    xax = ax.get_xaxis()
    xax.set_ticks(xticks)
    xax.set_ticklabels(xlbls)

    yticks = np.linspace(-90, 0, 5, endpoint=True)
    ylbls = [r'%d$^\circ$'%d for d in yticks]
    yax = ax.get_yaxis()
    yax.set_ticks(yticks)
    yax.set_ticklabels(ylbls)

    

    for ifield, field in enumerate(stats):
        print(ifield, field)
        ra = field[1]
        dec = field[2]
        hrs = nhrs[ifield]
        color = cmap(hrs/nhrs.max())
        circ = plt.Circle((ra, dec), 5.5/2., alpha=0.6, facecolor=color, edgecolor=color)
        height = np.sqrt(32.)
        width = abs(height/np.cos(np.radians(dec)))
        print(height, width)
        if dec < -30:
            rot = 0
        else:
            rot = 170.

        rect = plt.Rectangle((ra, dec), width, height, alpha=0.6, facecolor=color, edgecolor='white', angle=rot)
        ell = Ellipse((ra, dec), width, height, alpha=0.6, facecolor=color, edgecolor='white')
        #fig.show_circles([ra],[dec], radius=5.5/2., alpha=0.6, facecolor=color, edgecolor=color)
        #ax.add_artist(circ)
        #ax.add_artist(rect)
        ax.add_artist(ell)
        

    plt.plot(pos.ra.deg, pos.dec.deg, '+', ms=20, color='r', mew=2)
    
    pylab.xlabel('RA (J2000)')
    pylab.ylabel('Dec (J2000)')
    pylab.savefig(values.files[1]+'.pdf')

    fig = plt.figure(figsize=[1,2])
    ax = fig.add_axes([0.05, 0.05, 0.4, 0.9])
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label('Observing time (hrs)')
    fig.savefig(values.files[1]+'.cb.pdf')

    pylab.show()
    
    
    

if __name__ == '__main__':
    _main()
