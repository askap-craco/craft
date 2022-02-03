#!/usr/bin/env python
"""
Matplotlib plot utilities

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def subplots(*args, **kwargs):
    fig, axes =  pylab.subplots(*args, **kwargs)
    if not hasattr(axes, '__len__'):
        axes = np.array([axes])

    return fig, axes


def onpick(event):
    ''' USE LIKE THIS:
        
    pylab.gcf().canvas.mpl_connect('pick_event', onpick)
    '''
    thisline = event.artist
    xdata, ydata = thisline.get_data()
    ind = event.ind

    print(thisline.get_label(), xdata[ind], ydata[ind])

def addpick(fig=None):
    if fig is None:
        fig = pylab.gcf()

    fig.canvas.mpl_connect('pick_event', onpick)

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
