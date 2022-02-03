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
from .rtdata import FreddaRescaleData
from astropy.time import Time
from influxdb import InfluxDBClient


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def convert(v):
    if np.isnan(v):
        return None
    else:
        return v

def isok(v):
    return np.isfinite(v)

class Stats(object):
    def __init__(self, antennas, nbeams, time, values):
        self.antennas = antennas
        self.nbeams = nbeams
        self.data = {}
        self.time = time
        self.values = values

    def to_influx(self):
        s = ''
        statnames = list(self.data.keys())
        body = []
        for ia, a in enumerate(self.antennas):
            for b in range(self.nbeams):
                fields = {statname:self.data[statname][ia, b] for statname in statnames
                          if isok(self.data[statname][ia, b])}
                tags = {'antenna':a, 'beam':b}
                m = {
                    'measurement':'freddars',
                    'tags':tags,
                    'fields':fields,
                    'time':self.time
                    }
                body.append(m)
                
        return body

    def __iadd__(self, v):
        statname, d, summary = v
        assert d.shape == (len(self.antennas), self.nbeams), 'INvalid shape:{}'.format(d.shape)
        statkey = statname + '_' + summary
        self.data[statkey] = d
        if statkey == 'mean_std' and self.values.plot:
            #pylab.imshow(d)
            pylab.plot(d)
            pylab.title(statkey)
            pylab.show()
        
        return self
        

def plot(f, values):
    try :
        d = FreddaRescaleData(f)
    except:
        logging.exception('Could not parse data')
        return
    
    client = InfluxDBClient(database='craft', host='akingest01', username='craftwriter', password='craft')
    print('ANT', d.antennas)
    print('NBEAM', d.nbeams)
    
    for block in d.blocks(step=values.step):
        t = Time(block.mjd, format='mjd')
        print('TIME', t.isot, block.mjd, block.rsdata.tstart, block.rsdata.tsamp)
        stat = Stats(d.antennas, d.nbeams_per_antenna, t.isot, values)
        names = ['mean','std','kurt','scale','offset', 'decay_offset', 'nsamps', 'dm0']
        names = ['mean', 'std','kurt', 'scale', 'offset']
        for iname, name in enumerate(names):
            bdn = block[name] # shape = [nant, nbeam, nchan] nbeam includes both pols
            nzero = (bdn == 0).sum(axis=2)
            ninf = np.isinf(bdn).sum(axis=2)
            nnan = np.isnan(bdn).sum(axis=2)
            stat += name, bdn.max(axis=2), 'max'
            stat += name, bdn.min(axis=2), 'min'
            stat += name, bdn.mean(axis=2), 'mean'
            stat += name, bdn.std(axis=2), 'std'
            stat += name, np.median(bdn, axis=2), 'med'
            stat += name, nzero, 'nzero'
            stat += name, ninf, 'ninf'
            stat += name, nnan, 'nnan'

            # if this axis is the frequency axis - just guessing
            if bdn.shape[2] == len(d.freqs):
                maxfreqidx = np.argmax(bdn, axis=2)
                freqmax = d.freqs[maxfreqidx.flat]
                freqmax.shape = maxfreqidx.shape
                stat += name, freqmax, 'freqmax'
        
        body = stat.to_influx()
        client.write_points(body)
        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--step', type=int, default=16, help='Skip this many blocks per write')
    parser.add_argument('-i','--start', type=int, default=0, help='Start block number')
    parser.add_argument('-p','--plot', action='store_true', help='Do plot', default=False)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        plot(f, values)
    

if __name__ == '__main__':
    _main()
