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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def commasep(s):
    return map(int, s.split(','))

def floatcommasep(s):
    return map(float, s.split(','))

def next_divisor(x, n):
    for t in xrange(x, n/2):
        if n % t == 0:
            return t

MJD_UNIX_EPOCH = 40587.0

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('--step', type=int, help='time step', default=10)
    parser.add_argument('--ntimes', type=int, help='Numerb of samples per block', default=1024)
    parser.add_argument('-O','--outfile', help='Output file for influxdb data. Inhibits live loading')
    parser.add_argument('--show', action='store_true', help='show plots', default=False)
    parser.set_defaults(verbose=False, nxy="1,1")
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    influxout = None
    client = None
    if values.outfile:
        influxout = open(values.outfile, 'w')
    else:
        from influxdb import InfluxDBClient
        client = InfluxDBClient(host='akingest01', database='craft')

    for file in values.files:
        try:
            get_meas(filename, client, influxout)
        except:
            pass

    influxout.close()
        
def get_meas(filename, client, influxout, values):
    fullpath  = os.path.abspath(filename)
    pathbits = fullpath.split('/') # hacky way of getting sbid, scanid and ant
    print len(pathbits), pathbits
    sbid = pathbits[-5]
    scanid = pathbits[-4]
    ant = pathbits[-3]
    startid = pathbits[-2]
    filebits = filename.split('.')
    assert filename.endswith('.fil')
    b = int(filebits[-2])
    st = CraftStatMaker(filename, values)

    while True:
        s, unix_start = st.next_stat()
        nbeams, nstat = s.shape
        # nanoseconds since unix epoch
        unix_nanosec = int(np.round(unix_start * 1e9))
        outs = 'craftstat,sbid={},scanid={},ant={},beam={} '.format(sbid, scanid, ant, b)
        statbits = ['{}={}'.format(statname, stat) for (statname, stat) in zip(st.stat_names, s[b,: ])]
        outs += ','.join(statbits)
        outs += ' {}\n'.format(unix_nanosec)
        field_dict = {}
        for sname, v in zip(st.stat_names, s[b, :]):
            field_dict[sname] = v
            
        body = {'measurement':'craftstat',
                'tags':{'sbid':sbid,'scanid':scanid,'ant':ant,'beam':b},
                'time':unix_nanosec,
                'fields':field_dict
                }

        if influxout is not None:
            influxout.write(outs)

        if client is not None:
            client.write_points([body])

class CraftStatMaker(object):
    def __init__(self, fname, values):
        self.fname = fname
        self.values = values
        self.ntimes = values.ntimes
        self.tstep = values.step
        self.tstart = 0
        self.fft_freqs = np.array([38.5, 50, 100, 200, 300])
        self.spfile = sigproc.SigprocFile(self.fname)
        self.mjdstart = self.spfile.tstart
        self.tsamp = self.spfile.tsamp
        self.stat_names = ['bmean','bstd', 'dm0std']
        for f in self.fft_freqs:
            self.stat_names.append('f{}'.format(f))

    def next_data(self):
        self.tstart += self.tstep + self.ntimes
        self.spfile.seek_sample(self.tstart)
        assert self.spfile.nbits == 8
        d = np.fromfile(self.spfile.fin, dtype=np.uint8, count=self.ntimes*self.spfile.nchans)
        d.shape = (self.ntimes, 1, self.spfile.nchans)

        return  d

    def next_stat(self):

        bi = self.next_data()
        mjdtime = self.mjdstart + self.tstart * self.tsamp/3600./24.
        unix_start = (self.mjdstart - MJD_UNIX_EPOCH)*86400.
        unix_time = unix_start + self.tstart * self.tsamp 
        self.tstart += self.tstep * self.ntimes

        # Normalise to nominal 0 mean, unit stdev - HACK!
        beams -= 128
        beams /= 18

        stat = np.zeros((1, 3 + len(self.fft_freqs)))

        ntimes, nfreq = bi.shape
            # spectra
        beam_mean = bi.mean(axis=0)
        beam_std = bi.std(axis=0)
        bstd = beam_std.std()
# dm0
        dm0 = bi.mean(axis=1)
        dm0std = dm0.std()

        stat[i, 0] = beam_mean.mean()
        stat[i, 1] = bstd
        stat[i, 2] = dm0std

        dm0f = abs(np.fft.rfft(dm0, axis=0))**2
        ntimes, nchans = bi.shape

        idxs = np.rint(self.fft_freqs * float(ntimes) * self.tsamp).astype(int)
        stat[i, 3:] = dm0f[idxs]/dm0std
        fftfreqs  = np.arange(len(dm0f))/float(ntimes)/self.tsamp

            #pylab.plot(fftfreqs, dm0f)
            #for f in self.fft_freqs:
            #    pylab.axvline(f)
            #pylab.show()

        return stat, unix_time

if __name__ == '__main__':
    _main()
