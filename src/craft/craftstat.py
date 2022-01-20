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
from .craftobs import load_beams
from .plotutil import subplots
import matplotlib.gridspec as gridspec
from . import sigproc

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def commasep(s):
    return list(map(int, s.split(',')))

def floatcommasep(s):
    return list(map(float, s.split(',')))

def next_divisor(x, n):
    for t in range(x, n/2):
        if n % t == 0:
            return t

MJD_UNIX_EPOCH = 40587.0

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('--step', type=int, help='time step', default=300)
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
        client = InfluxDBClient(host='akingest01', database='craft', username='craftwriter', password='craft')

    for filename in values.files:
        try:
            print(filename)
            get_meas(filename, client, influxout, values)
        except:
            logging.exception('Blah exception in get_meas')

    if influxout is not None:
        influxout.close()
        
def get_meas(filename, client, influxout, values):
    fullpath  = os.path.abspath(filename)
    pathbits = fullpath.split('/') # hacky way of getting sbid, scanid and ant
    print(len(pathbits), pathbits)
    sbid = pathbits[-5]
    scanid = pathbits[-4]
    ant = pathbits[-3]
    startid = pathbits[-2]
    filebits = filename.split('.')
    assert filename.endswith('.fil')
    b = int(filebits[-2])
    st = CraftStatMaker(filename, values)

    while True:
        try:
            s, unix_start = st.stats()
        except StopIteration:
            break
        # nanoseconds since unix epoch
        unix_nanosec = int(np.round(unix_start * 1e9))
        outs = 'craftstat,sbid={},scanid={},ant={},beam={} '.format(sbid, scanid, ant, b)
        statbits = ['{}={}'.format(statname, stat) for (statname, stat) in s.items()]
        outs += ','.join(statbits)
        outs += ' {}\n'.format(unix_nanosec)
        field_dict = {'sbid':sbid,'scanid':scanid}
        for sname, v in s.items():
            field_dict[sname] = v

        body = {'measurement':'craftstat',
                'tags':{'ant':ant,'beam':b},
                'time':unix_nanosec,
                'fields':field_dict
                }

        if influxout is not None:
            influxout.write(outs)

        if client is not None:
            logging.debug("Writing data to client %s", str(body))
            client.write_points([body])

            
    logging.debug('Finished')

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
        self.stat_names = []
        for f in self.fft_freqs:
            self.stat_names.append('f{}'.format(f))

    def next_data(self):
        self.tstart += self.tstep + self.ntimes
        self.spfile.seek_sample(self.tstart)
        count = self.ntimes
        try:
            d = self.spfile[self.tstart:self.tstart+count]
            return d
        except:
            logging.debug("Finished at %d", self.tstart)
            raise StopIteration
            

    def stats(self):

        bi = self.next_data()
            
        mjdtime = self.mjdstart + self.tstart * self.tsamp/3600./24.
        unix_start = (self.mjdstart - MJD_UNIX_EPOCH)*86400.
        unix_time = unix_start + self.tstart * self.tsamp 
        self.tstart += self.tstep * self.ntimes

        # Normalise to nominal 0 mean, unit stdev - HACK!
        if self.spfile.nbits == 8:
            bi -= 128
            bi /= 18

        stat = np.zeros((3 + len(self.fft_freqs)))

        ntimes, nfreq = bi.shape
            # spectra
        beam_mean = bi.mean(axis=0)
        beam_std = bi.std(axis=0)

        freqs = np.arange(self.spfile.nchans)*self.spfile.foff + self.spfile.fch1
        bstd = beam_std.std()
# dm0
        dm0 = bi.mean(axis=1)
        dm0std = dm0.std()

        stat = {}
        stat['spec.mean.mean'] = beam_mean.mean()
        stat['spec.mean.max'] = beam_mean.max()
        stat['spec.mean.min'] = beam_mean.min()
        stat['spec.mean.maxfreq'] = freqs[beam_mean.argmax()]
        stat['spec.mean.minfreq'] = freqs[beam_mean.argmin()]

        stat['spec.std.mean'] = beam_std.mean()
        stat['spec.std.max'] = beam_std.max()
        stat['spec.std.min'] = beam_std.min()
        stat['spec.std.maxfreq'] = freqs[beam_std.argmax()]
        stat['spec.std.minfreq'] = freqs[beam_std.argmin()]
        stat['dm0.mean'] = dm0.mean()
        stat['dm0.std'] = dm0.std()

        

        dm0f = abs(np.fft.rfft(dm0, axis=0))**2/dm0std**2
        ntimes, nchans = bi.shape

        idxs = np.rint(self.fft_freqs * float(ntimes) * self.tsamp).astype(int)
        for n, i in zip(self.fft_freqs, idxs):
            if i < len(dm0f): # Sample rate might not be high enough to measure all freqs
                freq = '{:0.1f}'.format(n).replace('.','d')
                stat['fft.{}Hz'.format(freq)] = dm0f[i]

        stat['fft.max'] = dm0f[1:].max()
        stat['fft.mean'] = dm0f[1:].mean()
        stat['fft.std'] = dm0f[1:].std()
        fmax = dm0f[1:].max()
        fmax_idx  = np.argmax(dm0f[1:])
        fftfreqs  = np.arange(len(dm0f))/float(ntimes)/self.tsamp
        stat['fft.max.freq'] = fftfreqs[dm0f[1:].argmax() + 1]
            

        if self.values.show:
            f, axes = pylab.subplots(2,2)
            axes = axes.flatten()
            axes[0].imshow(bi.T, aspect='auto')
            axes[0].set_xlabel('Time')
            axes[0].set_xlabel('Channel')
            axes[1].semilogy(fftfreqs, dm0f)
            for f, i in zip(self.fft_freqs, idxs):
                #axes[1].axvline(f)
                if i < len(dm0f):
                    axes[1].semilogy(fftfreqs[i], dm0f[i], 'o')

            axes[1].semilogy(stat['fft.max.freq'], stat['fft.max'], 'x')

            axes[1].set_xlabel('Frequency Hz')
            axes[1].set_ylabel('Amp')
            axes[2].plot(freqs, beam_mean)
            axes[3].plot(freqs, beam_std)
            

            pylab.show()

        return stat, unix_time

if __name__ == '__main__':
    _main()
