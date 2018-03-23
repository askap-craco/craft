#!/usr/bin/env python
"""
Correlate vcraft files

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import vcraft
from calc11 import ResultsFile

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def print_delay(xx):
    xxang = np.angle(xx)
    punwrap = np.unwrap(xxang)
    f = np.arange(len(punwrap))
    gradient, phase = np.polyfit(f, punwrap, 1)
    delay = gradient/2./np.pi*len(punwrap)
    print 'Unwrapped phase = {} rad = {} deg, gradient={} rad per channel, delay={} samples nsamp={}' \
        .format(phase, np.degrees(phase), gradient, delay, len(punwrap))

class FitsOut(object):
    def __init__(self, fname, corr):
        ants = self.corr.ants

    def put_product(self, corr, a1, a2, c, p1, p2, xx):
        pass

class PlotOut(object):
    def __init__(self, corr):
        self.stuff = []


    def put_product(self, a1, a2, c, p1, p2, xx):
        print a1.antname, a2.antname, c, p1, p2, xx.shape
        if a1 != a2:
            self.stuff.append(xx)
            fig, (ax1, ax2, ax3, ax4) = pylab.subplots(4,1)
            xxang = np.angle(xx)
            print_delay(xx)
            ax1.plot(abs(xx))
            ax2.plot(np.degrees(xxang))
            for i in xrange(8):
                ax2.axvline(54*i)
                print_delay(xx[54*i:54*(i+1)])

            ax3.plot(np.fft.fftshift(abs(np.fft.fft(xx))), label='lag')
            ax4.imshow(np.angle(np.array(self.stuff)), aspect='auto')
            pylab.show()

    def finish(self):
        stuff = np.array(self.stuff)
        print 'STUFF', stuff
        pylab.imshow(np.angle(stuff))
        pylab.show()


class AntennaSource(object):
    def __init__(self, vfile):
        self.vfile = vcraft.VcraftFile(vfile)
        self.antname = self.vfile.hdr['ANT'][0].lower()
        self.antno = int(self.vfile.hdr['ANTENNA_NO'][0])
        self.mjdstart = float(self.vfile.hdr['TRIGGER_MJD'][0])
        self.trigger_frame = int(self.vfile.hdr['TRIGGER_FRAMEID'][0])
        self.hdr = self.vfile.hdr

    def do_f(self, corr):
        self.frparams = FringeRotParams(corr, self)
        # calculate sample start
        framediff = corr.refant.trigger_frame - self.trigger_frame
        geom_delay_samp = self.frparams.delay*corr.fs/1e6
        fixed_delay_samp = corr.get_fixed_delay_usec(self.antno)*corr.fs/1e6
        total_delay_samp = -geom_delay_samp + framediff + fixed_delay_samp
        whole_delay = int(np.round(total_delay_samp))
        frac_delay = (total_delay_samp - whole_delay)

        # get data
        nsamp = corr.nint*corr.nfft
        logging.debug('F %s sample delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d',
            self.antname, framediff, geom_delay_samp, fixed_delay_samp, total_delay_samp, whole_delay, frac_delay, nsamp)
        rawd = self.vfile.read(corr.curr_samp + whole_delay, nsamp)
        assert rawd.shape == (nsamp, corr.ncoarse_chan)
        self.data = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        d1 = self.data
        geo_delay_rate = self.frparams.delay_rate
        nfine = corr.nfft - 2*corr.nguard_chan
        for c in xrange(corr.ncoarse_chan):
            cfreq = self.vfile.freqs[c]
            cbw = corr.coarse_chanbw/2.
            freqs = np.linspace( - cbw, +cbw, nfine) + (cfreq - corr.f0)
            x1 = rawd[:, c].reshape(-1, corr.nfft)
            xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)
            xfguard = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine]
            phases = np.zeros((corr.nint, nfine))

            for i in xrange(corr.nint):
                delta_t = (frac_delay - i*geo_delay_rate)
                theta0 = 0
                phases[i, :] = 2*np.pi*freqs*delta_t + theta0

            phasor = np.exp(1j*phases)
            xfguard *= phasor
            # slice out only useful channels
            fcstart = c*nfine
            fcend = (c+1)*nfine
            self.data[:, fcstart:fcend, 0] = xfguard


    def get_data(self, chan, pol):
        return self.data[:, chan, pol]


class FringeRotParams(object):
    cols = ('U (m)', 'V (m)', 'W (m)', 'DELAY (us)')

    def __init__(self, corr, ant):
        mid_data = corr.frdata_mid[ant.antname]
        self.u,self.v,self.w,self.delay = [mid_data[c] for c in FringeRotParams.cols]
        self.delay_start = corr.frdata_start[ant.antname]['DELAY (us)']
        self.delay_end = corr.frdata_end[ant.antname]['DELAY (us)']
        self.delay_rate = (self.delay_end - self.delay_start)/float(corr.nint)
        self.ant = ant
        self.corr = corr


class Correlator(object):
    def __init__(self, ants, values):
        self.ants = ants
        self.refant = ants[0]
        self.calcresults = ResultsFile(values.calcfile)
        self.mjd0 = self.refant.mjdstart
        self.frame0 = self.refant.trigger_frame
        self.values = values
        self.nint = 64*64*8
        self.nfft = 64
        self.nguard_chan = 5
        self.fs = 1e6*32./27. # samples per second
        self.ncoarse_chan = 8
        self.coarse_chanbw = 1.0
        self.nfine_chan = self.ncoarse_chan*(64 - 2*self.nguard_chan)
        self.npol_in = 1
        self.f0 = self.ants[0].vfile.freqs.mean() # centre frequency for fringe rotation
        self.inttime_secs = self.nint*self.nfft/self.fs
        self.inttime_days = self.inttime_secs/86400.
        self.curr_intno = 10
        self.curr_samp = self.curr_intno*self.nint
        self.prodout = PlotOut(self)
        self.ant_delays = {}
        self.calcmjd()
        self.get_fr_data()
        self.parse_parset()

    def parse_parset(self):
        with open(self.values.parset, 'rU') as f:
            for line in f:
                if '=' not in line:
                    continue

                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                namebits = name.split('.')
                if line.startswith('common.antenna.ant') and namebits[3] == 'delay':
                    antno = int(namebits[2][3:])
                    delayns = float(value.replace('ns',''))
                    delayus = delayns/1e3
                    self.ant_delays[antno] = delayus

    def get_fixed_delay_usec(self, antno):
        return self.ant_delays[antno]

    def calcmjd(self):
        i = float(self.curr_intno)
        self.curr_mjd_start = self.mjd0 + self.inttime_days*(i + 0.0)
        self.curr_mjd_mid = self.mjd0 + self.inttime_days*(i + 0.5)
        self.curr_mjd_end = self.mjd0 + self.inttime_days*(i + 1.0)

    def next_integration(self):
        #self.curr_intno +=
        self.curr_intno += 1
        self.curr_samp += self.nint
        self.calcmjd()
        self.get_fr_data()

    def get_calc_results(self, mjd):
        res = self.calcresults.scans[0].eval_src0_poly_delta(mjd, self.refant.antname.lower())

        return res

    def get_fr_data(self):
        self.frdata_start = self.get_calc_results(self.curr_mjd_start)
        self.frdata_mid = self.get_calc_results(self.curr_mjd_mid)
        self.frdata_end = self.get_calc_results(self.curr_mjd_end)


    def do_f(self):
        for iant, ant in enumerate(self.ants):
            ant.do_f(self)

    def do_x(self):
        nant = len(self.ants)
        for ia1 in xrange(nant):
            for ia2 in xrange(ia1, nant):
                a1 = self.ants[ia1]
                a2 = self.ants[ia2]
                self.do_x_corr(a1, a2)

    def do_x_corr(self, a1, a2):
        for p1 in xrange(self.npol_in):
            for p2 in xrange(self.npol_in):
                d1 = a1.data[:, :, p1]
                d2 = a2.data[:, :, p2]
                xx = (d1 * np.conj(d2)).mean(axis=0)
                self.put_product(a1, a2, 0, p1, p2, xx)

    def put_product(self, a1, a2, c, p1, p2, xx):
        self.prodout.put_product(a1, a2, c, p1, p2, xx)


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--offset', type=int, help='Offset samples')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot', default=0)
    parser.add_argument('-n','--fft-size', type=int, help='FFT size per coarse channel', default=128)
    parser.add_argument('--calcfile', help='Calc file for fringe rotation')
    parser.add_argument('-p','--parset', help='Parset for delays')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    calcresults = ResultsFile(values.calcfile)
    a1 = AntennaSource(values.files[0])
    a2 = AntennaSource(values.files[1])

    corr = Correlator([a1,a2], values)
    try:
        while(True):
            corr.do_f()
            corr.do_x()
            corr.next_integration()
    except:
        #corr.prodout.finish()
        raise



    f1 = a1.vfile
    f2 = a2.vfile
    refant = a1
    mjdstart = refant.mjdstart
    frdata0 = calcresults.scans[0].eval_src0_poly_delta(mjdstart, refant.antname.lower())
    a1.update_frparams(frdata0)
    a2.update_frparams(frdata0)
    geom_delay = a2.fr_params['DELAY (us)']
    frdata0 = calcresults.scans[0].eval_src0_poly_delta(mjdstart+14./86400., refant.antname.lower())
    a1.update_frparams(frdata0)
    a2.update_frparams(frdata0)
    geom_delay2 = a2.fr_params['DELAY (us)']
    ddif_us = geom_delay2 - geom_delay
    phasediff = np.degrees(2*np.pi*f1.freqs[0]*1e6*ddif_us/1e6)

    print('mjdstart {} Geoemtric delay is {} at start {} at end diff={} = {} deg'.format(mjdstart, geom_delay, geom_delay2, geom_delay2 - geom_delay, phasediff))



    d1 = f1.read()
    d2 = f2.read()

    # truncate to identical nsamp
    nsamp = min(d1.shape[0], d2.shape[0])
    print 'SHAPE BEFORE', d1.shape, d2.shape, nsamp

    d1 = d1[:nsamp, :]
    d2 = d2[:nsamp, :]
    print 'SHAPE AFTER', d1.shape, d2.shape

    assert d1.shape == d2.shape
    nsamp, nchan = d1.shape
    fig, axes = pylab.subplots(5,1)
    d1ax, d2ax,lagax,pax,lax = axes.flatten()
    N = 4096
    Nc = N*512
    if values.offset is None:
        h = 'TRIGGER_FRAMEID'
        offset = int(f2.hdr[h][0]) - int(f1.hdr[h][0])
    else:
        offset = values.offset

    geom_delay = a2.fr_params['DELAY (us)']
    print 'trigger OFFSET IS {} samples. Geoemtric delay is {} us'.format(offset, geom_delay)

    c = values.channel
    d1ax.plot(d1[:N, c].real, label='real')
    d1ax.plot(d1[:N, c].imag, label='imag')
    d1ax.legend()
    d2ax.plot(d2[:N, c].real)
    d2ax.plot(d2[:N, c].imag)

    Nf = values.fft_size
    shortsamp = ((nsamp-offset)/Nf)*Nf

    assert f1.freqs[c] == f2.freqs[c]
    x1 = d1[offset:shortsamp+offset, c].reshape(-1, Nf)
    xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)

    x2 = d2[:shortsamp, c].reshape(-1, Nf)
    xf2 = np.fft.fftshift(np.fft.fft(x2, axis=1), axes=1)
    xx12 = xf1 * np.conj(xf2)
    xx11 = xf1 * np.conj(xf1)
    xx22 = xf2 * np.conj(xf2)

    print 'PRODUCT SIZE', xx12.shape, xx12.shape[0]*Nf, nsamp
    punwrap = np.unwrap(np.angle(xx12.mean(axis=0)))
    xx = np.arange(len(punwrap))
    gradient, phase = np.polyfit(xx, punwrap, 1)
    delay = 32./27.*gradient/2./np.pi*len(punwrap)
    print 'Unwrapped phase = {} rad, graidnet={} rad per channel, delay={} us'.format(phase, gradient, delay)

    lagax.plot(abs(xx11.mean(axis=0)), label='auto0')
    lagax.plot(abs(xx22.mean(axis=0)), label='auto1')
    lagax.plot(abs(xx12.mean(axis=0)), label='crossamp')
    lagax.legend(frameon=False)
    pax.plot(np.degrees(np.angle(xx12.mean(axis=0))), 'o')
    pax.set_ylabel('Cross phase (deg)')
    pax.set_xlabel('Channel')


    lax.plot(np.fft.fftshift(abs(np.fft.fft(xx12.mean(axis=0)))), label='lag')
    Navg = 128*8
    Nout = xx12.shape[0]/Navg
    xx12 = xx12[:Nout*Navg, :]
    xx12.shape = [Nout, Navg, -1 ]
    xx12a = xx12.mean(axis=1)

    pylab.figure()
    pylab.imshow(np.angle(xx12a), aspect='auto')

    pylab.figure()
    lagvt = np.fft.fftshift(abs(np.fft.fft(xx12a, axis=1)))
    pylab.imshow(lagvt)

    pylab.figure()
    pylab.plot(lagvt.T)

    pylab.figure()
    phasevst = np.angle(xx12a[:, xx12a.shape[1]/2])
    pylab.plot(np.degrees(phasevst))

    pylab.show()


if __name__ == '__main__':
    _main()
