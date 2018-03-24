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
from corruvfits import CorrUvFitsFile

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

CLIGHT=299792458.0

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

    def put_product(self, corr, a1, a2, xx):
        pass

class PlotOut(object):
    def __init__(self, corr):
        self.stuff = []


    def put_product(self, a1, a2, xxp):
        print a1.antname, a2.antname, xxp.shape
        xx= xxp[:,0]
        if a1 != a2 and False:
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

        (geom_delay_samp, geom_delay_rate) = corr.get_geometric_delay_delayrate(self)
        fixed_delay_samp = corr.get_fixed_delay_usec(self.antno)*corr.fs/1e6
        total_delay_samp = -geom_delay_samp + framediff + fixed_delay_samp
        whole_delay = int(np.round(total_delay_samp))
        frac_delay = (total_delay_samp - whole_delay)

        # get data
        nsamp = corr.nint*corr.nfft
        logging.debug('F %s sample delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d',
            self.antname, framediff, geom_delay_samp, fixed_delay_samp, total_delay_samp, whole_delay, frac_delay, nsamp)
        rawd = self.vfile.read(corr.curr_samp*corr.nfft + whole_delay, nsamp)
        assert rawd.shape == (nsamp, corr.ncoarse_chan)
        self.data = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        d1 = self.data
        nfine = corr.nfft - 2*corr.nguard_chan
        print self.vfile.freqs

        for c in xrange(corr.ncoarse_chan):
            cfreq = self.vfile.freqs[c]
            cbw = corr.coarse_chanbw/2.
            coarse_off = cfreq - corr.f0
            #freqs = np.linspace( -cbw, +cbw, nfine) + (cfreq - corr.f0)
            # half channel offset because DC bin is in the center of the FFT
            freqs = (np.arange(nfine, dtype=np.float) - float(nfine)/2. + 0)*corr.fine_chanbw
            x1 = rawd[:, c].reshape(-1, corr.nfft)
            xf1 = np.fft.fftshift(np.fft.fft(x1, axis=1), axes=1)
            xfguard = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine]
            phases = np.zeros((corr.nint, nfine))

            for i in xrange(corr.nint):
                delta_t = (frac_delay - i*geom_delay_rate)
                # oh man - hard to explain. Need to draw a picture
                theta0 = 2*np.pi*coarse_off*float(nfine)*delta_t*corr.fine_chanbw
                phases[i, :] = 2*np.pi*freqs*delta_t/corr.oversamp + theta0

            # If you plot the phases you're about to correct, after adding a artificial
            # 1 sample delay ad tryig to get rid of it with a phase ramp, it becaomes
            # blatetly clear what you should do
            phasor = np.exp(1j*phases)

            '''
            pylab.figure(10)
            pylab.plot(np.angle(phasor[0, :]))
            pylab.plot(np.angle(phasor[-1:, :]))
            '''

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
    def __init__(self, ants, sources, values):
        self.ants = ants
        self.values = values

        self.parse_parset()

        for ia, a in enumerate(self.ants):
            a.ia = ia
            a.antpos = self.get_ant_location(a.antno)

        refantname = self.parset['cp.ingest.tasks.FringeRotationTask.params.refant'].lower()

        self.refant = filter(lambda a:a.antname == refantname, ants)[0]
        self.calcresults = ResultsFile(values.calcfile)
        self.dutc = -37.0
        self.mjd0 = self.refant.mjdstart + self.dutc/86400.0
        self.frame0 = self.refant.trigger_frame
        self.nint = 1000
        self.nfft = 64
        self.nguard_chan = 5
        self.oversamp = 32./27.
        self.fs = 1e6*self.oversamp # samples per second
        self.ncoarse_chan = 8
        self.coarse_chanbw = 1.0
        self.nfine_per_coarse = self.nfft - 2*self.nguard_chan
        self.nfine_chan = self.ncoarse_chan*self.nfine_per_coarse
        self.fine_chanbw = self.coarse_chanbw / float(self.nfine_per_coarse)
        self.npol_in = 1
        self.npol_out = 1
        self.f0 = self.ants[0].vfile.freqs.mean() # centre frequency for fringe rotation
        self.inttime_secs = self.nint*self.nfft/self.fs
        self.inttime_days = self.inttime_secs/86400.
        self.curr_intno = 10
        self.curr_samp = self.curr_intno*self.nint
        self.prodout = PlotOut(self)
        self.calcmjd()
        self.get_fr_data()
        self.fileout = CorrUvFitsFile('test.fits', self.f0, self.fine_chanbw, \
            self.nfine_chan, self.npol_out, sources, ants)

        logging.debug('F0 %f FINE CHANNEL %f kHz num=%d', self.f0, self.fine_chanbw*1e3, self.nfine_chan)


    def parse_parset(self):
        self.parset = {}
        with open(self.values.parset, 'rU') as f:
            for line in f:
                if '=' not in line:
                    continue

                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                self.parset[name] = value


    def get_ant_location(self, antno):
        key = 'common.antenna.ant{}.location.itrf'.format(antno)
        value = self.parset[key]
        location = map(float, value.replace('[','').replace(']','').split(','))
        return location

    def get_fixed_delay_usec(self, antno):
        key = 'common.antenna.ant{}.delay'.format(antno)
        value = self.parset[key]
        delayns =  float(value.replace('ns',''))
        delayus = delayns/1e3

        return delayus


    def get_uvw(self, ant1, ant2):
        fr1 = FringeRotParams(self, ant1)
        fr2 = FringeRotParams(self, ant2)
        uvw = np.array([fr1.u - fr2.u, fr1.v - fr2.v, fr1.w - fr2.w])/CLIGHT

        return uvw

    def get_geometric_delay_delayrate(self, ant):
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.refant)

        delay = fr1.delay - fr2.delay
        delayrate = fr1.delay_rate - fr2.delay_rate

        return (delay*self.fs/1e6, delayrate*self.fs/1e6)


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
        #res = self.calcresults.scans[0].eval_src0_poly_delta(mjd, self.refant.antname.lower())
        res = self.calcresults.scans[0].eval_src0_poly(mjd)

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
        npolout = self.npol_out
        xx = np.empty([self.nfine_chan, npolout], dtype=np.complex64)
        for p1 in xrange(self.npol_in):
            for p2 in xrange(self.npol_in):
                d1 = a1.data[:, :, p1]
                d2 = a2.data[:, :, p2]
                pout = p2 + p1*self.npol_in

                xx[:,pout] = (d1 * np.conj(d2)).mean(axis=0)
        self.put_product(a1, a2, xx)

    def put_product(self, a1, a2, xx):
        self.prodout.put_product(a1, a2, xx)
        uvw = self.get_uvw(a1, a2)
        self.fileout.put_data(uvw, self.curr_mjd_mid, a1.ia, a2.ia,
            self.inttime_secs, xx)

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
    sources = [{'name':'M87','ra':180.,'dec':15.}]

    a1 = AntennaSource(values.files[0])
    a2 = AntennaSource(values.files[1])

    corr = Correlator([a1,a2], sources, values)
    try:
        while(True):
            #fq = corr.fileout.fq_table()
            #an = corr.fileout.an_table(corr.ants)
            #su = corr.fileout.su_table(sources)
            corr.do_f()
            corr.do_x()
            corr.next_integration()
    finally:
        corr.fileout.close()




if __name__ == '__main__':
    _main()
