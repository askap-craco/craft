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
from astropy.coordinates import SkyCoord
import multiprocessing
import signal
import warnings

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

CLIGHT=299792458.0

def print_delay(xx):
    xxang = np.angle(xx)
    punwrap = np.unwrap(xxang)
    f = np.arange(len(punwrap)) - len(punwrap)/2
    gradient, phase = np.polyfit(f, punwrap, 1)
    '''
    pylab.figure(20)
    pylab.plot(f, xxang)
    pylab.plot(f, punwrap, 'x')
    pylab.plot(f, np.polyval((gradient, phase), f))
    pylab.show()
    '''
    delay = gradient/2./np.pi*len(punwrap)
    delayns = delay*27./32.*1e3*(54./len(punwrap))
    print 'Unwrapped phase = {} rad = {} deg, gradient={} rad per channel, delay={}samples ={} ns nsamp={}' \
        .format(phase, np.degrees(phase), gradient, delay, delay*27./32.*1e3, len(punwrap))

    return (delay, np.degrees(phase))

class FitsOut(object):
    def __init__(self, fname, corr):
        ants = self.corr.ants

    def put_product(self, corr, a1, a2, xx):
        pass

class PlotOut(object):
    def __init__(self, corr):
        self.stuff = []
        self.corr = corr
        self.delayout = open(corr.values.outfile.replace('.fits','')+'.delays.parset', 'w')

    def put_product(self, a1, a2, xxp):
        if self.corr.values.show:
            if a1 != a2:
                xx= xxp[:,0]
                self.stuff.append(xx)
                self.last_xx = (a1, a2, xx)
                if len(self.stuff) >0:
                    self.plot_stuff(a1, a2)

    def plot_stuff(self, a1, a2):

        (a1, a2, xx) = self.last_xx
        fig, (ax1, ax2, ax3, ax4, ax5,ax6) = pylab.subplots(6,1)
        xxang = np.angle(xx)
        print_delay(xx)
        ax1.plot(abs(xx))
        ax2.plot(np.degrees(xxang), 'o')
        nf = self.corr.nfine_out_per_coarse
        for i in xrange(self.corr.ncoarse_chan):
            ax2.axvline(nf*i, c='r')
            ax1.axvline(nf*i, c='r')
            print 'PLOTSTUFF', nf, i, nf*i, nf*(i+1), self.corr.ncoarse_chan, xx.shape
            print_delay(xx[nf*i:nf*(i+1)])

        lag = np.fft.fftshift(abs(np.fft.fft(xx)))
        lagidx = np.argmax(lag)
        lagoff = lagidx - len(xx)/2.0
        lagns= lagoff/self.corr.full_bw*1e3
        lagsn = lag[lagidx]/np.std(lag[0:lagidx-100])
        lagx = np.arange(len(lag)) - len(lag)/2.0
        ax3.plot(lagx/self.corr.full_bw*1e3, lag,label='lag')
        ax3.axvline(0, c='r')
        ax3.set_xlabel('Delay (ns)')
        angxx = np.angle(np.array(self.stuff))
        ax4.imshow(angxx, aspect='auto')
        phase_vt = np.degrees(angxx[:, 10:50].mean(axis=1))
        t = np.arange(len(phase_vt))
        rate, offset = np.polyfit(t, phase_vt, 1)
        fit_std = (np.polyval((rate, offset), t) - phase_vt).std()
        ax5.plot(phase_vt)
        ax6.imshow(abs(np.array(self.stuff)), aspect='auto')
        print 'Phase rate={} offset={} std={} deg lagoff={}samples = {}ns lag S/N={}'.format(rate, offset, fit_std, lagoff, lagns, lagsn)
        fig.suptitle('a1 {} a2 {}'.format(a1.antname, a2.antname))
        curr_delay = self.corr.get_fixed_delay_usec(a2.antno)*1e3
        if a1 != a2:
            self.delayout.write('#{}={} S/N={} delay={}ns\n'.format(a1.antname, a2.antname, lagsn, lagns))
            self.delayout.write('common.antenna.ant{}.delay={}ns\n'.format(a2.antno, curr_delay+lagns))
            self.delayout.flush()

        pylab.show()

    def finish(self):
        stuff = np.array(self.stuff)
        pylab.imshow(np.angle(stuff))
        pylab.show()


class AntennaSource(object):
    def __init__(self, vfile):
        self.vfile = vfile
        self.antname = self.vfile.hdr['ANT'][0].lower()
        self.antno = int(self.vfile.hdr['ANTENNA_NO'][0])
        self.mjdstart = self.vfile.start_mjd
        self.trigger_frame = self.vfile.start_frameid
        self.hdr = self.vfile.hdr
        self.init_geom_delay_us = None
        self.all_geom_delays = []
        self.all_mjds = []
        print 'antenna {} {}'.format(self.antname, self.vfile.freqconfig)

    def do_f(self, corr):
        self.frparams = FringeRotParams(corr, self)
        # calculate sample start
        framediff_samp = corr.refant.trigger_frame - self.trigger_frame
        framediff_us = framediff_samp / corr.fs
        (geom_delay_us, geom_delay_rate_us) = corr.get_geometric_delay_delayrate_us(self)
        self.all_geom_delays.append(geom_delay_us)
        self.all_mjds.append(corr.curr_mjd_mid)

        #logging.debug'ALL GEOM DELAYS', self.all_geom_delays, type(geom_delay_us)
        #print 'ALL MJDs', self.all_mjds
        # test with 0 delay rate
        #geom_delay_us = self.all_geom_delays[0]
        
        geom_delay_samp = geom_delay_us * corr.fs
        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
        fixed_delay_samp = fixed_delay_us*corr.fs
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))
        total_delay_us = total_delay_samp / corr.fs
        whole_delay_us = whole_delay / corr.fs

        frac_delay_samp = (total_delay_samp - whole_delay)
        frac_delay_us = frac_delay_samp *corr.fs

        # get data
        nsamp = corr.nint*corr.nfft
        logging.debug('F %s sample delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg',
            self.antname, framediff_samp, geom_delay_samp, fixed_delay_samp,
            total_delay_samp, whole_delay, frac_delay_samp, nsamp,
            360.*frac_delay_us*corr.f0)
        logging.debug('F %s us delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg rate=%e',
                self.antname, framediff_us, geom_delay_us, fixed_delay_us,
                total_delay_us, whole_delay_us, frac_delay_us, nsamp,
                      360.*frac_delay_us*corr.f0, geom_delay_rate_us)

        sampoff = corr.curr_samp*corr.nfft + whole_delay
        rawd = self.vfile.read(sampoff, nsamp)
        assert rawd.shape == (nsamp, corr.ncoarse_chan), 'Unexpected shape from vfile: {} expected ({},{})'.format(rawd.shape, nsamp, corr.ncoarse_chan)
        self.data = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        d1 = self.data
        nfine = corr.nfft - 2*corr.nguard_chan

        for c in xrange(corr.ncoarse_chan):
            cfreq = corr.freqs[c]
            freqs = (np.arange(nfine, dtype=np.float) - float(nfine)/2.0)*corr.fine_chanbw
            if corr.sideband == -1:
                freqs = -freqs
                
            x1 = rawd[:, c].reshape(-1, corr.nfft)
            xf1 = np.fft.fft(x1, axis=1)
            xf1 = np.fft.fftshift(xf1, axes=1)
            xfguard = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine:] # scale because oterhwise it overflows
            delta_t = -fixed_delay_us + geom_delay_us
            phases = delta_t * (freqs + cfreq)
            logging.debug('PHASOR %s[%s] chan=%s freq=%sfixed=%f us geom=%f us delta_t %s us coff*fixed = %f deg coff*geom = %f deg',
                          self.antname, self.ia, c, cfreq, fixed_delay_us, geom_delay_us, delta_t, cfreq*fixed_delay_us*360., cfreq*geom_delay_us*360.)


            # If you plot the phases you're about to correct, after adding a artificial
            # 1 sample delay ad tryig to get rid of it with a phase ramp, it becaomes
            # blatetly clear what you should do
            phasor = np.exp(np.pi*2j*phases, dtype=np.complex64)
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
        self.u,self.v,self.w,self.delay = map(float, [mid_data[c] for c in FringeRotParams.cols])
        self.delay_start = float(corr.frdata_start[ant.antname]['DELAY (us)'])
        self.delay_end = float(corr.frdata_end[ant.antname]['DELAY (us)'])
        self.delay_rate = (self.delay_end - self.delay_start)/float(corr.nint)
        self.ant = ant
        self.corr = corr

    def __str__(self):
        s = 'FR {} uvw=({},{},{}) m = {} us'.format(self.ant.antname, self.u, self.v, self.w, self.delay)
        return s

    __repr__ = __str__

class Correlator(object):
    def __init__(self, ants, sources, values):
        self.running = True
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self.ants = ants
        self.values = values
        self.pool = None
        if self.values.num_threads > 1:
            self.pool = multiprocessing.Pool(processes=values.num_threads)

        self.parse_parset()

        for ia, a in enumerate(self.ants):
            a.ia = ia
            a.antpos = self.get_ant_location(a.antno)

        refantname = self.parset['cp.ingest.tasks.FringeRotationTask.params.refant'].lower()

        #self.refant = filter(lambda a:a.antname == refantname, ants)[0]
        self.refant = ants[0]
        self.calcresults = ResultsFile(values.calcfile)
        self.dutc = 0
        self.mjd0 = self.refant.mjdstart + self.dutc/86400.0
        self.frame0 = self.refant.trigger_frame
        self.nint = values.nint
        self.nfft = 64*values.fft_size
        self.nguard_chan = 5*values.fft_size
        self.oversamp = 32./27.
        self.fs = self.oversamp # samples per microsecnd
        self.ncoarse_chan = len(self.refant.vfile.freqs)
        self.sideband = -1
        self.coarse_chanbw = 1.0
        self.nfine_per_coarse = self.nfft - 2*self.nguard_chan
        self.nfine_chan = self.ncoarse_chan*self.nfine_per_coarse
        self.fine_chanbw = self.coarse_chanbw / float(self.nfine_per_coarse)
        self.full_bw = self.fine_chanbw * self.nfine_chan
        self.fscrunch = values.fscrunch
        assert self.fscrunch >= 1
        assert self.nfine_per_coarse % self.fscrunch == 0, 'Fsrunch must yield an integer number of fine channels per coarse channel'
        self.nfine_out_per_coarse = self.nfine_per_coarse / self.fscrunch
        self.nfine_out_chan = self.nfine_out_per_coarse*self.ncoarse_chan
        self.out_chanbw = self.coarse_chanbw / float(self.nfine_out_per_coarse)
        self.npol_in = 1
        self.npol_out = 1
        #self.f0 = self.ants[0].vfile.freqs.mean() # centre frequency for fringe rotation
        self.f0 = self.ants[0].vfile.freqs[0]
        self.freqs = self.ants[0].vfile.freqs
        self.fmid = self.freqs.mean()
        self.inttime_secs = float(self.nint*self.nfft)/(self.fs*1e6)
        self.inttime_days = self.inttime_secs/86400.
        self.curr_intno = 0
        self.curr_samp = self.curr_intno*self.nint + 1000
        self.prodout = PlotOut(self)
        self.calcmjd()
        self.get_fr_data()
        self.fileout = CorrUvFitsFile(values.outfile, self.fmid, self.sideband*self.out_chanbw, \
                                      self.nfine_out_chan, self.npol_out, self.mjd0, sources, ants, self.sideband)

        logging.debug('F0 %f FINE CHANNEL %f kHz num=%d freqs=%s', self.f0, self.fine_chanbw*1e3, self.nfine_chan, self.freqs)

    def exit_gracefully(self, signum, frame):
        self.running = False

    def parse_parset(self):
        self.parset = {}
        with open(self.values.parset, 'rU') as f:
            for line in f:
                if '=' not in line or line.startswith('#'):
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

    def get_geometric_delay_delayrate_us(self, ant):
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.refant)

        delay = fr1.delay - fr2.delay
        delayrate = fr1.delay_rate - fr2.delay_rate

        return (delay, delayrate)


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
            if not self.running:
                raise KeyboardInterrupt()

            ant.do_f(self)

    def do_x(self):
        nant = len(self.ants)
        for ia1 in xrange(nant):
            for ia2 in xrange(ia1, nant):
                if not self.running:
                    raise KeyboardInterrupt()
                    
                a1 = self.ants[ia1]
                a2 = self.ants[ia2]
                self.do_x_corr(a1, a2)

    def do_x_corr(self, a1, a2):
        npolout = self.npol_out
        xx = np.empty([self.nfine_chan, npolout], dtype=np.complex64)
        #np.seterr(all='raise')
        rfidelay = self.values.rfidelay
        for p1 in xrange(self.npol_in):
            for p2 in xrange(self.npol_in):
                d1 = a1.data[:, :, p1]
                d2 = a2.data[:, :, p2]
                pout = p2 + p1*self.npol_in
                ntimesi = d1.shape[0] - rfidelay
                ntimes = float(ntimesi)

                try:
                    for c in xrange(self.nfine_chan):
                        #xx[:,pout] = (d1 * np.conj(d2)).mean(axis=0)
                        # vdot conjugates the first argument
                        # this is equivalent to (d1 * conj(d2)).mean(axis=0)
                        if self.sideband == -1:
                            xx[c, pout] = np.vdot(d2[:ntimesi, c], d1[rfidelay:, c])/ntimes
                        else:# take complex conjugate if inverted
                            xx[c, pout] = np.vdot(d1[:ntimesi, c], d2[rfidelay:, c])/ntimes
                            
                except Exception, e:
                    print 'Error', e
                    import ipdb
                    ipdb.set_trace()


        if self.fscrunch > 1:
            chans = np.arange(self.nfine_chan, step=self.fscrunch)
            xx = np.add.reduceat(xx, chans, axis=0)/float(self.fscrunch)

        assert xx.shape == (self.nfine_out_chan, self.npol_out)
            
        self.put_product(a1, a2, xx)

    def put_product(self, a1, a2, xx):
        self.prodout.put_product(a1, a2, xx)
        uvw = self.get_uvw(a1, a2)
        self.fileout.put_data(uvw, self.curr_mjd_mid, a1.ia, a2.ia,
            self.inttime_secs, xx)


def parse_delays(values):
    delayfile = values.calcfile.replace('.im','.hwdelays')
    delays = {}
    if os.path.exists(delayfile):
        with open(delayfile, 'rU') as dfile:
            for line in dfile:
                bits = line.split()
                if not line.startswith('#') and len(bits) == 2:
                    delays[bits[0].strip()] = -int(bits[1])

        logging.info('Loaded %s delays from %s', len(delays), delayfile)
    else:
        logging.info('No delays loaded. %s does not exist', delayfile)


    return delays

def load_sources(calcfile):
    calc_input = calcfile.replace('.im','.calc')
    d = {}
    for line in open(calc_input, 'rU'):
        if len(line) == 0 or line.startswith('#'):
            continue
        bits = line.split(':')
        if len(bits) != 2:
            continue
        
        k,v = bits

        d[k.strip()] = v.strip()

    assert d['NUM SOURCES'] == '1'
    name = d['SOURCE 0 NAME']
    # ra/dec in radians
    ra = float(d['SOURCE 0 RA'])
    dec = float(d['SOURCE 0 DEC'])
    pos = SkyCoord(ra, dec, unit=('rad','rad'), frame='icrs')
    sources = [{'name':name,'ra':pos.ra.deg,'dec':pos.dec.deg}]

    return sources



def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--outfile', help='Output fits file', default='corr.fits')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot', default=0)
    parser.add_argument('-n','--fft-size', type=int, help='Multiple of 64 channels to make channels- default=1', default=1)
    parser.add_argument('-t','--num-threads', type=int, help='Number of threads to run with', default=1)
    parser.add_argument('--calcfile', help='Calc file for fringe rotation')
    parser.add_argument('-p','--parset', help='Parset for delays')
    parser.add_argument('--show', help='Show plot', action='store_true', default=False)
    parser.add_argument('-i','--nint', help='Number of fine spectra to average', type=int, default=128)
    parser.add_argument('-f','--fscrunch', help='Frequency average by this factor', default=1, type=int)
    parser.add_argument('--rfidelay', type=int, help='Delay in fine samples to add to second component to make an RFI data set', default=0)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    calcresults = ResultsFile(values.calcfile)
    sources = load_sources(values.calcfile)
    # hacking delays
    delaymap = parse_delays(values)
    antennas = [AntennaSource(mux) for mux in vcraft.mux_by_antenna(values.files, delaymap)]
    #antennas = [AntennaSource(vcraft.VcraftFile(f)) for f in values.files]
    corr = Correlator(antennas, sources, values)
    try:
        while(corr.running):
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
