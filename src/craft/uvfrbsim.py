#!/usr/bin/env python
"""
Create simulated FRB in visibility data

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from scipy import constants
from .cmdline import strrange
from astropy.coordinates import SkyCoord
from . import fdmt
from . import simfrb
import subprocess
import shutil
import time
from astropy.io import fits

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def runmir(task, **kwargs):
    cmd = [str(task)]
    argarr = [str(k)+'='+str(v) for k,v, in kwargs.items()]
    cmd.extend(argarr)
    logging.info('Running command: %s', ' '.join(cmd))
    subprocess.check_call(cmd, shell=False)

def uvgen(values):
    sourcefile = values.outfile+'.frbpos'
    with open(sourcefile, 'w+') as fout:
        # Assumign FRB relpos in 2 numbers offset in arcsec
        dra, ddec = list(map(float, values.frb_relpos.split(',')))
        fout.write('{} {} {}\n'.format(values.frb_amp, dra, ddec))
        
    # make source file
    fch1 = values.fch1 # GHz
    foff = values.foff # GHz
    nchan = values.nchan
    freqs = fch1 + np.arange(nchan)*foff
    centerfreq_ghz = freqs.mean()
    width_mhz = nchan*foff*1e3
    inttime=values.tint
    totalbw_ghz = width_mhz/1e3
    duration_sec = values.duration*values.tint/1e3
    duration_hr = duration_sec/3600.
    ha_start=0
    ha_end = ha_start + duration_hr
    outfile = values.outfile+'.mir'
    if os.path.exists(outfile):
        if values.clobber:
            shutil.rmtree(outfile)
        else:
            logging.info('File already exists')
            return
        
    args = {'source':sourcefile,
            'ant':values.antfile,
            'telescop':'askap',
            'baseunit':'3.33564',
            'corr':'{nchan},1,0,{width_mhz}'.format(**locals()),
            'time':'20JUN10:00:00:00',
            'freq':'{centerfreq_ghz:0.3f},{totalbw_ghz:0.3f}'.format(**locals()),
            'harange':'{ha_start:0.6f},{ha_end:0.6f}'.format(**locals()),
            'radec':values.phase_center,
            'inttime':inttime/1e3,
            'lat':-26.6790,
            'out':outfile,
            'stokes':values.stokes}
    
    return runmir('uvgen', **args)

def export_fits(values):
    outfile = values.outfile
    if os.path.exists(outfile):
        if values.clobber:
            os.remove(outfile)
        else:
            logging.info('File already exists')
            return

    args = {'in': values.outfile+'.mir',
            'out': outfile,
            'op': 'uvout'
            }

    return runmir('fits', **args)


def antno(antstr):
    return int(antstr[2:])


def pointsim(amp, lm, p, telnames, freqs, noiseamp=0):
    '''
    Simulates a point source in visibilities

    :amp: Amplitude of the point source
    :lm: (l, m) tuple of the point source position in radians
    :p: calc11 polynomial evaluated at a given mjd
    :telnames: list of telescope names to evaluate the visibility at
    :freqs: np array of frequencies (GHz)
    :noiseamp: Additinal gaussian noise amplitude to generate per channnel
    
    :returns: dictionary keyed by (antenna1, antenna2) with the value
    of the visibilities for all frequencies
    '''
    nant = len(telnames)
    nbl = nant*nant-1//2
    nf = len(freqs)
    lambdas = constants.c / (freqs*1e9)
    l, m = lm
    d = {}
    for ia1, a1 in enumerate(telnames):
        for ia2, a2 in enumerate(telnames[ia1+1:]):
            # UVW im meters
            u_m = p[a1]['U (m)'] - p[a2]['U (m)']
            v_m = p[a1]['V (m)'] - p[a2]['V (m)']
            w_m = p[a1]['W (m)'] - p[a2]['W (m)']
            
            # uvw in wavelengths - these are vectors becuase lambdas is a vector
            u = u_m/lambdas
            v = v_m/lambdas
            w = w_m/lambdas
            
            # TMS equation 3.7 - don't include 1/sqrt(1 - l*l - m*m) denomiator term for point sources
            # so says Dan Mitchell anyway.
            vis = amp*np.exp(-2j*np.pi*(u*l + v*m + w*(np.sqrt(1.0 - l*l - m*m) - 1)))
            if noiseamp > 0:
                vishape = vis.shape
                noise = noiseamp*(np.random.randn(*vishape) + 1j*np.random.randn(*vishape))
                vis += noise
            
            d[(a1, a2)] = vis
            
    return d

def modify_data(amps, h0, nant, values):
    t = 0
    currdate = h0.data[0]['DATE']
    nrows = h0.data.shape[0]
    logging.info('Data shape is %s', h0.data[0]['DATA'].shape)
    ntimes, _ = amps.shape

    f1 = values.fch1
    Nchan = values.nchan
    Npol = 2 # ASKAP has XX and YY only
    chanbw = values.foff

    if np.isinf(values.frb_sn):
        noiseamp = 0
    else:
        nbl = nant*(nant-1)/2
        # S/N = A/noisermsinonechannel/sqrt(Nbl,Npol,Nchan,Ntime)
        noiseamp = values.frb_amp/values.frb_sn/np.sqrt(Nchan*nbl*Npol)

    firstbl = h0.data[0]['BASELINE']
    for irow in range(nrows):
        row = h0.data[irow]
        blid = row['BASELINE']
        if row['DATE'] != currdate or (blid == firstbl and irow != 0):
            diff = row['DATE'] - currdate
            currdate = row['DATE']
            t += 1
            logging.info('Time change by %f milliseconds or new baseline irow=%s t=%d', diff*86400*1e3, irow, t)

        # modify amplitude of real and imaginary parts
        if t >= ntimes:
            amp = 0
        else:
            amp = amps[t, :]


        if noiseamp == 0:
            noise_r = 0
            noise_i = 0
        else:
            noise_r = np.random.randn(Nchan)*noiseamp
            noise_i = np.random.randn(Nchan)*noiseamp

        dreal = row['DATA'][0,0,:,0,0]
        dimag = row['DATA'][0,0,:,0,1]

        if False:
            logging.debug('Dreal %s, amp=%s, noiseamp=%s', dreal, amp, noiseamp)
            pylab.plot(dreal)
            pylab.plot(dimag)
            pylab.plot(dreal*amp)
            pylab.plot(dimag*amp)
            pylab.ylim(-0.1, +1.1)
            pylab.show()

        # real part
        row['DATA'][0,0,:,0,0] = dreal*amp + noise_r

        # imag part
        row['DATA'][0,0,:,0,1] = dimag*amp + noise_i


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Generate simulated FRB in UV data', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--fch1', type=float, help='Center frequency of first channel (GHz)', default=1.0)
    parser.add_argument('--foff', type=float, help='Bandwidth of a channel (GHz)', default=1e-3)
    parser.add_argument('--nchan', type=int, help='Number of channels', default=256)
#    parser.add_argument('-f', '--calcfile', required=True, help='Calc11 im file for UV coordinates and stuff')
    parser.add_argument('-o','--outfile', help='Output data file', default='frb')
    parser.add_argument('--antfile', help='Antenna file', default='askap-ak1-ak30.ant')
    parser.add_argument('--clobber', action='store_true', help='Clobber output files')
    #parser.add_argument('--ignore-ant', help='Antennas numbers to ignore', default='31-36', type=strrange)
    #parser.add_argument('--include-ant', help='Antenna numbers to include (overrides --ignore-ant)', type=strrange)
    #parser.add_argument('--format', help='Output data format. [uvfits, raw]', default='uvfits')
    parser.add_argument('--tint', type=float, help='Integration time (ms)', default=0.876)
    parser.add_argument('--duration', type=float, help='Duration (integraions)', default=256)
    parser.add_argument('--phase_center', help='Phase center of observations. hh:mm:ss,dd:mm:ss format, or as decimal hours and decimal degrees', default='0,-30')
    parser.add_argument('--stokes', help='Sotkes to produce. I or XX,YY', default='I')
    
    parser.add_argument('--frb_tstart', type=float, help='FRB start time (ms) from 0', default=0)
    parser.add_argument('--frb_dm', type=float, help='FRB DM (pc/cm3)', default=0)
    parser.add_argument('--frb_idm', type=float, help='FRB DM (samples) - overrides --frb_dm')
    parser.add_argument('--frb_amp', type=float, help='FRB amplitude', default=1)
    parser.add_argument('--frb_sn', type=float, help='FRB S/N in the image', default=np.inf)
    parser.add_argument('--frb_relpos', help='FRB relative position - dra,dec in arcseconds', default='0,0')
    parser.add_argument('--sim-method', help='Simulation method. see functions in simfrb.py', choices=('mkfrb','mkfrb2','mkfrb_fdmt'), default='mkfrb_fdmt')
    parser.add_argument('--show', help='Show plots', action='store_true')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Simulating FRB. values=%s', values)

    fch1 = values.fch1
    Nchan = values.nchan
    Npol = 2 # ASKAP has XX and YY only
    chanbw = values.foff

    #freqs = f1 + np.arange(Nchan)*chanbw
    #lambdas = constants.c / (freqs*1e9)
    #rfile = 'SB6637neweop_alltell.im'
    #rfile = values.calcfile
    #cfile = calc11.ResultsFile(rfile)
    #mjdstart = cfile.scans[0].first_mjd
    #if values.include_ant:
    #    telnames = filter(lambda t: antno(t) in values.include_ant, cfile.telnames)
    #else:
    #    telnames = filter(lambda t: antno(t) not in values.ignore_ant, cfile.telnames)

    #logging.info('Simulating FRB. Start mjd=%0.5f telnames=%s', mjdstart,  ','.join(telnames))

    noiserms = 0 # Need to make independant noise for every baseline/pol - we'll do that shortly.

    fend = values.fch1 + values.foff*(Nchan-1)
    ddm = simfrb.calc_delta_dm(fch1, values.foff, Nchan, values.tint)

    if values.frb_idm:
        idm = values.frb_idm
        dm = simfrb.idm2dm(fch1, values.foff, Nchan, values.tint, idm)
    else:
        dm = values.frb_dm
        idm = simfrb.dm2idm(fch1, values.foff, Nchan, values.tint, dm)

        
    tdelay = 4.15*dm*(fch1**-2 - fend**-2)
    logging.info('DM is %s idm=%s fch1=%s tdelay=%sms', dm, idm, fch1, tdelay)

    simfunc = getattr(simfrb, values.sim_method)
    amps = simfunc(values.fch1, \
            values.foff, \
            Nchan, \
            values.tint, \
            dm, \
            values.frb_amp, \
            values.frb_tstart, \
            noiserms, \
            values.duration)


    if values.show:
        pylab.imshow(amps.T, aspect='auto', origin='lower', interpolation='nearest')
        pylab.xlabel('Time')
        pylab.ylabel('Frequency')
        pylab.show()


    uvgen(values)
    export_fits(values)

    fitsfile = values.outfile
    hdul = fits.open(fitsfile)
    h0 = hdul[0]
    # fix header
    logging.info('Fixing header %s', h0.header['HISTORY'][16])
    h0.header['HISTORY'][16] = ''
    nant = hdul[1].data.shape[0]

    try:
        modify_data(amps, h0, nant, values)
    finally:
        if os.path.exists(fitsfile):
            os.remove(fitsfile)

        logging.info('Writing output to %s', fitsfile)
        hdul.writeto(fitsfile, output_verify='fix+warn')

        
                        
    

if __name__ == '__main__':
    _main()
