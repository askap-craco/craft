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
import glob
import ephem
import aces.mro as mro
from . import sigproc
from . import craftobs
from .crafthdr import DadaHeader
import datetime

def opentle(f):
    sats = {}
    for iline, line in enumerate(open(f, 'rU')):
        line = line.strip()
        if iline % 3 == 0:
            name = line
        elif iline % 3 == 1:
            line2 = line
        else:
            line3 = line
            sats[name] = (name, line2, line3)

    return sats

class Trajectory(object):
    def __init__(self, sat, satlines, ntimes, nbeams):
        self.sat = sat
        self.satlines = satlines
        self.times = np.zeros(ntimes)
        self.az = np.zeros(ntimes)
        self.alt = np.zeros(ntimes)
        self.ra = np.zeros(ntimes)
        self.dec = np.zeros(ntimes)
        self.bsepmin = np.zeros(ntimes)
        self.bsepbeam = np.zeros(ntimes)
        self.bsep = np.zeros((ntimes, nbeams))
        self.nt = 0

    def append(self, time, az, alt, ra, dec, bsepmin, bsepbeam, bsep):
        t = self.nt
        self.times[t] = time
        self.az[t] =az
        self.alt[t] = alt
        self.ra[t] = ra
        self.dec[t] = dec
        self.bsepmin[t] = bsepmin
        self.bsepbeam[t] = bsepbeam
        self.bsep[t,:] = bsep
        self.nt += 1

        

def get_trajectories(satellites, times, beam_bodies, septhresh_deg=None):
    trajectories = {}
    for satname, satlines in satellites.items():
        traj = Trajectory(satname, satlines, len(times), len(beam_bodies))
        sat = ephem.readtle(*satlines)

        for it, t in enumerate(times):
            askap = mro.observer()
            askap.date = t
            sat.compute(askap)
            bsep = np.zeros(len(beam_bodies))
            for ibeam, bbody in enumerate(beam_bodies):
                bbody.compute(askap)
                bsep[ibeam] = ephem.separation(bbody, sat)

            sepmin = bsep.min()

            traj.append(t, sat.az, sat.alt, sat.ra, sat.dec, bsep.min(), bsep.argmin(), bsep)

        if septhresh_deg is None or traj.bsep.min() < np.radians(septhresh_deg):
            trajectories[satname] = traj

    return trajectories

def make_bodies_from_arrays(ras, decs):
    bodies = []
    for ra, dec in zip(ras, decs):
        b = ephem.FixedBody()
        b._ra = np.radians(ra)
        b._dec = np.radians(dec)
        bodies.append(b)
        
    return bodies

class Constellation(object):
    def __init__(self, tlename, septhresh, freqs):
        self.tlename = tlename
        self.septhresh = septhresh
        self.freqs = freqs

    def load_tle(self, dir):
        tlefile = os.path.join(dir, self.tlename)
        self.tlefile = tlefile
        self.tlefile_mtime = os.path.getmtime(tlefile)
        self.tle = opentle(tlefile)

    def calc_trajectories(self, times, beam_bodies):
        traj = get_trajectories(self.tle, times, beam_bodies, self.septhresh)
        return traj

    def get_tle_mtime(self):
        return self.tlefile_mtime

def plot(times, beam_bodies, trajectories):

    pylab.figure()
    times_mins = (times - min(times))*24.*60.


    for satname, (const, traj) in trajectories.items():
        ra = np.degrees(traj.ra)
        dec = np.degrees(traj.dec)
        bsepmin = np.degrees(traj.bsepmin)
        pylab.plot(ra, dec, label=satname)

    pylab.legend(frameon=False)

    for ib, b in enumerate(beam_bodies):
        pylab.text(np.degrees(b._ra), np.degrees(b._dec), str(ib), ha='center', va='center')
        pylab.plot(np.degrees(b._ra), np.degrees(b._dec))
    

    xlim = pylab.xlim()
    pylab.xlim(xlim[1], xlim[0])
    pylab.xlabel('Ra (deg)')
    pylab.ylabel('Dec (deg)')

    fig, axes = pylab.subplots(1, 2)
    #fig.set_size_inches([12,4])
    for satname, (const, traj) in trajectories.items():
        bsepmin = np.degrees(traj.bsepmin)
        bsepbeam = traj.bsepbeam
        axes[0].plot(times_mins, bsepmin, label=satname)
        axes[1].plot(times_mins, bsepbeam)
        
    axes[0].legend(frameon=False)
    axes[0].set_xlabel('time (min)')
    axes[0].set_ylabel('Minimum beam separation (deg)')
    axes[1].set_ylabel('Beam number with minimum separation')

def rotate_beams(ra0, dec0, ras, decs, rot_deg=45.0):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    center = SkyCoord(ra0, dec0, unit='deg')
    newframe = center.skyoffset_frame(rotation=rot_deg*u.deg)
    beams = [SkyCoord(r, d, unit='deg') for (r, d) in zip(ras, decs)]
    new_beams = [b.transform_to(newframe) for b in beams]

    new_ras = [b.lon.degree + ra0 for b in new_beams]
    new_decs = [b.lat.degree + dec0 for b in new_beams]

    return new_ras, new_decs

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-d', '--duration', type=float, help='Duration to simulate (seconds)', default=1200.)
    parser.add_argument('-t', '--tdelta', type=float, help='Delta t to simulate (seconds)', default=10.)
    parser.add_argument('--tle-dir', help='Directory to find TLEs in. Defaults to env["TLE_DIR"]', default=os.environ['TLE_DIR'])
    parser.add_argument('-r','--rotate-beams', type=float, help='rotate beams to fix a header bug by a number of degrees')
    parser.add_argument('-s', '--show', action='store_true', default=False, help='Show plots')
    parser.add_argument('-T', '--septhresh', type=float, help='Separation threshold (degrees)', default=5.0)
    parser.add_argument(dest='hdrfile', nargs=1)
    parser.set_defaults(verbose=False, rotate_beams=0.0)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    hdrfile = values.hdrfile[0]
    print('Loading header', hdrfile)
    hdr = DadaHeader.fromfile(hdrfile)
    beamdir = os.path.dirname(os.path.abspath(hdrfile))


    print('Loading beam dir', beamdir)
    beams, filfiles = craftobs.load_beams(beamdir, 128,128, return_files=True)

    if 'TSTART' in hdr:
        mjdstart = float(hdr['TSTART'])
        duration = values.duration/3600./24.
    else:
        mjdstart = filfiles[0].tstart
        duration_sec = filfiles[0].file_size_elements*filfiles[0].tsamp
        duration = duration_sec/3600./24.

    if values.duration is not None:
        duration = values.duration/3600./24.


    ras = list(map(float, hdr['BEAM_RA'][0].split(',')))[0:36]
    decs = list(map(float, hdr['BEAM_DEC'][0].split(',')))[0:36]
    if values.rotate_beams != 0:
        logging.info('Rotating beams by %f', values.rotate_beams)
        ra0 = float(hdr['RA'][0])
        dec0 = float(hdr['DEC'][0])
        ras, decs = rotate_beams(ra0, dec0, ras, decs, values.rotate_beams)

    beam_bodies = make_bodies_from_arrays(ras, decs)
    nbeams = len(beam_bodies)

    # pyephem epoch is 1899 December 31 12:00 according to http://rhodesmill.org/pyephem/tutorial.html
    # I need to subtract 2 days from this to make sense. Wierd
    # i.e. mjd 57681 = 20 Oct 2016
    tdelt = values.tdelta/3600./24.
    mjd0 = 15021.5 - 2.0
    start_time = ephem.Date(mjdstart - mjd0)
    print('start time', start_time, float(start_time))
    times = np.arange(float(start_time), float(start_time)+ duration, tdelt)
    times_mjd = times + mjd0

    # GPS: http://www.navipedia.net/index.php/GPS_Signal_Plan
    gps = Constellation('gps-ops.txt', # tle name
                        values.septhresh, # separation threshold
                        ((1575.42, 30.0, 'L1'), 
                         (1227.60, 30.0, 'L2'), 
                         (1176.45, 30.0, 'L5')))

    #Galileo http://www.navipedia.net/index.php/Galileo_Signal_Plan

    galileo = Constellation('galileo.txt',
                            values.septhresh,
                            ((1575.42, 40.0, 'E1'),
                             (1278.75, 40.0, 'E6'),
                             (1191.795, 90.0, 'E5')))

    # Beidou signal plan: http://www.navipedia.net/index.php/BeiDou_Signal_Plan
    beidou = Constellation('beidou.txt', 
                           values.septhresh,
                           (((1561.098 + 1589.742)/2., 50.0, 'B1'),
                            (1207.14, 30.0, 'B2'),
                            (1268.52, 40.0, 'B3')))

    # Glonass: http://www.navipedia.net/index.php/GLONASS_Signal_Plan
    glonass = Constellation('glo-ops.txt', values.septhresh,
                            (((1598.0625 + 1605.375)/2.0, 18.0, 'L1'),
                             ((1242.9375 + 1248.625)/2.0, 18.0, 'L2'),
                             (1201.0, 20.0, 'L3')))

    constellations = [gps, galileo, beidou, glonass]

    all_traj = {}

    for const in constellations:
        const.load_tle(values.tle_dir)
        trajectories = const.calc_trajectories(times, beam_bodies)
        for satname, traj in trajectories.items():
            all_traj[satname] = (const, traj)

    if values.show:
        plot(times, beam_bodies, all_traj)


    flagfile = os.path.join(beamdir, 'satellites.flags')
    flout = open(flagfile, 'w')
    flout.write('# CRAFT Flagfile written by {}\n'.format(sys.argv[0]))
    flout.write('# Cmdline: {}\n'.format(' '.join(sys.argv)))
    flout.write('# Local Now:{}\n'.format(datetime.datetime.now().isoformat()))
    flout.write('# UTC Now:{}\n'.format(datetime.datetime.utcnow().isoformat()))
    flout.write('# Beam dir:{}\n'.format(beamdir))
    flout.write('# TLE dir: {}\n'.format(values.tle_dir))
    for const in constellations:
        flout.write('# TLEFILE {}. {} satellites. Last modified {}\n'.
                    format(const.tlename, len(const.tle), datetime.datetime.fromtimestamp(const.get_tle_mtime()).isoformat()))

    flout.write('# {} trajectories are within the threshold\n'.format(len(all_traj)))
    for satname, (const, traj) in all_traj.items():
        bsepmin = np.degrees(traj.bsepmin)
        bsepbeam = traj.bsepbeam
        minidx = bsepmin.argmin()
        tlelines = const.tle[satname]
        flout.write('#\n#{} minimum beam separation {:0.2f} deg at {}. Beams {}\n'.format(satname, bsepmin.min(), times_mjd[minidx], sorted(set(bsepbeam))))
        flout.write('#\n#{0}\n#{1}\n#{2}\n'.format(*tlelines))

        flout.write('# flags\n')
        flout.write('# mjdstart, mjdend, freq_mhz_start, freq_mhz_end, beam_start, beam_end\n')
        
        for ibeam in range(nbeams):
            for fcent, fwidth, bandname in const.freqs:
                beam_seps = traj.bsep[:, ibeam]
                # find start and end time when it's too close
                sep_times = times[beam_seps < np.radians(const.septhresh)]
                if len(sep_times) == 0:
                    continue
                    
                start_time = sep_times[0]
                end_time = sep_times[-1]
                start_mjd = start_time + mjd0
                end_mjd = end_time + mjd0
                fstart = fcent - fwidth/2.0
                fend = fcent + fwidth/2.0
                flout.write('{} {} {} {} {} {}\n'.format(start_mjd, end_mjd, fstart, fend, ibeam, ibeam))

    
    flout.flush()
    flout.close()
    if values.show:
        pylab.show()

if __name__ == '__main__':
    _main()




