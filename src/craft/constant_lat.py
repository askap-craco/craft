#!/usr/bin/env python
"""
Makes a parset wth constant latitude

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from aces import mro

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    lats = [0,-5,5] # degrees
    step = 5. # degrees
    fields = []
    mro = EarthLocation(lat=-26.98*u.deg, lon=116.53*u.deg, height=360*u.m)
    utcoffset = 8*u.hour  # Eastern Daylight Time
    time = Time('2012-7-12 23:00:00') - utcoffset
    galann = open('gal.ann', 'w')
    eqann = open('eq.ann', 'w')

    for lat in lats:
        longstep = step/np.cos(np.radians(lat)) # degrees
        longs = np.arange(-180, 180, longstep)
        radius = step/2.

        for lon in longs[0:-1]:
            cgal = SkyCoord(lon, lat, frame='galactic', unit='deg')
            ceq = cgal.fk5
            if ceq.dec.degree > 40.:
                continue

            galann.write('CIRCLE {} {} {}\n'.format(lon, lat, radius)) # Galactic
            eqann.write('CIRCLE {} {} {}\n'.format(ceq.ra.degree, ceq.dec.degree, radius)) # equatorial
            fields.append((ceq, cgal))


    pset = open('constant_lat.parset','w')
    pset.write('# Made by constant_lat.py- KB 2016\n')
    pset.write('common.enable_cp = false\n')
    pset.write('common.target.src%d.autoant = true\n')
    pset.write('common.target.src%d.duration = 600\n')
    pset.write('common.target.src%d.sky_frequency = 1278.5\n')

    #pset.write('common.target.src%d.footprint.name=square_6x6\n')
    #pset.write('common.target.src%d.footprint.pitch=0.9\n')
    #pset.write('common.target.src%d.footprint.pa=45\n')
    targnames = ['src{}'.format(i+1) for i in range(len(fields))]
    pset.write('common.targets = [[{}]]\n'.format(','.join(targnames)))
    pset.write('standard.repeat = 999\n')
    
    pset.write('common.antennas = [ant2,ant4,ant5,ant10,ant12,ant13,ant14,ant16,ant24,ant27,ant28,ant30]\n')

    for ifield, field in enumerate(fields):
        target = ifield+1
        ceq, cgal = field
        if ceq.dec.degree > -26.98:
            pa = 170 # limits are +- in fcm tos.opl.limits
        else:
            pa = 0
        pset.write('common.target.src{}.field_name = G{:3.0f}{:+2.0f}\n'.format(target, cgal.l.degree, cgal.b.degree))
        pset.write('common.target.src{}.field_direction = [{}, {}, J2000]\n'.format(target,
                                                                                    ceq.ra.to_string(u.hour, sep=':'),
                                                                                    ceq.dec.to_string(u.deg, sep=':')))
        pset.write('common.target.src{}.pol_axis = [pa_fixed, {}]\n'.format(target, pa))
        
        

    midnight = Time('2016-07-08 00:00:00')
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    nup = np.zeros(len(delta_midnight))
    times = midnight + delta_midnight
    frame = AltAz(obstime=times, location=mro)
    for field, cgal in fields:
        altaz = field.transform_to(frame)
        print(altaz.alt.degree)
        #pylab.plot(times.plot_date, altaz.alt.degree)
        print(field, altaz)
        nup[altaz.alt.degree > 15.0] += 1


    pylab.ylabel('Nummber of fields above horizon')
    pylab.xlabel('Date (UT)')
    pylab.plot(times.plot_date, nup)
    pylab.savefig('nup.png')


if __name__ == '__main__':
    _main()
