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
import aplpy
from astropy.coordinates import SkyCoord
from astropy import units as u
from crafthdr import DadaHeader
from aces.footprint_class import Footprint

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show figure')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    psr = SkyCoord(ra=251.20534,dec=-45.98597, unit='deg', frame='icrs')
    psrname ='1644'
    nearest_beams = []
    fp_shape = 'square_6x6'
    fp_pitch = 0.9 # degrees
    fp_pa = 45 # degrees

    for f in values.files:
        pylab.clf()
        hdr = DadaHeader.fromfile(f)
        beamra = np.array(map(float, hdr['BEAM_RA'][0].split(',')))
        nbeams = len(beamra)/2
        beamra = beamra[0:nbeams]
        beamdec = np.array(map(float, hdr['BEAM_DEC'][0].split(','))[0:nbeams])
        antra = float(hdr['RA'][0])
        antdec = float(hdr['DEC'][0])
        ant_parangle = float(hdr['PAR_ANGLE'][0])
        antpos = SkyCoord(ra=antra, dec=antdec, frame='icrs', unit='deg')
        beampos = [SkyCoord(ra=bra, dec=bdec, frame='icrs', unit='deg') for (bra, bdec) in zip(beamra, beamdec)]

        fp = Footprint.named(fp_shape, np.radians(fp_pitch), np.radians(fp_pa + ant_parangle))
        fp.setRefpos(np.radians([antra, antdec]))
        fp_ras = np.degrees([p.ra for p in fp.positions])
        fp_decs = np.degrees([p.dec for p in fp.positions])
        pylab.plot(beamra, beamdec, 'o', alpha=0.5)
        pylab.plot(fp_ras, fp_decs, 'v', alpha=0.5)
        
        for ibeam, (bra, bdec) in enumerate(zip(beamra, beamdec)):
            pylab.text(bra, bdec, str(ibeam))

        # todo: spherical angle
        distances = [psr.separation(bm).deg for bm in beampos]
        nearest_beam = np.argmin(distances)
        nearest_beams.append((beamra[nearest_beam], beamdec[nearest_beam]))
        
        print 'antcoord', antpos.to_string('hmsdms'), antra, antdec, 'parang %0.1f'%ant_parangle,'Nearest beam', nearest_beam, 'Distance', distances[nearest_beam]*3600, 'or', psr.separation(beampos[nearest_beam]).arcsec, 'arcsec', 'beam RA/dec', beampos[nearest_beam].to_string('hmsdms'), psr.to_string('hmsdms')

        pylab.plot([antra], [antdec], 'x')
        pylab.plot(psr.ra.deg, psr.ra.deg, '^')
        pylab.text(psr.ra.deg, psr.dec.deg, psrname)
        pylab.title(f)
        #pylab.savefig('%s.png' % f)

        if values.show:
            pylab.show()

    pylab.figure()
    nearest_beams = np.array(nearest_beams)

    fig, axes = pylab.subplots(2,1)
    axes = axes.flatten()
    ra_err = (nearest_beams[:, 0] - psr.ra.deg)*60
    dec_err = (nearest_beams[:, 1] - psr.dec.deg)*60
    axes[0].plot(ra_err, marker='x')
    axes[0].set_ylabel('RA Offset (arcmin)')
    axes[1].plot(dec_err, marker='x')
    axes[1].set_ylabel('Dec offset (arcmin)')
    axes[1].set_xlabel('Beam number')

    pylab.figure()
    pylab.plot(ra_err, dec_err, 'x')
    for ibeam, (rae, dee) in enumerate(zip(ra_err, dec_err)):
        pylab.text(rae, dee, str(ibeam))

    pylab.xlabel('RA Offset (arcmin)')
    pylab.ylabel('Dec Offset (arcmin)')
                 
    pylab.show()

    

if __name__ == '__main__':
    _main()
