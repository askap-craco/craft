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
from astropy.coordinates import SkyCoord
from astropy import units as u
from .crafthdr import DadaHeader
from aces.footprint_class import Footprint

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show figure')
    parser.add_argument('--psr', help='plot a pulsar. optinos: vela, 1644')
    parser.add_argument(dest='files', nargs='+')
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    psrs = {'1644':SkyCoord(ra=251.205355749999,dec=-45.98597222222222, unit='deg', frame='icrs'),
            'vela': SkyCoord(' 08:35:20.61149 -45:10:34.8751'.replace(':',' '), unit=('hour','deg'), frame='icrs')
            }
    if values.psr is not None:
        psr = psrs[values.psr]
        psrname = values.psr
    else:
        psr = None

    nearest_beams = []
    fp_shape = 'closepack36'
    fp_pitch = 0.9 # degrees
    fp_pa = 60 # degrees
    psr_beam_numbers = []

    for f in values.files:
        pylab.clf()
        hdr = DadaHeader.fromfile(f)
        beamra = np.array(list(map(float, hdr['BEAM_RA'][0].split(','))))
        nbeams = len(beamra)/2
        beamra = beamra[0:nbeams]
        beamdec = np.array(list(map(float, hdr['BEAM_DEC'][0].split(',')))[0:nbeams])
        beampos = [SkyCoord(ra=bra, dec=bdec, frame='icrs', unit='deg') for (bra, bdec) in zip(beamra, beamdec)]

        ra = float(hdr['RA'][0])
        dec = float(hdr['DEC'][0])
        ant_parangle = float(hdr['PAR_ANGLE'][0])
        pos = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')

        antra = float(hdr['ANT_RA'][0])
        antdec = float(hdr['ANT_DEC'][0])
        antpos = SkyCoord(ra=antra, dec=antdec, unit='deg', frame='icrs')


        fp = Footprint.named(fp_shape, np.radians(fp_pitch), np.radians(fp_pa +ant_parangle))
        fp.setRefpos(np.radians([ra, dec]))
        fp_ras = np.degrees([p.ra for p in fp.positions])
        fp_decs = np.degrees([p.dec for p in fp.positions])
        pylab.plot(beamra, beamdec, 'o', alpha=0.5)
        #pylab.plot(fp_ras, fp_decs, 'v', alpha=0.5)
        
        for ibeam, (bra, bdec) in enumerate(zip(beamra, beamdec)):
            pylab.text(bra, bdec, str(ibeam))

        if psr is not None:
            # todo: spherical angle
            distances = [psr.separation(bm).arcsec for bm in beampos]
            nearest_beam = np.argmin(distances)
            nearest_beams.append((beamra[nearest_beam], beamdec[nearest_beam]))

            for ib, bm in enumerate(beampos):
                print('Beam {} {} separation {} nearest? {}'.format(ib, bm.to_string('hmsdms'), distances[ib], ib ==nearest_beam))
        
            print(f, 'antcoord', antpos.to_string('hmsdms'), antra, antdec, 'parang %0.1f'%ant_parangle,'Nearest beam', nearest_beam, 'Distance', distances[nearest_beam], 'or', psr.separation(beampos[nearest_beam]).arcsec, 'arcsec', 'beam RA/dec', beampos[nearest_beam].to_string('hmsdms'), psr.to_string('hmsdms'))
            try:
                psrbeamno = int(hdr['FIELD_NAME'][0].split('_')[2].replace('beam',''))
                psr_beam_numbers.append(psrbeamno)
            except:
                psr_beam_numbers.append(None)
            
            pylab.plot(psr.ra.deg, psr.dec.deg, '^')
            pylab.text(psr.ra.deg, psr.dec.deg, psrname)


        pylab.plot([antra], [antdec], 'x', label='Ant position')
        pylab.plot([ra],[dec], '+', label='Field position')
        pylab.title(f)
        #pylab.savefig('%s.png' % f)

        if values.show:
            pylab.show()

    pylab.figure()
    nearest_beams = np.array(nearest_beams)

    fig, axes = pylab.subplots(2,1)
    axes = axes.flatten()
    ra_err = (nearest_beams[:, 0] - psr.ra.deg)*3600/np.cos(psr.dec.rad)
    dec_err = (nearest_beams[:, 1] - psr.dec.deg)*3600
    axes[0].plot(ra_err, marker='x')
    axes[0].set_ylabel('RA Offset (arcsec)')
    axes[1].plot(dec_err, marker='x')
    axes[1].set_ylabel('Dec offset (arcsec)')
    axes[1].set_xlabel('Beam number')

    pylab.figure()
    pylab.plot(ra_err, dec_err, 'x')
    for ibeam, (rae, dee, psr_beam_no) in enumerate(zip(ra_err, dec_err, psr_beam_numbers)):
        pylab.text(rae, dee, str(psr_beam_no))

    pylab.xlabel('RA Offset (arcsec)')
    pylab.ylabel('Dec Offset (arcsec)')
    pylab.savefig('offsetsxy.png')
                 
    pylab.show()

    

if __name__ == '__main__':
    _main()
