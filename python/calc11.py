#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from collections import OrderedDict

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class CalcFile(object):
    def __init__(self):
        self.cards = OrderedDict()

    def __iadd__(self, v):
        self.cards[v[0]] = v[1]
        return self

    def addtel(self, cname, vkey=None, index=None, default=None):
        antdata = self.ant_data[self.telno]
        if default:
            value = default
        else:
            value = antdata.get(vkey, default)

        print antdata
        assert value is not None, 'Unknown value for ant {} key {} {}'.format(self.telno, cname, vkey)
        if index is not None:
            bits = value.replace('[','').replace(']','').split(',')
            value = bits[index].strip()

        cardkey = 'TELESCOPE {} {}'.format(self.itel, cname)
        self += (cardkey, value)

    def add_antdata(self, ant_data, antnos):
        self += ('NUM TELESCOPES', len(antnos))
        self.ant_data = ant_data
        for itel, telno in enumerate(antnos):
            teld = self.ant_data[telno]
            self.itel = itel
            self.telno = telno
            self.teld = teld
            self.addtel('NAME', 'name')
            self.addtel('MOUNT', default='AZEL')
            self.addtel('OFFSET (m)', default='0.00000')
            self.addtel('X (m)', 'location.itrf', 0)
            self.addtel('Y (m)', 'location.itrf', 1)
            self.addtel('Z (m)', 'location.itrf', 2)
            self.addtel('SHELF', default='NONE')

    def writeto(self, foutname):
        with open(foutname, 'w') as fout:
            for k,v in self.cards.iteritems():
                kcol = k +':'
                fout.write('{:<20}{}\n'.format(kcol,v))

class Poly(object):
    def __init__(self, polyid):
        self.polyid = polyid
        self.source0antpolys = {}
        self.mjd = None
        self.sec = None

    def update(self, name, value):
        if 'MJD' in name:
            self.mjd = int(value)
        if 'SEC' in name:
            self.sec = int(value)
        if name.startswith('SRC 0'):
            self.mjdfull = float(self.mjd) + float(self.sec)/86400.
            namebits = name.split()
            antid = int(namebits[3])
            polyname = ' '.join(namebits[4:])
            antpolys = self.source0antpolys.get(antid, {})
            self.source0antpolys[antid] = antpolys
            antpolys[polyname] = map(float, value.split())

    def __str__(self):
        return 'POLY ID={} mjd={} sec={} ants={}'.format(self.polyid, self.mjd, self.sec, str(self.source0antpolys.keys()))

    __repr__ = __str__

class Source(object):
    def __init__(self, srcid):
        self.src = srcid
        self.polys = []

class Scan(object):
    def __init__(self, scanid, resfile):
        self.resfile = resfile
        self.scanid = scanid
        self.sources = []
        self.curr_poly = None
        self.polys = []

    def update(self, name, value):
        namebits = name.split()
        if 'POINTING SRC' in name:
            self.pointing_source = value
        if 'NUM PHS CTRS' in name:
            self.num_phase_centers = int(value)
            assert self.num_phase_centers == 1 # limitation for now
        if 'PHS CTR 0 SRC' in name:
            self.phase_center_source = value
        if namebits[2] == 'POLY':
            polyid = int(namebits[3])
            if self.curr_poly is None or self.curr_poly.polyid != polyid:
                self.curr_poly = Poly(polyid)
                self.polys.append(self.curr_poly)

            self.curr_poly.update(name, value)
        if namebits[0] == 'SRC':
            self.curr_poly.update(name, value)

    def get_poly(self, mjd):
        the_poly = min(self.polys, key=lambda p: abs(mjd - p.mjdfull))
        return the_poly

    def eval_src0_poly(self, mjd):
        poly = self.get_poly(mjd)
        secs = mjd - poly.mjdfull
        ant_results = {}
        for ant, polys in poly.source0antpolys.iteritems():
            antname = self.resfile.telnames[ant]
            ant_results[antname] = {}
            for polyname, pcoeff in polys.iteritems():
                ant_results[antname][polyname] = np.polyval(pcoeff[::-1], secs)

        return ant_results

class ResultsFile(OrderedDict):
    def __init__(self, fname):
        OrderedDict.__init__(self)
        self.fname = fname
        self.scans = []
        curr_scan = None

        with open(self.fname, 'rU') as f:
            for line in f:
                name, value = [s.strip() for s in line.split(':')]
                if name.startswith('SCAN'):
                    namebits = name.split()
                    scanid = int(namebits[1])
                    if curr_scan is None or scanid != curr_scan.scanid:
                        curr_scan = Scan(scanid, self)
                        self.scans.append(curr_scan)
                    curr_scan.update(name, value)
                elif name.startswith('SRC'):
                    curr_scan.update(name, value)
                else:
                    self[name] = value

        self.num_telsecopes = int(self['NUM TELESCOPES'])
        self.telnames = []
        for itel in xrange(self.num_telsecopes):
            self.telnames.append(self['TELESCOPE {} NAME'.format(itel)])

        self.num_scans = int(self['NUM SCANS'])

    def get_fringe_params(self, scanid, srcname, antname, mjd):
        polynames = ('DELAY (us)','U (m)', 'V (m)', 'W (m)')
        polyid = 0
        for polyname in polynames:
            rfile.scans[0].polys[0].source0antpolys.iteritems()


def plot_polys(rfile):
    polynames = ('DELAY (us)','U (m)', 'V (m)', 'W (m)')
    polyid = 0
    for polyname in polynames:
        pylab.figure()
        p0 = None
        print  rfile.scans[0].polys[0]
        for ant, polys in rfile.scans[0].polys[0].source0antpolys.iteritems():
            pcoeff = polys[polyname]
            t = np.arange(120.)
            pvalue = np.polyval(pcoeff[::-1], t)
            if p0 is None:
                p0 = pvalue
            pylab.plot(t, pvalue - p0, label=str(ant))
        pylab.title(polyname)
    pylab.show()



def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = ResultsFile(values.files[0])
    mjd = 58154 + 39360./86400.
    results = f.scans[0].eval_src0_poly(mjd)
    print 'TELNAMES', f.telnames
    print results
    plot_polys(f)



if __name__ == '__main__':
    _main()
