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
import warnings
from astropy.utils import iers
from astropy import units as u
from astropy.time import Time
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import subprocess

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def leapSecondsAt(t: Time):
    '''
    Returns the leap second value at a particular Astropy time
    '''
    dtaiutc = (t.tai.datetime - t.utc.datetime).seconds
    return dtaiutc

def getpoly(polyvals, ant,key):
    return np.array([p[ant][key] for p in polyvals])

def get_uvw(polyvals, ant):
    u = getpoly(polyvals,ant,'U (m)')
    v = getpoly(polyvals,ant,'V (m)')
    w = getpoly(polyvals,ant,'W (m)')
    return np.array([u,v,w])

def get_ant_uvw(telnames, polyvals):
    return np.array([get_uvw(polyvals, ant) for ant in telnames]).transpose(0,2,1)

def uvw_ant_to_baselines(uvw):
    nant = uvw.shape[0]
    bl = []
    for iant1 in range(nant):
        for iant2 in range(iant1+1,nant):
            bl.append(uvw[iant1,...] - uvw[iant2,...])
    return np.array(bl)


class CalcFile(object):
    calc_template =  '''
JOB ID:             374055
JOB START TIME:     {start.utc.mjd}
JOB STOP TIME:      {stop.utc.mjd}
DUTY CYCLE:         1.000
OBSCODE:            CRAFTFR
DIFX VERSION:       trunk
SUBJOB ID:          0
SUBARRAY ID:        0
VEX FILE:           /fred/oz002/adeller/askap/frb190608/0407/data/c1_f0/craftfrb.vex
START MJD:          {start.utc.mjd}
START YEAR:         {start.utc.datetime.year}
START MONTH:        {start.utc.datetime.month}
START DAY:          {start.utc.datetime.day}
START HOUR:         {start.utc.datetime.hour}
START MINUTE:       {start.utc.datetime.minute}
START SECOND:       {start.utc.datetime.second}
SPECTRAL AVG:       1
TAPER FUNCTION:     UNIFORM
NUM SOURCES:        1
SOURCE 0 NAME:      {name}
SOURCE 0 RA:        {skycoord.ra.rad}
SOURCE 0 DEC:       {skycoord.dec.rad}
SOURCE 0 CALCODE:    
SOURCE 0 QUAL:      0
NUM SCANS:          1
SCAN 0 IDENTIFIER:  No0001
SCAN 0 START (S):   0
SCAN 0 DUR (S):     {duration}
SCAN 0 OBS MODE NAME:askap
SCAN 0 UVSHIFT INTERVAL (NS):2000000000
SCAN 0 AC AVG INTERVAL (NS):2000000
SCAN 0 POINTING SRC:0
SCAN 0 NUM PHS CTRS:1
SCAN 0 PHS CTR 0:   0
NUM SPACECRAFT:     0
IM FILENAME:        /fred/oz002/adeller/askap/frb190608/0407/data/c1_f0/craftfrb_374055.im
FLAG FILENAME:      /fred/oz002/adeller/askap/frb190608/0407/data/c1_f0/craftfrb_374055.flag
'''
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

        logging.debug('%s %s %s', vkey, self.telno, value)
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

    def add_eops(self, eoptab=None, toff=3):
        '''
        Add entries from the given EOP table to this CalcFile
        if eoptab is None, it gets the astropy.ieors.earth_origientation_table()
        and uses the START MJD from the table to populate the CalcFile
        :eoptab: Table of EOP values -= we use 'MJX', 'PM_x', PM_y'
        :toff: number of days either side of START_MJD to find entries for if eoptab is None
        '''


        neops = self.cards.get('NUM EOPS', 0)
        assert neops == 0, f'CalcFile already has {neops} eops'
        
        if eoptab is None:
            if 'START MJD' not in self.cards.keys():
                raise ValueError('No START MJD in file. Cant add eops')
                
            start = Time(float(self.cards['START MJD']), scale='utc', format='mjd')
            eoptab = iers.earth_orientation_table.get()
            eoptab = eoptab[(start.mjd -toff < eoptab['MJD'].value) & (eoptab['MJD'].value < start.mjd+toff)]
        
        assert len(eoptab) > 0, 'Not enough EOP table entries'

        self += ('NUM EOPS', int(len(eoptab)))
        for irow, row in enumerate(eoptab):
            mjd = Time(row['MJD'], scale='utc', format='mjd')
            self += (f'EOP {irow} TIME (mjd)', str(int(mjd.value)))
            self += (f'EOP {irow} TAI_UTC (sec)',str(leapSecondsAt(mjd)))
            self += (f'EOP {irow} UT1_UTC (sec)', str(row['UT1_UTC'].to(u.second).value))
            self += (f'EOP {irow} XPOLE (arcsec)', str(row['PM_x'].to(u.arcsecond).value))
            self += (f'EOP {irow} YPOLE (arcsec)', str(row['PM_y'].to(u.arcsecond).value))

        return self

    def __str__(self):
        s = ''
        for k,v in self.cards.items():
            kcol = k.strip() +':'
            s += '{:<20}{}\n'.format(kcol,v) # DIFX looks at first 20 rows for keyword

        return s
        

    def writeto(self, foutname):
        with open(foutname, 'w') as fout:
            fout.write(str(self))

    @staticmethod
    def fromstring(s):
        f = CalcFile()
        for line in s.split('\n'):
            bits = line.split(':')
            if line.strip().startswith('#'):
                continue
            
            if len(bits) == 2:
                f += (bits[0].strip(), bits[1].strip())

        return f

    @staticmethod
    def from_time_src(start:Time, stop:Time, skycoord:SkyCoord, name:str):
        '''
        Creates basic calcfile from the given start, stop and source coordinates
        uses the built-in calc_template string, so it'll have some weird old stuff
        but I'm too scard to work out what to remove.
        
        :start: Astropy start time
        :stop: Astropy stop time
        :name: source name
        :skycoord: Astropy SkyCoord
        '''
        duration = (stop.datetime - start.datetime).seconds
        s = CalcFile.calc_template.format(start=start, stop=stop, name=name, skycoord=skycoord, duration=duration)
        return CalcFile.fromstring(s)


def run_calc(calc_file:CalcFile, filename:str):
    '''
    Runs difxcalc on the given calcfile with the name supplied.
    Returns a calc ResultsFile from the result
    :calc_file: CalcFile with all the bits in it.
    :filename: Name of the file to write the file to. must end in '.calc'
    e.g. myscan.calc. difxcalc will write the results to a file with the same root and ending in .im
    '''
    assert filename.endswith('.calc'), f'Invalid calc file name {filename}'
    outfile = filename.replace('.calc','.im')
    try:
        os.remove(outfile)
    except FileNotFoundError:
        pass

    calc_file.writeto(filename)
    subprocess.call(['difxcalc',filename])

    return ResultsFile(outfile)


class Poly(object):
    def __init__(self, polyid):
        self.polyid = polyid
        self.source0antpolys = {}
        self.mjd = None
        self.sec = None

    @property
    def mjdfull(self):
        return self.mjd + self.sec/86400.

    def update(self, name, value):
        if 'MJD' in name:
            self.mjd = int(value)
        if 'SEC' in name:
            self.sec = int(value)
        if name.startswith('SRC 0'):
            namebits = name.split()
            antid = int(namebits[3])
            polyname = ' '.join(namebits[4:])
            antpolys = self.source0antpolys.get(antid, {})
            self.source0antpolys[antid] = antpolys
            antpolys[polyname] = list(map(float, value.split()))

    def __str__(self):
        return 'POLY ID={} mjd={} sec={} ants={}'.format(self.polyid, self.mjd, self.sec, str(list(self.source0antpolys.keys())))

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

    @property
    def first_mjd(self):
        return min(self.polys, key=lambda p:p.mjdfull).mjdfull

    @property
    def last_mjd(self):
        return max(self.polys, key=lambda p:p.mjdfull).mjdfull # COuld add interval but don't care right at the moment

    def update(self, name, value):
        namebits = name.split()
        if 'POINTING SRC' in name:
            self.pointing_source = value # source name = e.g. M87
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
        #the_poly = min(self.polys, key=lambda p: abs(mjd - p.mjdfull))
        after_polys = [p for p in self.polys if p.mjdfull <= mjd]
        if len(after_polys) == 0:
            the_poly = max(self.polys, key=lambda p:p.mjdfull)
            warnings.warn('Past last polynomial. mjd={} last mjd={}'.format(mjd, the_poly.mjdfull))
        else:
            the_poly = max(after_polys, key=lambda p:p.mjdfull)

    
        return the_poly

    def eval_src0_poly(self, mjd):
        ''' evaluates polynomials for given mjd for source 0
        :returns: Dictionry of dictionaries. First key: antenna name. Second key: polynomialname
        '''
        poly = self.get_poly(mjd)
        secs = (mjd - poly.mjdfull)*86400.
        
        if secs < 0.:
            warnings.warn('ERR Dodgey offset. mjd={} polymjd ={} secoffset={}'.format(mjd, poly.mjdfull, secs))
        assert secs >= 0
            
        ant_results = {}
        for ant, polys in poly.source0antpolys.items():
            antname = self.resfile.telnames[ant]
            ant_results[antname] = {}
            for polyname, pcoeff in polys.items():
                ant_results[antname][polyname] = np.polyval(pcoeff[::-1], secs)

            elevation = ant_results[antname]['EL GEOM']
            if abs(elevation) < 2:
                warnings.warn('ERR: Evaluating geoemtry at very low elevation. Be careful! Elevation={} mjd={} poly={}'.format(elevation, mjd, polys))

        ant_results['secoff'] = secs
        ant_results['poly'] = poly

        return ant_results

    def eval_src0_poly_delta(self, mjd, refant):
        res = self.eval_src0_poly(mjd)
        ref_results = res[refant]
        resdelta = {}
        for ant, polydata in res.items():
            resdelta[ant] = {}
            for polyname, value in polydata.items():
                resdelta[ant][polyname] = value - ref_results[polyname]

        return resdelta

    def eval_src0_poly_table(self, time:Time):
        '''
        Returns a numpy array of the values at the given astropy time
        Instead of the rediculous dictionary version it returns an astropy QTable
        Units should be correctly set
        column names should have units stripped and be lower case
        UTC of time is used as that's what used in calc11
        
        '''
        polyvals = [self.eval_src0_poly(m) for m in time.utc.value]
        telnames = self.resfile.telnames
        p0 = polyvals[0]
        v0 = list(p0.values())[0]
        colnames = list(v0.keys())
        d  = QTable()
        d['antennas'] = telnames
        for cname in colnames:
            if cname.endswith('(us)'):
                unit = u.microsecond
            elif cname.endswith('(m)'):
                unit = u.meter
            elif cname == 'AZ' or cname == 'EL GEOM':
                unit = u.degree
            else:
                raise ValueError(f'Unknown unit for field {cname}')
        
        
            colname = cname.split(' ')[0].lower()
            d[colname] = np.array([getpoly(polyvals, ant, cname) for ant in telnames])*unit

        return d
        

class ResultsFile(OrderedDict):
    def __init__(self, fname):
        OrderedDict.__init__(self)
        self.fname = fname
        self.scans = []
        curr_scan = None

        with open(self.fname, 'rU') as f:
            for line in f:
                if ':' not in line:
                    continue
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
        for itel in range(self.num_telsecopes):
            self.telnames.append(self['TELESCOPE {} NAME'.format(itel)])

        self.num_scans = int(self['NUM SCANS'])

    def get_fringe_params(self, scanid, srcname, antname, mjd):
        polynames = ('DELAY (us)','U (m)', 'V (m)', 'W (m)')
        polyid = 0
        for polyname in polynames:
            iter(r.scans[0].polys[0].source0antpolys.items())


def plot_polys(rfile, tmax):
    polynames = ('DELAY (us)','U (m)', 'V (m)', 'W (m)')
    polyid = 0
    for polyname in polynames:
        pylab.figure()
        p0 = None
        print(r.scans[0].polys[0])
        for ant, polys in rfile.scans[0].polys[0].source0antpolys.items():
            pcoeff = polys[polyname]
            t = np.arange(tmax)
            pvalue = np.polyval(pcoeff[::-1], t)
            if p0 is None:
                p0 = pvalue
            pylab.plot(t, pvalue - p0, label=str(ant))
        pylab.title(polyname)
    pylab.show()


def plot_polys_range(rfile, mjdstart, tmax):
    '''
    :rfile: Calcfile
    :mjdstart: MJD
    :tmax: seconds
    '''
    toff = np.arange(tmax, step=60.)
    mjds = toff/86400. + mjdstart
    values = []
    for imjd, mjd in enumerate(mjds):
        p = rfile.scans[0].eval_src0_poly(mjd) 
        values.append(p)
        if imjd == 563 or imjd == 560:
            print('MJD set', mjd)

        
    polynames = ('DELAY (us)','U (m)', 'V (m)', 'W (m)')
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = pylab.subplots(6,1, sharex=True)
    lines = []
    labels = []
    print('At MJD', mjdstart)
    
    for ia1, a1 in enumerate(rfile.telnames):
        for ia2, a2 in enumerate(rfile.telnames[ia1:]):
            lbl = '{}-{}'.format(a1, a2)
            labels.append(lbl)
            u = [val[a1]['U (m)'] - val[a2]['U (m)'] for val in values]
            v = [val[a1]['V (m)'] - val[a2]['V (m)'] for val in values]
            w = [val[a1]['W (m)'] - val[a2]['W (m)'] for val in values]
            el = [val[a1]['EL GEOM'] for val in values]
            secoffs = [val['secoff'] for val in values]
            delays = [val[a1]['DELAY (us)'] - val[a2]['DELAY (us)']for val in values]
            u,v,w,secoffs, delays = list(map(np.array, (u,v,w, secoffs, delays)))
            print('{} u={} v={} w={} el={} secoff={} delay={:0.9f}us'.format(lbl, u[0], v[0], w[0], el[0], secoffs[0], float(delays[0])))
            if np.any(abs(u) > 1e5):
                bad_times = np.where(abs(u) > 1e5)[0]
                print('BAD TIMES',bad_times, mjds[bad_times])
                print(values[bad_times[0]]['poly'])
                
            pylab.figure(2)
            line, = pylab.plot(u,v)
            pylab.plot(-u,-v, color=line.get_color())

            ax1.plot(toff, u)
            ax2.plot(toff, v)
            ax3.plot(toff, w)
            ax4.plot(toff, secoffs)
            ax5.plot(toff, el)
            l, = ax6.plot(toff, delays)
            lines.append(l)



    ax1.set_ylabel('U (m)')
    ax2.set_ylabel('V (m)')
    ax3.set_ylabel('W (m)')
    ax4.set_ylabel('Secoff')
    ax5.set_ylabel('Elevation')
    ax6.set_ylabel('Delay (us)')
    ax6.set_xlabel('Seconds from MJD {:.5f}'.format(mjdstart))
    fig.legend(lines, labels, 'upper right')

    fig.savefig(rfile.fname + '.uvt.png')
    pylab.xlabel('U (m)')
    pylab.ylabel('V (m)')
    pylab.grid(True)
    pylab.savefig(rfile.fname + '.uv.png')
    pylab.figlegend(lines, labels, 'upper right')
    pylab.show()

    

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--mjdstart', type=float, help='MJD start to plot')
    parser.add_argument('--nhrs', type=float, help='Number of hours to plot', default=1.0)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = ResultsFile(values.files[0])
    #mjd = 58154 + 39360./86400.
    #results = f.scans[0].eval_src0_poly(mjd)
    #plot_polys(f,10)
    mjdstart = f.scans[0].first_mjd
    mjdend = f.scans[0].last_mjd
    print('File {} starts at mjd {} and ends at {}'.format(values.files[0], mjdstart, mjdend))
    if values.mjdstart is None:
        mjd = mjdstart
    else:
        mjd = values.mjdstart

    if values.nhrs is None:
        nhrs = (mjdend - mjdstart)*24.
    else:
        nhrs = values.nhrs

    print('PLotting mjds starting from', mjd)
    plot_polys_range(f, mjd, 3600*values.nhrs)
    




if __name__ == '__main__':
    _main()
