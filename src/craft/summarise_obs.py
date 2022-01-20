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
from .plot_fredda_cand import find_files
from .crafthdr import DadaHeader
import glob
import re
from astropy.coordinates import SkyCoord
from aces.footprint_class import Footprint
from aces.skypos import skypos

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def find_dada_header(f):
    ''' FInds the dada header  sitting next to a logfile'''
    d = os.path.dirname(f)
    hdrs = glob.glob(os.path.join(d, 'ak*.hdr.fixed'))
    assert len(hdrs) == 1, 'Couldnt find uniqe hdr in %s. Got %s' % (d, hdrs)
    
    return DadaHeader.fromfile(hdrs[0])

def mysep(pos1, pos2):
    sep =np.arccos(np.sin(pos1.dec.rad)*np.sin(pos2.dec.rad) + np.cos(pos1.dec.rad)*np.cos(pos2.dec.rad)*np.cos(pos1.ra.rad - pos2.ra.rad))
    return sep

class FreddaStats(object):
    fields = ('ncand','nblocks','nsamples','nseconds','nreal_time','freq_flag','freq_total','dm0_flag','dm0_total', 'nhits')
    def __init__(self):
        for f in FreddaStats.fields:
            setattr(self, f, 0)

    @staticmethod
    def fromfile(fname):
        stats = FreddaStats()
        
        with open(fname, 'rU') as f:
            for line in f:
                bits = line.split()
                if line.startswith('Found'):
                    stats.ncand = int(bits[1])
                elif line.startswith('Processed'):
                    stats.nblocks, stats.nsamples = list(map(int, (bits[1], bits[4])))
                    stats.nseconds = float(bits[7])
                    stats.nreal_time = float(bits[-3].replace('x',''))
                elif line.startswith('Freq auto-flagged'):
                    stats.freq_flag, stats.freq_total = list(map(int, bits[2].split('/')))
                elif line.startswith('DM0 auto-flagged'):
                    stats.dm0_flag, stats.dm0_total = list(map(int, bits[2].split('/')))
                else:
                    pass

        stats.nhits = 1

        return stats

    def __iadd__(self, other):
        for f in FreddaStats.fields:
            me = getattr(self, f)
            them = getattr(other, f)
            total = me + them
            setattr(self, f, total)

        return self
        
    def __str__(self):
        s = ' '.join('%s=%s' % (field, getattr(self, field)) for field in FreddaStats.fields)
        return s

    def values_aslist(self):
        l = [getattr(self, field) for field in FreddaStats.fields]
        return l

    __repr__ = __str__
        
class TargetParset(object):
    def __init__(self, f):

        self.sources = {}
        self._parse_sources(f)
        self._add_field_directions()

    def _parse_sources(self, pset):
        sources = self.sources
        with open(pset, 'rU') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if not '=' in line:
                    continue

                key, value = line.split('=')
                key = key.strip()
                keyparts = key.strip().split('.')
                if key.startswith('common.target.src'):
                    src = keyparts[2]

                    if src not in sources:
                        sources[src] = {'src':src}

                    srckey = keyparts[3]
                    sources[src][srckey] = value.strip()

        return sources

    def _add_field_directions(self):
        for src, data in self.sources.items():
            if 'field_direction' not in data:
                continue

            fdstr = data['field_direction']
            fdbits = fdstr.replace('[','').replace(']','').split(',')
            rastr, decstr, epochstr = fdbits
            skycoord = SkyCoord(rastr, decstr, unit='hour, deg')
            sp = skypos(rastr, decstr)
            data['skycoord'] = skycoord
            data['skypos'] = sp

            print(fdstr, rastr, decstr, skycoord, sp)

            pastr = data['pol_axis']
            pol_type, pol_angle = pastr.replace('[','').replace(']','').split(',')
            data['pol_type'] = pol_type
            data['pol_angle'] = float(pol_angle)

            fp = Footprint.named('square_6x6', np.radians(0.9), np.radians(float(pol_angle)+45))
            fp.setRefpos(sp)
            data['footprint'] = fp

    def get_nearest_source(self, pos):
        sources = [s for s in list(self.sources.values()) if 'skycoord' in s]
        nearest = min(sources, key=lambda p: mysep(p['skycoord'], pos))
        #for p in sources:
        #print p['field_name'], p['skycoord'], pos, pos.separation(p['skycoord']), np.degrees(mysep(pos, p['skycoord']))
        #sources = [s for s in self.sources.values() if 'skypos' in s]

        #str(pos) # need this to populate pos.ras and pos.decs

        #for p in sources:
        #    print pos, p['field_name'], p['skypos'], p['skypos'].distDeg(pos.ra, pos.dec)

        #nearest = min(sources, key=lambda p: p['skypos'].dist(pos.ra, pos.dec))

        return nearest

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

    all_stats = {}
    for log_file in values.files:
        try :
            hdr = find_dada_header(log_file)
            stats =  FreddaStats.fromfile(log_file)
            field = hdr['SOURCE'][0]
            
            if field not in all_stats:
                all_stats[field] = (stats, hdr)
            else:
                s = all_stats[field][0]
                s += stats
                print(log_file, field, s)
        except:
            logging.exception('COuld not load log file %s', log_file)
            raise


    fbybeam = open('stats_bybeam.txt','w')

    total_stats = FreddaStats()

    with open('stats_byfield.txt', 'w') as fout:
        fout.write('# name, ra, dec, %s\n' % (','.join(FreddaStats.fields)))
        for field in sorted(all_stats.keys()):
            (stats, hdr) = all_stats[field]
            total_stats += stats
            statlist = list(map(str, stats.values_aslist()))
            ra = float(hdr['RA'][0])
            dec = float(hdr['DEC'][0])
            fname = hdr['SOURCE'][0]
            beam_ra = list(map(float, hdr['BEAM_RA'][0].split(',')))
            beam_dec = list(map(float, hdr['BEAM_DEC'][0].split(',')))
            dout = list(map(str, [fname, ra, dec]))
            dout.extend(statlist)
            fout.write(' '.join(dout))
            fout.write('\n')

            for bra, bdec in zip(beam_ra, beam_dec):
                dout = list(map(str, [fname, ra, dec, bra, bdec]))
                dout.extend(statlist)
                fbybeam.write(' '.join(dout))
                fbybeam.write('\n')

        fout.write('# TOTALS %s\n' % total_stats)
        fbybeam.write('# TOTALS %s\n' % total_stats)
        
        print('TOTALS', total_stats)

if __name__ == '__main__':
    _main()
