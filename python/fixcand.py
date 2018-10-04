#!/usr/bin/env python
"""
Fix candidate files - make MJD referenced to infinity

Copyright (C) CSIRO 2018
"""
import matplotlib as mpl
import numpy as np
import logging
import sys

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-b','--ref-bandwidth', type=float, help='Reference bandwidth', default=336.0)
    parser.add_argument('-f', '--ref-freq', type=float, help='Low end of the band that s the reference frequency of the MJD (MHz)', default=(1488.-336.))
    parser.add_argument('-e','--freq-offset', type=float, help='Frequency offset to be applied to DM (MHz)', default=0.0)
    parser.add_argument('-s', '--suffix', default='.finf', help='Suffix for the output files')
    parser.add_argument('-t','--toff', type=float, help='Milliseconds to adjust MJD by', default=0.0)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        fix(f, values)

def fix(f, values):
    ref_freq = values.ref_freq/1e3 # Ghz
    ref_bw = values.ref_bandwidth/1e3 # GHz
    ftop = ref_freq + ref_bw
    fbot = ref_freq
    foff = values.freq_offset
    new_file = f + values.suffix
    fout = open(new_file, 'w')
    for line in open(f, 'rU'):
        if line.startswith('#'):
            fout.write(line)
            fout.write('# processed with ' + ' '.join(sys.argv) + '\n')
            continue
            
        line = line.strip()
        bits = line.split()
        dm = float(bits[5])
        mjd = float(bits[7])
        offset_ms = 4.15*dm*(0 - ref_freq**-2) + values.toff
        offset_days = offset_ms/ 86400. / 1e3
        assert offset_days < 0
        bits[7] = '%0.15f'%(mjd + offset_days)

        # time delay in milliseconds at the reference frequency in ms
        dt = 4.15*dm*(ftop**-2 - fbot**-2)
        dmnew = dt/4.15/((ftop + foff)**-2 - (fbot + foff)**-2)
        bits[5] = '%0.2f' % dmnew
        
        fout.write(' '.join(bits) + '\n')

    fout.close()

if __name__ == '__main__':
    _main()
