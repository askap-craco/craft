#!/usr/bin/env python
"""
conver uvfits file to filterbank
Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
#from . import sigproc
#from astropy.io import fits
#from craco import uvfits_meta
from craft import uvfits
from craft.sigproc import SigprocFile
from craft.craco import bl2ant
import IPython

def get_bl_selector(values, baseline_order):
    bl_selector = []
    data_selector = []
    outnames = []
    if len(values.ants) > 0:
        myants = values.ants
    else:
        myants = np.arange(1, 37, 1)

    for ii, blid in enumerate(baseline_order):
        a1, a2 = bl2ant(blid)
        if a1 in myants or a2 in myants:
            if a1 != a2:
                if values.no_cross:
                    continue
                outnames.append(f"{values.uv.strip('.uvfits')}_cross_{a1:02g}_{a2:02g}.fil")

            else:
                if values.no_autos:
                    continue
                outnames.append(f"{values.uv.strip('.uvfits')}_auto_{a1:02g}_{a2:02g}.fil")

            bl_selector.append(blid)
            data_selector.append(ii)

    return np.array(bl_selector), np.array(data_selector), np.array(outnames)

def get_tsamp_uvfits(f):
    '''
    Opens the uvfits file and reads the INTTIM from the first visrow
    '''
    blocker = f.fast_raw_blocks(nt=1)
    return next(blocker)[0][0]['INTTIM']    #in seconds

def get_fil_header(f):
    '''
    Composes the information needed for a filterbank header from a uvfits file
    '''
    hdr = f.header
    start_date_jd = hdr['PZERO4']
    start_mjd = start_date_jd - 2400000.5
    tsamp = get_tsamp_uvfits(f)

    fcent = hdr['CRVAL4'] / 1e6 # MHz
    foff = hdr['CDELT4'] / 1e6 # MHz
    ref_channel = hdr['CRPIX4'] # pixels
    nchan = hdr['NAXIS4']

    fch1 = fcent - foff*(nchan - ref_channel)
    hdr = {'fch1':fch1, 'foff':foff,'tsamp':tsamp, 'tstart':start_mjd, 'nbits':32, 'nifs':1, 'nchans':nchan, 'src_raj':0.0, 'src_dej':0.0}
    return hdr

def open_outfiles(outnames, hdr):
    outfiles = []
    for o in outnames:
        outfiles.append(SigprocFile(o, 'w', hdr))

    return outfiles

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Convert uvfits file to filterbank', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='show plots')
    parser.add_argument('-uv', '--uv', type=str, help="Path to the uvfits file to read from", required=True)
    parser.add_argument('-mf', '--metadata-file', type=str, help="Path to the metadata file if needed (optional)", default=None)
    parser.add_argument('--nt', type=int, help="nt samples to read at a time", default=256)
    parser.add_argument('--mode', type=str, help="Mode for converting complex to real - ['real' or 'abs']", default='abs')
    parser.add_argument('--no-cross', action='store_false', default=False)
    parser.add_argument('--no-autos', action='store_false', default=False)
    parser.add_argument('--ants', nargs='+', type=int, help="List of antennas to write out (def:all)", default = [])
    parser.add_argument('--maxfiles', help='max files to open', default=np.inf, type=int)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    log = logging.getLogger(__name__)
    log.debug(f"Values - {values}")
    #f = uvfits_meta.open(values.uv, metadata_file = values.metadata_file)
    log.info(f"Openning the uvfits file - {values.uv}")
    f = uvfits.open(values.uv)
    log.info(f"Getting baselines to read")
    bl_mask, data_mask, fil_outnames = get_bl_selector(values, f.raw_baseline_order)

    log.debug(f"Selected {len(bl_mask)} baselines. Blids are - {bl_mask}")

    log.info(f"Restricting the outfiles to {values.maxfiles}")
    if values.maxfiles < len(fil_outnames):
        bl_mask = bl_mask[:values.maxfiles]
        data_mask = data_mask[:values.maxfiles]
        fil_outnames = fil_outnames[:values.maxfiles]
    
    log.info("Creating the filterbank header")
    fil_hdr = get_fil_header(f)
    log.debug(f"Filterbank header is - \n{fil_hdr}")

    log.info(f"Opening {len(fil_outnames)} filterbank files")
    fouts = open_outfiles(fil_outnames, fil_hdr)

    try:
        for iblk, blk in enumerate(f.fast_time_blocks(nt = values.nt)):
            log.debug(f"iblk -> {iblk}")        
            for ii in range(len(fouts)):
                #IPython.embed()
                block = np.squeeze(blk[0].filled(fill_value = 0))
                data = block[data_mask][ii]
                if values.mode == 'real':
                    spec = data.real
                elif values.mode == 'abs':
                    spec = np.abs(data)
                elif values.mode == 'imag':
                    spec = data.imag
                else:
                    raise ValueError("values.mode is invalid = {valus.mode}, expected [real, imag, abs]")
                
                spec.T.astype(np.float32).tofile(fouts[ii].fin)
    except KeyboardInterrupt:
        pass
    finally:
        for fout in fouts:
            fout.fin.close()
        
    
if __name__ == '__main__':
    _main()
