#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import numpy as np
import logging
import os
import astropy.io.fits as pyfits
from astropy.io.fits import Column as Col
from astropy.time import Time
from astropy import units as u
import warnings
import sys
import datetime
import scipy

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

parnames = ('UU','VV','WW','DATE','BASELINE','FREQSEL','SOURCE','INTTIM','DATA')
dtype_to_bitpix = {np.dtype('>i2'):16,
                   np.dtype('>i4'):32,
                   np.dtype('>f4'):-32,
                   np.dtype('>f8'):-64}

class CorrUvFitsFile(object):
    def __init__(self, fname, fcent, foff, nchan, npol, mjd0, sources, antennas, sideband=1, telescop='ASKAP', instrume='VCRAFT', origin='CRAFT', output_dtype=np.dtype('>f4'), bmax=None, time_scale=1.0*u.day):
        '''
        Make a correlator UV fits file
        :fname: file name
        :fcent: Center frequency (Hz)
        :foff: Channel offset (Hz):
        :nchan: numberof channels (Hz)
        :npol: Number of polarisations
        :mjd0: MJD of first sample
        :sources: List of osurces (format?)
        :antenas: List of antenas
        :sideband: (Not sure)
        :telescop: Put into header
        :instrume: put into header
        :origin: put into header
        :output_dtype: Dtype of output. According to fits standard: https://archive.stsci.edu/fits/fits_standard/node39.html - valid values are (np.int16, np.int32, np.float32, np.float64)
        :bmax: Astropy unit of distance that is the longest baseline in teh array. If set, UVW values are scale for dtype=int16 and dtype=int32 to get decent dynamic range
        :time_scale: Set to an astropy unit that converts to days, it will be used to specify the scale in FITS headers for the DATE and INTTIM columns. Can be used so the version on disk is in units of "integrations". THis might not be very sensible. If you set this to the integration time, you can use put_data(..t=integration_number) and get sensile output
        '''

        self.bmax = bmax  # maximum baseline 

        if bmax is None:
            uvw_scale = 1.0
        else:
            tmax = bmax.to(u.meter).value/scipy.constants.c # maximum baseline - seconds
            print('tmax', tmax, 'bmax', bmax, 'bmax to meter', bmax.to(u.meter))
            if output_dtype == np.dtype('int16'):
                vmax = 1<<15 # largerst value of a 16 bit signed number, if that's what we want to use for uvw
                uvw_scale = vmax/tmax
            elif output_dtype == np.dtype('int32'):
                vmax = 1<<31
                uvw_scale = vmax/tmax
            else:
                uvw_scale = 1.0
            
        self.uvw_scale = uvw_scale
        
        self.dshape = [1,1,1,nchan, npol, 3]
        hdr = pyfits.Header()
        self.hdr = hdr
        self.nchan = nchan
        self.foff = foff
        self.fcent = fcent
        self.npol = npol
        self.mjd0 = mjd0
        self.jd0 = self.mjd0 + 2400000.5
        self.bandwidth = foff*nchan
        self.sources = sources
        self.antennas = antennas
        self.fname = fname
        self.no_if = 1
        self.sideband = sideband

        bitpix = dtype_to_bitpix.get(output_dtype, None)
        if bitpix is None:
            raise ValueError(f'Invalid output_dtype={output_dtype}. Must be one of {dtype_to_bitpix.keys()}')

        print(f'Bitpix {bitpix}')
        self.output_dtype = output_dtype
        self.bitpix = bitpix

        hdr['SIMPLE'] = True
        hdr['BITPIX'] = bitpix
        self.set_axes()
        hdr['EXTEND'] = (True, 'Tables will follow')
        hdr['BLOCKED'] = (True, 'File may be blocked')
        hdr['GROUPS'] = (True, 'Random group UV data')
        hdr['PCOUNT'] = (8, 'Number of random parameters')
        hdr['GCOUNT'] = (1, 'Number of groups (rows) - updated when file closed')
        hdr['EPOCH'] = 2e3
        hdr['BSCALE'] = 1.0
        hdr['BZERO'] = 0.0
        hdr['BUNIT'] = 'UNCALIB'
        refchan = float(nchan)/2. + 0.5 # half a channel because it's centered on DC 
        # CTYPES
        self.add_type(2, ctype='COMPLEX', crval=1.0, cdelt=1.0, crpix=1.0)
        self.add_type(3, ctype='STOKES', crval=-5.0, cdelt=-1.0, crpix=1.0)
        self.add_type(4, ctype='FREQ', crval=fcent*1e6, cdelt=foff*1e6, crpix=refchan, crota=0.0)
        self.add_type(5, ctype='IF', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)
        self.add_type(6, ctype='RA', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)
        self.add_type(7, ctype='DEC', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)

        self.time_scale = time_scale
        tscale = 1./time_scale.to(u.day).value
        uvwscale = (1/self.uvw_scale)
        log.info("UVW Scale: %s header scale %s", self.uvw_scale, uvwscale)
        log.info("Time scale: %s header %s", self.time_scale, tscale)

        #ptypes
        self.add_type(1, ptype='UU', pscal=uvwscale, pzero=0.0)
        self.add_type(2, ptype='VV', pscal=uvwscale, pzero=0.0)
        self.add_type(3, ptype='WW', pscal=uvwscale, pzero=0.0)
        self.add_type(4, ptype='DATE', pscal=tscale, pzero=self.jd0, comment='Day number')
        #self.add_type(5, ptype='DATE', pscal=1.0, pzero=0.0, comment='Day fraction')
        self.add_type(5, ptype='BASELINE', pscal=1.0, pzero=0.0)
        self.add_type(6, ptype='FREQSEL', pscal=1.0, pzero=0.0)
        self.add_type(7, ptype='SOURCE', pscal=1.0, pzero=0.0)
        self.add_type(8, ptype='INTTIM', pscal=tscale, pzero=0.0)

        #self.first_time.format = 'fits'
        hdr['OBJECT'] = 'MULTI'
        #hdr['DATE_OBS'] = self.first_time.value Miriad gets confused by this 
        hdr['TELESCOP'] = telescop
        hdr['INSTRUME'] = instrume
        hdr['ORIGIN'] = origin
        hdr['OBSERVER'] = ''
        hdr['SORTORD'] = 'TB'
        hdr['SPECSYS'] = 'TOPOCENT'

        histstr = 'Created on {} by {}'.format(datetime.datetime.now().isoformat(), ' '.join(sys.argv))
        hdr['HISTORY'] = histstr

        self.fout = open(fname, 'w+b')
        self.fout.write(bytes(hdr.tostring(), 'utf-8'))
        self.ngroups = 0
        dt = self.output_dtype
        assert dt.byteorder == '>', 'Byte order of FITS output ata type must be big endian'
        # aaah, craparooney - dthe data type for the whole row has to be the same. This means you can't
        # just have 165 bit data and 32 bit UVWs, which means it's rubbisharooney, unless I can
        # be bothered ot do somethign with BZERO and BSCALE (Maybe?)
        self.dtype = np.dtype([('UU', dt), ('VV', dt), ('WW', dt), \
            ('DATE', dt), ('BASELINE', dt), \
            ('FREQSEL', dt), ('SOURCE', dt), ('INTTIM', dt), \
            ('DATA', dt, (1, 1, 1, nchan, npol, 3))])

        log.debug('Dtype size is %s', self.dtype.itemsize)

    def fq_table(self):
        cols = [Col('FRQSEL','1J', array=[1]),
                Col('IF FREQ','1D','HZ', array=[0]), #??
                Col('CH WIDTH','1E','HZ', array=[self.foff*1e6]),
                Col('TOTAL BANDWIDTH','1E','HZ', array=[self.bandwidth*1e6]),
                Col('SIDEBAND','1J', array=[self.sideband]), #??
                Col('RXCODE','1J',array=[0]) # ?? Not in AIPSMEM117, but AIPS/FRING1 complains if it isn't there
        ]
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.header['EXTNAME'] = 'AIPS FQ'
        tbhdu.header['NO_IF'] = (self.no_if, 'hard coded')

        return tbhdu

    def an_table(self, antennas):
        nant = len(antennas)

        cols = [
            Col('ANNAME', '8A', array=[a.antname for a in antennas]),
            Col('STABXYZ','3D','METERS', array=[a.antpos for a in antennas]),
            Col('ORBPARM','1D', array=[]),
            Col('NOSTA','1J', array=np.arange(1,nant+1, dtype=np.int32)),
            Col('MNTSTA','1J', array=np.ones(nant, dtype=np.int32)),
            Col('STAXOF','1E','METERS', array=np.zeros(nant)),
            Col('POLTYA','1A', array=['X' for i in range(nant)]),
            Col('POLAA','1E','DEGREES', array=np.zeros(nant)),
            Col('POLCALA','1E', array=np.zeros(nant)),
            Col('POLTYB','1A', array=['Y' for i in range(nant)]),
            Col('POLAB','1E','DEGREES', array=np.zeros(nant)),
            Col('POLCALB','1E', array=np.zeros(nant)),
            Col('DIAMETER','1E', array=np.ones(nant)*12.)
        ]
        h = pyfits.BinTableHDU.from_columns(cols)
        hdr = h.header
        hdr['EXTNAME']= 'AIPS AN'
        hdr['ARRAYX'] = 0
        hdr['ARRAYY'] = 0
        hdr['ARRAYZ'] = 0
        hdr['GSTIA0'] = (1.764940880935E+02 , 'hard coded. ????')
        hdr['DEGPDY'] = (3.609856473692E+02, 'hard coded. ????')
        hdr['FREQ'] = self.fcent*1e6
        #hdr['RDATE'] = (self.first_time.value, 'hard coded. ????') # miriad gets confused by this
        hdr['POLARX'] = 0
        hdr['POLARY'] = 0
        hdr['UT1UTC'] = 0
        hdr['IATUTC'] = (3.700000000000E+01, 'hard coded')
        hdr['TIMSYS'] = 'UTC'
        hdr['ARRNAM'] = 'ASKAP'
        hdr['NUMORB'] = 0
        hdr['NOPCAL'] = 0
        hdr['NO_IF'] = (self.no_if, 'hard coded')
        hdr['POLTYPE'] = ''
        hdr['XYZHAND'] = 'RIGHT'
        hdr['FRAME'] = 'ITRF'

        return h

    def su_table(self, sources):
        ns = len(sources)
        zeros = np.zeros(ns)
        cols = [
            Col('ID. NO.', '1J', array=np.arange(1,ns+1, dtype=np.int32)),
            Col('SOURCE','50A','METERS', array=[s['name'] for s in sources]),
            Col('QUAL','1J', array=zeros),
            Col('CALCODE','4A', array=None),
            Col('IFLUX','1E','JY', array=zeros),
            Col('QFLUX','1E','JY', array=zeros),
            Col('UFLUX','1E','JY', array=zeros),
            Col('VFLUX','1E','JY', array=zeros),
            Col('FREQOFF','1D', array=zeros),
            Col('BANDWIDTH','1D','HZ',array=np.ones(ns)*self.bandwidth),
            Col('RAEPO','1D','DEGREES', array=[s['ra'] for s in sources]),
            Col('DECEPO','1D','YEARS', array=[s['dec'] for s in sources]),
            Col('EPOCH','1D','DEGREES', array=np.ones(ns)*2000.),
            Col('RAAPP','1D','DEGREES', array=[s['ra'] for s in sources]),
            Col('DECAPP','1D','DEGREES', array=[s['dec'] for s in sources]),
            Col('LSRVEL','1D','M/SEC', array=zeros),
            Col('RESTFREQ','1D','HZ', array=zeros),
            Col('PMRA', '1D', 'DEG/DAY', array=zeros),
            Col('PMDEC','1D','DEG/DAY', array=zeros)
        ]
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.header['EXTNAME'] = 'AIPS SU'
        tbhdu.header['NO_IF'] = self.no_if
        tbhdu.header['FREQID'] = 1

        return tbhdu

    def close(self):
        if self.fout is None:
            return
        
        fout = self.fout
        hdr = self.hdr
        currbytes = fout.tell()
        n_extra_bytes = 2880 - currbytes % 2880
        if n_extra_bytes == 2880:
            n_extra_bytes = 0

        fout.write(bytes(n_extra_bytes))

        # update headdr
        fout.seek(0, 0)
        hdr['GCOUNT'] = self.ngroups
        fout.write(bytes(hdr.tostring(), 'utf-8'))
        assert fout.tell() % 2880 == 0
        fout.flush()
        fout.close()

        hdu = pyfits.open(self.fname, 'append')
        fq = self.fq_table()

        hdu.append(fq)

        an = self.an_table(self.antennas)
        hdu.append(an)

        su = self.su_table(self.sources)
        hdu.append(su)
        hdu.close()
        self.fout = None

    def add_type(self, typeno, comment=None, **args):
        hdr = self.hdr
        for k,v, in args.items():
            hdr['{}{}'.format(k,typeno).upper()] = (v, comment)

    def put_data(self, uvw, mjd ,ia1, ia2, inttim, data, weights=None, source=1, t=None):
        '''
        :uvw: UVW in meters
        :mjd: MJD of integration - you must have time_scale=1.0day
        :t: If Not none, this is put into the file raw for the date, instead of mjd0.
        :source: 1-indexed source id
        '''
        assert ia1 >= 0
        assert ia2 >= 0
        assert 1<= source <= len(self.sources), f'Invaoid source ID={source}. Source table has {len(self.sources)}'

        visdata_all = np.recarray(1, dtype=self.dtype)
        visdata = visdata_all[0]
        if weights is None:
            weights = 7.71604973e-05 #???

        jd = mjd + 2400000.5
        day = np.floor(jd)
        dayfrac = jd - day

        visdata['UU'], visdata['VV'], visdata['WW'] = uvw*self.uvw_scale
        #visdata['DATE'] = jd - self.jd0
        if t is None:
            assert self.time_scale == 1*u.day
            visdata['DATE'] = mjd - self.mjd0
        else:
            visdata['DATE'] = t
        #visdata['_DATE'] = dayfrac
        visdata['BASELINE'] = (ia1 + 1)*256 + ia2 + 1
        visdata['INTTIM'] = inttim
        visdata['FREQSEL'] = 1
        visdata['SOURCE'] = source

        d = visdata['DATA']
        npol = data.shape[1]
        if np.iscomplexobj(data):
            d[0,0,0,:,:,0] = data.real
            d[0,0,0,:,:,1] = data.imag
        else:
            d[0,0,0,:,:,0] = data[...,0]
            d[0,0,0,:,:,1] = data[...,1]
            
        d[0,0,0,:,:,2] = weights
        self.fout.write(visdata_all.tobytes())
        self.ngroups += 1


    def set_axes(self):
        hdr = self.hdr
        hdr['NAXIS'] = len(self.dshape) + 1
        hdr['NAXIS1'] = (0, 'Random axis group')
        for ia, a in enumerate(self.dshape[::-1]):
            hdr['NAXIS{}'.format(ia+2)] = a


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

    fname= values.files[0]
    f = CorrUvFitsFile(fname, 850e6, 18e3)
    f.put_data()
    f.close()
    fin = pyfits.open(fname)
    fin.info()
    print(fin[0].data[0].data.shape)


if __name__ == '__main__':
    _main()
