#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import numpy as np
import logging
import os
import pyfits
from pyfits import Column as Col

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

parnames = ('UU','VV','WW','DATE','_DATE','BASELINE','FREQSEL','SOURCE','INTTIM','DATA')

class CorrUvFitsFile(object):
    def __init__(self, fname, fcent, foff, nchan, npol, sources, antennas):
        self.dshape = [1,1,1,nchan, npol, 3]
        hdr = pyfits.Header()
        self.hdr = hdr
        self.nchan = nchan
        self.foff = foff
        self.fcent = fcent
        self.npol = npol
        self.bandwidth = foff*nchan
        self.sources = sources
        self.antennas = antennas
        self.fname = fname

        hdr['SIMPLE'] = True
        hdr['BITPIX'] = -32
        self.set_axes()
        hdr['EXTEND'] = (True, 'Tables will follow')
        hdr['BLOCKED'] = (True, 'File may be blocked')
        hdr['GROUPS'] = (True, 'Random group UV data')
        hdr['PCOUNT'] = (9, 'Number of random parameters')
        hdr['GCOUNT'] = (1, 'Number of groups (rows) - updated when file closed')
        hdr['EPOCH'] = 2e3
        hdr['BSCALE'] = 1.0
        hdr['BZERO'] = 0.0
        hdr['BUNIT'] = 'UNCALIB'
        # CTYPES
        self.add_type(2, ctype='COMPLEX', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)
        self.add_type(3, ctype='STOKES', crval=-5.0, cdelt=-1.0, crpix=1.0, crota=0.0)
        self.add_type(4, ctype='FREQ', crval=fcent*1e6, cdelt=foff*1e6, crpix=(float(nchan)/2. + 1), crota=0.0)
        self.add_type(5, ctype='IF', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)
        self.add_type(6, ctype='RA', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)
        self.add_type(7, ctype='DEC', crval=1.0, cdelt=1.0, crpix=1.0, crota=0.0)

        #ptypes
        self.add_type(1, ptype='UU', pscal=1.0, pzero=0.0)
        self.add_type(2, ptype='VV', pscal=1.0, pzero=0.0)
        self.add_type(3, ptype='WW', pscal=1.0, pzero=0.0)
        self.add_type(4, ptype='DATE', pscal=1.0, pzero=0.0, comment='Day number')
        self.add_type(5, ptype='DATE', pscal=1.0, pzero=0.0, comment='Day fraction')
        self.add_type(6, ptype='BASELINE', pscal=1.0, pzero=0.0)
        self.add_type(7, ptype='FREQSEL', pscal=1.0, pzero=0.0)
        self.add_type(8, ptype='SOURCE', pscal=1.0, pzero=0.0)
        self.add_type(9, ptype='INTTIM', pscal=1.0, pzero=0.0)

        hdr['OBJECT'] = 'MULTI'
        hdr['DATE_OBS'] = '2018-03-19T11:02:30.909121'
        hdr['TELESCOP'] = 'ASKAP'
        hdr['INSTRUME'] = 'VCRAFT'
        hdr['OBSERVER'] = ''
        hdr['SORTORD'] = 'TB'
        hdr['SPECSYS'] = 'TOPOCENT'
        hdr['ORIGIN'] = 'craftcor'

        self.fout = open(fname, 'w+')
        self.fout.write(hdr.tostring())
        self.ngroups = 0
        self.dtype = [('UU', '>f4'), ('VV', '>f4'), ('WW', '>f4'), \
            ('DATE', '>f4'), ('_DATE', '>f4'), ('BASELINE', '>f4'), \
            ('FREQSEL', '>f4'), ('SOURCE', '>f4'), ('INTTIM', '>f4'), \
            ('DATA', '>f4', (1, 1, 1, nchan, npol, 3))]

    def fq_table(self):
        cols = [Col('FRQSEL','1J', array=[1]),
                Col('IF FREQ','1D','HZ', array=[0]), #??
                Col('CH WIDTH','1E','HZ', array=[self.foff*1e6]),
                Col('TOTAL BANDWIDTH','1E','HZ', array=[self.bandwidth*1e6]),
                Col('SIDEBAND','1J', array=[1]) #??
        ]

        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.header['EXTNAME'] = 'AIPS FQ'

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
            Col('POLTYA','1A', array=['X' for i in xrange(nant)]),
            Col('POLAA','1E','DEGREES', array=np.zeros(nant)),
            Col('POLCALA','1E', array=np.zeros(nant)),
            Col('POLTYB','1A', array=['Y' for i in xrange(nant)]),
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
        hdr['RDATE'] = ('2018-03-19T11:02:30.909121', 'hard coded. ????')
        hdr['POLARX'] = 0
        hdr['POLARY'] = 0
        hdr['UT1UTC'] = 0
        hdr['IATUTC'] = (3.700000000000E+01, 'hard coded')
        hdr['TIMSYS'] = 'UTC'
        hdr['ARRNAM'] = 'ASKAP'
        hdr['NUMORB'] = 0
        hdr['NOPCAL'] = 0
        hdr['POLTYPE'] = ''
        hdr['XYZHAND'] = 'RIGHT'
        hdr['FRAME'] = 'ITRF'

        return h

    def su_table(self, sources):
        ns = len(sources)
        zeros = np.zeros(ns)
        cols = [
            Col('ID. NO.', '1J', array=np.arange(1,ns+1, dtype=np.int32)),
            Col('SOURCE','20A','METERS', array=[s['name'] for s in sources]),
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
        tbhdu.header['NO_IF'] = 1
        tbhdu.header['FREQID'] = 1

        return tbhdu

    def close(self):
        fout = self.fout
        hdr = self.hdr
        currbytes = fout.tell()
        n_extra_bytes = 2880 - currbytes % 2880
        if n_extra_bytes == 2880:
            n_extra_bytes = 0

        fout.write('\x00'*n_extra_bytes)

        # update headdr
        fout.seek(0, 0)
        hdr['GCOUNT'] = self.ngroups
        fout.write(hdr.tostring())
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

    def add_type(self, typeno, comment=None, **args):
        hdr = self.hdr
        for k,v, in args.iteritems():
            hdr['{}{}'.format(k,typeno).upper()] = (v, comment)

    def put_data(self, uvw, mjd, ia1, ia2, inttim, data, weights=None):
        assert ia1 >= 0
        assert ia2 >= 0
        visdata_all = np.recarray(1, dtype=self.dtype)
        visdata = visdata_all[0]
        if weights is None:
            weights = 7.71604973e-05 #???


        jd = mjd + 2400000.5
        day = np.floor(jd)
        dayfrac = jd - day

        visdata['UU'], visdata['VV'], visdata['WW'] = uvw
        visdata['DATE'] = day
        visdata['_DATE'] = dayfrac
        visdata['BASELINE'] = (ia1 + 1)*256 + ia2 + 1
        visdata['INTTIM'] = inttim
        visdata['FREQSEL'] = 1
        visdata['SOURCE'] = 1

        d = visdata['DATA']
        npol = data.shape[1]
        d[0,0,0,:,:,0] = data.real
        d[0,0,0,:,:,1] = data.imag
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
    print fin[0].data[0].data.shape


if __name__ == '__main__':
    _main()
