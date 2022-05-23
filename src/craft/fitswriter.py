#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import numpy as np
import os
import sys
import logging
import astropy.io.fits as fits

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class FitsTableWriter:
    '''
    A generic class that writes infinitely-long FITS tables efficiently and quickly

    It works by using astropy.io.fits to generate a dummy file with headers
    Then it writes blocks of data ta the end of the file.
    When the file is closed, it updates the header to include the number
    of blocks written.
    It's a bit nasty, but it works.

    Note: FITS is insane - it requires big-endian data. To be compatible with fits
    set byteswap=True (the default). If you want maximum efficience, set byteswap=False
    and byteswap when you read it in

    The header card 'BYTESWAP' will be 'TRUE' if byteswap was True, otherwise 'FALSE'
    

    @see corruvfits
    '''
    def __init__(self, fname, dtype, byteswap=True, header={}):
        self.fname = fname
        self.dtype = dtype
        self.byteswap = byteswap

        # create a dummy table with 1 element and write it to disk
        # this generates all teh headers that we'll need for later
        # without us haveing to write too much code
        single_element = np.ones(1, dtype=dtype)
        hdu = fits.BinTableHDU(data=single_element)
        hdu.header['BYTESWAP'] = 'TRUE' if byteswap  else 'FALSE'

        for k,v in header.items():
            hdu.header[k] = v
        
        for c in hdu.columns:
            c.bzero = 0 # astropy.io.fits is ***ING* insane and sets bzero to some insane value because it's insane. That's INSANE
            #pass

        hdu.writeto(fname, overwrite=True)

        file_size = os.path.getsize(fname)

        # now we open the file a a normal python file and seek to where the table starts
        # WE'll get the headers for the primary extension and the table header
        self.prihdr = fits.getheader(fname)
        self.tabhdr = fits.getheader(fname, 1)
        header_length = len(self.prihdr.tostring()) + len(self.tabhdr.tostring())
        self.fout = open(fname, 'r+b')
        self.fout.seek(header_length)

        self.ngroups = 0

        # now we're ready to push data into the file

    def write(self, data):
        '''
        Writes the given data to the file
        Flattens it first - but it can have any shape
        '''
        assert data.dtype == self.dtype

        # FITS file format must be big-endian. There's no way to chagne this. IN ... SANE ...
        if self.byteswap:
            dout = data.copy().byteswap()
        else:
            dout = data

        self.fout.write(dout.tobytes())
        self.ngroups += len(data)

    def close(self):
        '''
        Closes the file - writes zeros to make the file a multipl eof 2880 bytes
        as required by the FITS standard and hten updates the header accordingly
        '''

        fout = self.fout
        # Write zeros to make the whole thing a multiple of 2880 bytes
        currbytes = fout.tell()
        n_extra_bytes = 2880 - currbytes % 2880
        if n_extra_bytes == 2880:
            n_extra_bytes = 0

        fout.write(bytes(n_extra_bytes))
        assert fout.tell() % 2880 == 0
        file_size = fout.tell()
        #assert file_size == os.path.getsize(self.fname)

        # update extension header
        hdr = self.tabhdr
        hdr['NAXIS2'] = self.ngroups

        # seek to the start of the table extension header
        start = len(self.prihdr.tostring())
        fout.seek(start,0)
#        assert fout.tell() == start
        fout.write(hdr.tostring().encode('utf-8'))
        fout.close()
#        assert file_size == os.path.getsize(self.fname)


    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return True

        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
