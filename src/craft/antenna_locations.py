#!/usr/bin/env python
"""
List of antenna locations that are pre-installed with the package

Copyright (C) CSIRO 2020
"""
import os
import glob
import logging


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

'''
Return full directory location of pre-installed antenna locations files
'''
locations_path = os.path.realpath(os.path.join(os.path.dirname(__file__), 'data','antenna_locations'))


'''
Returns a list of paths of all antenna locations
'''
antenna_files = glob.glob(os.path.join(locations_path, '*.ant'))

def get_antenna_file(filename):
    '''
    Tries to find an antenna file name. Returns the path to the antenan file.
    Will try to find the file first as specifed on the commad line. If it isn't there,
    it will try to find the path inside the pre-installed area
    
    @param filename - name of the antenna file
    @returns path to the file name
    @throws FileNotFoundError if that filename couldn't be found in the local filesystem, or the pre-installed area
    '''
    if os.path.isfile(filename):
        fullpath = os.path.realpath(filename)
    else:
        fullpath = os.path.join(locations_path, filename)

    if not os.path.isfile(fullpath):
        raise FileNotFoundError(f'Could not find antenna file {filename} checked {locations_path}')

    return fullpath

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='*', help='If you dont specify, itll list all pre-installed ones. Otherwise, itl print the path to the one  you want')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.verbose:
        print(f'Antenna locations directory:{locations_path}')

    if len(values.files) == 0:
        for p in antenna_files:
            print(p)
    else:
        for p in values.files:
            afile = get_antenna_file(p)
            if values.verbose:
                print(f'Path {p} resolves to {afile}')
            else:
                print(afile)
            
    

if __name__ == '__main__':
    _main()
