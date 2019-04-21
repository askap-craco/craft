#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData

################################################################################
# General python imports
################################################################################
import argparse, sys, os, math

pid = os.getpid()

################################################################################
# Set up AIPS stuff
################################################################################
try:
    aipsver = os.environ['AIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', help="AIPS user number", type=int)
parser.add_argument('-D', '--disk', default=1, help="AIPS disk", type=int)
parser.add_argument('-b', '--bp', default=1, help="Bandpass table to apply, Default 1", type=int)
parser.add_argument('-f', '--fg', default=1, help="Flag table to apply, Default 1", type=int)
parser.add_argument('-C', '--cl', default=3, help="CL table to apply, Default 3", type=int)
parser.add_argument("-c", "--calib", default=False, action="store_true", help="Apply Calibration")
parser.add_argument('aipsfile', help="AIPS file to export")
parser.add_argument('fitsfile', help="UVFITS file to write")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

################################################################################
# Some useful functions
################################################################################
def run_fittp(aipsdata, fitsfile):
    # Write an AIPS file out as uvfits

    if not fitsfile[0] == '/':
        fitsfile = os.getcwd() + '/' + fitsfile
    
    if os.path.exists(fitsfile):
        print "Warning: {} already exists. Aborting".format(fitsfile)
        sys.exit()

    fittp = AIPSTask("fittp")
    fittp.indata = aipsdata
    fittp.dataout = fitsfile
    fittp()

def run_splat(aipsdata, outdata, clver=None, bpver=None, fgver=None):
    splat = AIPSTask('splat', version = aipsver)
    splat.indata = aipsdata
    splat.aparm[5] = 1
    if fgver is not None:
        splat.flagver = fgver
    if bpver is not None:
        splat.doband = 1
        splat.bpver = bpver
    if clver is not None:
        splat.docalib = 1
        splat.gainuse = clver
    splat.outdata = outdata
    splat()
    

################################################################################
# Main code
################################################################################


(name, aipsclass, seq) = args.aipsfile.split('.')
aipsdata = AIPSUVData(name, aipsclass, aipsDisk, int(seq))
if not aipsdata.exists():
    print args.aipsfile, "does not exists! Aborting"
    sys.exit()

if args.calib: # Need to apply calibration
    tempaips = AIPSUVData("TMPCAL{}".format(pid), "UVDATA", aipsDisk, 1)
    if tempaips.exists(): loadUV.zap()
    run_splat(aipsdata, tempaips, args.cl, args.bp, args.fg)
    run_fittp(tempaips, args.fitsfile)
    tempaips.zap()
else:
    run_fittp(aipsdata, args.fitsfile)

