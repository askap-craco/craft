#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask
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
    aipsver = os.environ['PSRPIAIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', help="AIPS user number", type=int)
parser.add_argument('-d', '--disk', default=1, help="AIPS disk", type=int)
parser.add_argument('fitsfile', help="Input FITS data")
parser.add_argument('aipsfile', help="Ouput AIPS data")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata):
    # Load FITS file using FITLD

    if not os.path.exists(uvfitsfile):
        print(uvfitsfile + " does not exist! Aborting.")
        sys.exit()
    infits = os.path.abspath(uvfitsfile)

    tmp = aipsdata.split('.')
    outname = tmp[0]
    if len(tmp)<2:
        AIPSclass = 'UVDATA'
    else:
        AIPSclass = tmp[1]

    if len(tmp)<3:
        seq = 1
        while True:
            dataout = AIPSUVData(outname, AIPSclass, aipsDisk, seq)
            if not dataout.exists(): break
            seq += 1
    else:
        seq = tmp[2]
        dataout = AIPSUVData(outname, AIPSclass, aipsDisk, seq)
        if not dataout.exists():
            print("Error: {} already exists".format(aipsdata))
            sys.exit()
            
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.datain = infits
    fitld.outdata = dataout
    fitld.digicor = 1
    fitld.ncount = 1
    fitld()

fitld_uvfits(args.fitsfile, args.aipsfile)
