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
parser.add_argument('-o', '--outfile', default="data", help="Output AIPS file")
parser.add_argument('-f', '--fitsfileoutname', help="Output FITS file, default doesn't write")
parser.add_argument('-a', '--antlist', help="Force antenna list")
#parser.add_argument('-s', '--specav', default=18, help="Spectral Averaging", type=int)
#parser.add_argument("-m", "--nomerge", default=False, action="store_true", help="Don't Merge IFs")
parser.add_argument("-d", "--delete", default=False, action="store_true", help="Delete AIPS file if writing FITS out")
parser.add_argument('fitsfile', nargs='+', help="Input FITS data")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

fitsfileoutname = args.fitsfileoutname
if fitsfileoutname is not None:
    if not fitsfileoutname[0] == '/':
        fitsfileoutname = os.getcwd() + '/' + args.fitsfileoutname
        if os.path.exists(fitsfileoutname):
            print "Warning: {} already exists. Aborting".format(fitsfileoutname)
            sys.exit()

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata, antlist=None):
    # Load FITS file

    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()

    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 1
    fitld.ncount = 1
    fitld.doconcat = 1
    fitld.outdata = aipsdata
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits

    if antlist is not None:
        antname = [x.strip() for x in antlist.split(',')]
        fitld.antname = AIPSList(antname)
    
    fitld()
    os.system("rm -f " + tempinfits)

################################################################################
# Main code
################################################################################

card_uvdata = []

single = len(args.fitsfile)==1

seq = 1
while True:
    outdata = AIPSUVData(args.outfile, "UVDATA", aipsDisk, seq)
    if not outdata.exists(): break
    seq += 1


for fits in args.fitsfile:

    if not os.path.exists(fits):
        print "{} does not exist. Aborting".format(fits)
        sys.exit()
    
    # Load up the FITS file into AIPS
    fitld_uvfits(fits, outdata, args.antlist)

if fitsfileoutname is not None:
    if not fitsfileoutname[0] == '/':
        fitsfileoutname = os.getcwd() + '/' + fitsfileoutname
    fittp = AIPSTask("fittp")
    fittp.indata = outdata
    fittp.dataout = fitsfileoutname
    fittp()
    if args.delete:
        outdata.zap()
