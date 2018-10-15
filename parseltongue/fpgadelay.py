#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

################################################################################
# General python imports
################################################################################
import argparse, sys, os, math

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
parser.add_argument('-o', '--outfile', help="Save delays to file")
parser.add_argument('fitsfile', help="FITS file to load")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

avgIf = True

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata):
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 1
    fitld.ncount = 1
    fitld.outdata = aipsdata
    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits
    fitld()
    os.system("rm -f " + tempinfits)

def fring_fpga(aipsdata):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = aipsdata
    fring.refant = 1
    fring.solint = 0.5
    fring()

    

################################################################################
# Main code
################################################################################

snversion = 1
inputfitsfile = args.fitsfile

uvdata = AIPSUVData("TMPDELAY", "UVDATA", 1,1)

if uvdata.exists(): uvdata.zap()

# Load up the FITS file into AIPS
fitld_uvfits(inputfitsfile, uvdata)

fring_fpga(uvdata)

# Make a list of antennas
maxanid = 0
antable = uvdata.table('AN', 1)
for row in antable:
    if row.nosta > maxanid:
        maxanid = row.nosta
annames = []
for i in range(maxanid):
    annames.append("")
for row in antable:
    annames[row.nosta-1] = row.anname.strip()

# Read the AIPS delays and phases
delays = {}
delaysN = {}

sntable = uvdata.table('SN', 1)
num_if = sntable.keywords['NO_IF']
num_pol = sntable.keywords['NO_POL']
num_snterms = num_if*num_pol

fpga = os.path.basename(os.getcwd())

if num_pol>1:
    print "Only support single polarisation"
    exit()

for row in sntable:
    ant = annames[row.antenna_no-1]
    if not ant in delays:
        delays[ant] = [0.0]*num_if
        delaysN[ant] = 0
            
    for i in range(num_if):
        if abs(row.delay_1[i])>1 or (num_pol > 1 and abs(row.delay_1[i])>1): continue

        delays[ant][i] += row.delay_1[i]
    delaysN[ant] += 1

if avgIf:
    for ant in delays:
        for i in range(1,num_if):
            delays[ant][0] += delays[ant][i]
        delaysN[ant] *= num_if

    for ant in delays:
        if args.outfile is not None:
            delayfile =  open(args.outfile, 'a+')
        if delaysN[ant]>0:
            FPGAdelay = delays[ant][0]/delaysN[ant]*1e9/6750.0
            if abs(FPGAdelay>0.1):
                print "{} {:3.0f}".format(ant, FPGAdelay)
                if args.outfile is not None:
                    delayfile.write("{} {} {:3.0f} {:3.0f}\n".format(ant, fpga, FPGAdelay, FPGAdelay))
        if args.outfile is not None:
            delayfile.close()

else:

    for ant in delays:
        if delaysN[ant]==0:
            print ant, 0
        else:
            print ant,
            for i in range(num_if):
                print "{:4.1f}".format(delays[ant][i]/delaysN[ant]*1e9/6750.0),",",
            print
