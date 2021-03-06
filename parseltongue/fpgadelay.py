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
    fring.dparm[9] = 1
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
delays1 = {}
delays2 = {}
delaysN = {}

sntable = uvdata.table('SN', 1)
num_if = sntable.keywords['NO_IF']
num_pol = sntable.keywords['NO_POL']
num_snterms = num_if*num_pol

fpga = os.path.basename(os.getcwd())

for row in sntable:
    ant = annames[row.antenna_no-1]
    if not ant in delays1:
        delays1[ant] = [0.0]*num_if
        if num_pol>1: delays2[ant] = [0.0]*num_if
        delaysN[ant] = 0
            
    for i in range(num_if):
        #if abs(row.delay_1[i])>1 or (num_pol > 1 and abs(row.delay_2[i])>1): continue

        delays1[ant][i] += row.delay_1[i]
        if num_pol>1: delays2[ant][i] += row.delay_2[i]
    delaysN[ant] += 1

if avgIf:
    for ant in delays1:
        #print ant
        delaysN[ant] *= num_if
        for i in range(1,num_if):
            delays1[ant][0] += delays1[ant][i]
        if num_pol>1:
            delaysN[ant] *= 2
            for i in range(num_if):
                delays1[ant][0] += delays2[ant][i]
        #print delays1[ant][0]

    for ant in sorted(delays1):
        if args.outfile is not None:
            delayfile =  open(args.outfile, 'a+')
        if delaysN[ant]>0:
            FPGAdelay = delays1[ant][0]/delaysN[ant]*1e9/6750.0
            if abs(FPGAdelay)>0.1:
                print "{} {:3.0f}".format(ant, FPGAdelay)
                if args.outfile is not None:
                    delayfile.write("{} {} {:3.0f}\n".format(ant, fpga, FPGAdelay))
        if args.outfile is not None:
            delayfile.close()

else:
    for ant in sorted(delays1):
        if delaysN[ant]==0:
            print ant, 0
        else:
            print ant,
            for i in range(num_if):
                print "{:4.1f}".format(delays1[ant][i]/delaysN[ant]*1e9/6750.0),",",
                if num_pol>1:
                    print "{:4.1f}".format(delays2[ant][i]/delaysN[ant]*1e9/6750.0),",",
            print
