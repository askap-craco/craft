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
parser.add_argument('-r', '-refant', '--refant', default=1, help="Referance antenna", type=int)
parser.add_argument('aipsfile', help="AIPS data to fring")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

################################################################################
# Some useful functions
################################################################################
def run_fring(aipsdata, refant=1):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = aipsdata
    fring.refant = refant
    fring.solint = 0.5
    fring.dparm[9] = 1
    fring()

(name, aipsclass, seq) = args.aipsfile.split('.')

aipsfile = AIPSUVData(name, aipsclass, aipsDisk, int(seq))

if not aipsfile.exists():
    print args.aipsfile, "does not exists! Aborting"
    sys.exit()

run_fring(aipsfile, args.refant)
