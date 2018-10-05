#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
#from AIPS import AIPS, AIPSDisk
from AIPS import AIPS
#from AIPSTask import AIPSTask, AIPSList
#from AIPSTask import AIPSTask
#from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSData import AIPSUVData
#from AIPSTV import AIPSTV

################################################################################
# General python imports
################################################################################
import argparse, sys, os, math
#import json

################################################################################
# Set up AIPS stuff
################################################################################
try:
    aipsver = os.environ['PSRPIAIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702
snversion = 1

parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', help="AIPS user number", type=int)
parser.add_argument('-s', '--sn', help="SN table version", type=int)
parser.add_argument('-o', '--outfile', help="Save delays to file")
parser.add_argument("-a", "--av", default=False, action="store_true", help="Average IFs")
parser.add_argument('aipsfile', help="AIPS  file ")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user
if args.sn is not None: snversion = args.user
avgIf = args.av

aipsDisk = 1

aipsfile = args.aipsfile

(name, aipsclass, seq) = aipsfile.split('.')

uvdata = AIPSUVData(name, aipsclass, aipsDisk, int(seq))

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

sntable = uvdata.table('SN', snversion)
num_if = sntable.keywords['NO_IF']
num_pol = sntable.keywords['NO_POL']
num_snterms = num_if*num_pol

for row in sntable:
    ant = annames[row.antenna_no-1]
    if not ant in delays1:
        delays1[ant] = [0.0]*num_if
        delays2[ant] = [0.0]*num_if
        delaysN[ant] = 0

    if num_if==1:
        rowDelay1 = [row.delay_1]
        if num_pol>1:
            rowDelay2 = [row.delay_2]
    else:
        rowDelay1 = row.delay_1
        if num_pol>1:
            rowDelay2 = row.delay_2
    for i in range(num_if):
        if abs(rowDelay1[i])>1 or (num_pol > 1 and abs(rowDelay2[i])>1): continue
        delays1[ant][i] += rowDelay1[i]
        if num_pol > 1:
            delays2[ant][i] += rowDelay2[i]
    delaysN[ant] += 1

if avgIf:
    for ant in delays1:
        for i in range(1,num_if):
            delays1[ant][0] += delays1[ant][i]
        delaysN[ant] *= num_if

    for ant in sorted(delays1):
        if delaysN[ant]==0:
            print ant, 0
        else:
            print ant, delays1[ant][0]/delaysN[ant]*1e9

else:

    for ant in sorted(delays1):
        if delaysN[ant]==0:
            print ant, 0
        else:
            print ant,
            for i in range(num_if):
                print delays1[ant][i]/delaysN[ant]*1e9,",",
            print
