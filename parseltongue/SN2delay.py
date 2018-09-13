#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
#from AIPSTV import AIPSTV

################################################################################
# General python imports
################################################################################
import sys, os, math
#import json

################################################################################
# Set up AIPS stuff
################################################################################
try:
    aipsver = os.environ['PSRPIAIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

# Check invocation
if not len(sys.argv) == 2:
    print "Usage: %s <AIPSFILE>" % sys.argv[0]
    sys.exit()

snversion = 1    # Should get from getopt style argument
aipsDisk = 1

aipsfile = sys.argv[1]

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
        delays1[ant] = 0.0
        delays2[ant] = 0.0
        delaysN[ant] = 0 
    
    for i in range(num_if):
        if abs(row.delay_1[i])>1 or (num_pol > 1 and abs(row.delay_1[i])>1): continue
        
        delays1[ant] += row.delay_1[i]
        if num_pol > 1:
            delays2[ant] += row.delay_2[i]
        delaysN[ant] += 1
       
for ant in delays1:
    if delaysN[ant]==0:
        print ant, 0
    else:
        print ant, delays1[ant]/delaysN[ant]*1e9
