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
import sys, os, math
import json

################################################################################
# Set up AIPS stuff
################################################################################
try:
    aipsver = os.environ['PSRPIAIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata, sources):
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 0
    fitld.ncount = 1
    if sources != "":
        fitld.sources[1:len(sources)+1] = sources
    fitld.outdata = aipsdata
    if aipsdata.exists():
        if interaction.yesno("Delete existing UV dataset " + aipsdata.name + \
                             "? (No will abort pipeline)"):
            aipsdata.zap()
        else:
            sys.exit()
    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits
    fitld()
    os.system("rm -f " + tempinfits)

################################################################################
# Main code
################################################################################

# Check invocation
if not len(sys.argv) == 3:
    print "Usage: %s <fits file> <json delay/phase file>" % sys.argv[0]
    sys.exit()

snversion = 1
inputfitsfile = sys.argv[1]
outputfile = sys.argv[2]

uvdata = AIPSUVData("CRAFTPOS", "UVSRT", 1, 3)
#uvdata = AIPSUVData("JUNKSN", "JUNKSN", 1, 1)
#if uvdata.exists():
#    uvdata.zap()

# Load up the FITS file into AIPS
#fitld_uvfits(inputfitsfile, uvdata, [])

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
phases1 = {}
delays2 = {}
phases2 = {}
sntable = uvdata.table('SN', snversion)
num_if = sntable.keywords['NO_IF']
num_pol = sntable.keywords['NO_POL']
num_snterms = num_if*num_pol
for row in sntable:
    if not annames[row.antenna_no-1] in delays1.keys():
        delays1[annames[row.antenna_no-1]] = {}
        phases1[annames[row.antenna_no-1]] = {}
        if num_pol > 1:
            delays2[annames[row.antenna_no-1]] = {}
            phases2[annames[row.antenna_no-1]] = {}
    delays1[annames[row.antenna_no-1]][str(row.time)] = []
    phases1[annames[row.antenna_no-1]][str(row.time)] = []
    if num_pol > 1:
        delays1[annames[row.antenna_no-1]][str(row.time)] = []
        phases1[annames[row.antenna_no-1]][str(row.time)] = []
    for i in range(num_if):
        delays1[annames[row.antenna_no-1]][str(row.time)].append(row.delay_1[i])
        phases1[annames[row.antenna_no-1]][str(row.time)].append(math.atan2(row.real1[i], row.imag1[i]))
        if num_pol > 1:
            delays2[annames[row.antenna_no-1]][str(row.time)].append(row.delay_2[i])
            phases2[annames[row.antenna_no-1]][str(row.time)].append(math.atan2(row.real2[i], row.imag2[i]))

with open(outputfile, 'w') as output:
     output.write(json.dumps(delays1,sort_keys=True) + "\n") # use `json.loads` to do the reverse
     output.write(json.dumps(phases1,sort_keys=True) + "\n") # use `json.loads` to do the reverse
     output.write(json.dumps(delays2,sort_keys=True) + "\n") # use `json.loads` to do the reverse
     output.write(json.dumps(phases2,sort_keys=True) + "\n") # use `json.loads` to do the reverse

print "Done"
