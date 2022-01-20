#!/bin/env python3

import sys, os, glob

if len(sys.argv) != 3:
    print(("Usage: {} <SSID1> <SSID2>".format(sys.argv[0])))
    sys.exit()

SSID1path = sys.argv[1]
SSID2path = sys.argv[2]

SSID1 = os.path.basename(SSID1path)
SSID2 = os.path.basename(SSID2path)

outdir = "{}-{}".format(SSID1,SSID2)

# Check directories exists

if not os.path.exists(SSID1path):
    print(("{} does not exist! Exiting".format(SSID1)))
    sys.exit()
if not os.path.isdir(SSID1path):
    print(("{} is not a directory. Exiting".format(SSID1)))
    sys.exit()
if not os.path.exists(SSID2path):
    print(("{} does not exist! Exiting".format(SSID2)))
    sys.exit()
if not os.path.isdir(SSID2path):
    print(("{} is not a directory. Exiting".format(SSID2)))
    sys.exit()


ants1 = [os.path.basename(i) for i in glob.glob("{}/ak*".format(SSID1path))]
ants1.sort()
ants2 = [os.path.basename(i) for i in glob.glob("{}/ak*".format(SSID2path))]
ants2.sort()

if not ants1 == ants2:
    print(("Error: {} and {} do not contain the same antennas!!".format(SSID1, SSID2)))
    sys.exit()

# Check the antennas themselves are directories

for s in [SSID1path, SSID2path]:
    for a in ants1:
        if not os.path.isdir("{}/{}".format(s,a)):
            print(("{}/{} is not a directory. Exiting".format(s,a)))
            sys.exit()

# Create a new directory locally and add symlinks

try:
    os.mkdir(outdir)
except OSError:  # May could continue if directory already exists
    print(("Creation of the directory {} failed".format(outdir)))
    sys.exit()

#try:
#    os.chdir(outdir)
#except OSError:
#    print("Could not change directory to {}".format(outdir))
#    sys.exit()


def getBeam(ssidpath, ant):
    beam = glob.glob("{}/{}/beam*".format(ssidpath, ant))
    if len(beam)==0:
        print(("Error no beams in {}/{}".format(ssidpath, a)))
        sys.exit()
    if len(beam)>1:
        print(("Error multiple beams in {}/{}".format(ssidpath, a)))
        sys.exit()
    return(os.path.basename(beam[0]))

def makeSymlink(ssidpath,outdir,ant,beam):
    if os.path.isabs(ssidpath):
        os.symlink("{}/{}/{}".format(ssidpath,ant,beam), "{}/{}/{}".format(outdir,ant,beam))
    else:
        os.symlink("../../{}/{}/{}".format(ssidpath,ant,beam), "{}/{}/{}".format(outdir,ant,beam))

for a in ants1:
    try:
        os.mkdir("{}/{}".format(outdir,a))
    except OSError: 
        print(("Creation of the directory {} failed".format(a)))
        sys.exit()

    beam1 = getBeam(SSID1path, a)
    beam2 = getBeam(SSID2path, a)

    makeSymlink(SSID1path,outdir,a,beam1)
    makeSymlink(SSID2path,outdir,a,beam2)


