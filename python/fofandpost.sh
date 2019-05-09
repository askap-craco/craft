#!/bin/bash

fname=$1
fredfof_ics.py $@
plot_fredda_cand.py $fname -o ${fname}.png
postslack.py -m 'Scan summary' scan=SB0${SBID}/${SCANNNAME}/ICS/${CID} -i ${fname}.png

