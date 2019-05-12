#!/bin/bash

echo Doing $0 $@ scan=SB0${SBID}/${SCANNAME}/ICS/${CID} alias=${SB_ALIAS}
fname=$1
fredfof_ics.py $@
plot_fredda_cand.py $fname.fof -o ${fname}.fof.png
postslack.py -t 'Scan summary' scan=SB0${SBID}/${SCANNAME}/ICS/${CID} alias=${SB_ALIAS} -i ${fname}.png

