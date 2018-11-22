#!/bin/bash

hdrdir=/data/FRIGG_2/ban115/voltages/0407/20181112222901
for c in {1..7} ; do
    solvefordelays.py -p $hdrdir/ak*/beam00/*.hdr -d 20181112222901_c${c}_beam00_i1024_f1.uvspec.real -b 0  -i 7034

    # delays are all in card 7 - grep and sed it to the desired card
    grep _c7_ hwdelays_SB7034*.real.txt | sed s/_c7_/c${c}/ >> hwdelays.txt

done
