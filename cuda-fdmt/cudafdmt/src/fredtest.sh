#!/bin/bash

trap 'kill $(jobs -p)' EXIT

DADA_KEY=babb
infile=/DATA/SERPENS_2/ban115/craft/rtfredda/co25.dada
hdr=/DATA/SERPENS_2/ban115/craft/rtfredda/co25*.hdr
echo Infile $infile hdr=$hdr
cat $hdr
export DADA=$HOME/psrdada-install
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
$DADA/bin/dada_db -d  -k $DADA_KEY > /dev/null 2>&1
$DADA/bin/dada_db -a 32768 -b 8128512 -n 16 -k $DADA_KEY
#$DADA/bin/dada_install_header -k $DADA_KEY -H $hdr
#echo Header installed $hdr $DADA_KEY
#../Debug\ test_cuda/cudafdmt -t 512 -d 1024 $DADA_KEY &
echo dada_diskdb -k $DADA_KEY -f $infile -z
dada_diskdb -k $DADA_KEY -f $infile -z &
dada_dbmonitor -k $DADA_KEY
wait

