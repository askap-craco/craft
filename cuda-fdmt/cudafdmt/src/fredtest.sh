#!/bin/bash

trap 'kill $(jobs -p)' EXIT

DADA_KEY=babb
infile=`ls *.dada`
hdr=`ls co*.hdr`
echo Infile $infile hdr=$hdr
#cat $hdr
#export DADA=$HOME/psrdada-install
export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdata
export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
#cudafdmt=$HOME/craft/craft/cuda-fdmt/cudafdmt/Debugtest_cuda/cudafdmt
cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt
#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/src/cudafdmt


ls -l $cudafdmt
# 84 samples/block
block_size=8128512

# 512 samples/block
let block_size=72*336*4*512
$DADA/bin/dada_db -d  -k $DADA_KEY > /dev/null 2>&1
$DADA/bin/dada_db -a 32768 -b $block_size -n 16 -k $DADA_KEY 
#$DADA/bin/dada_install_header -k $DADA_KEY -H $hdr
#echo Header installed $hdr $DADA_KEY
$cudafdmt -t 512 -d 512 $DADA_KEY -p -r 1 -D -r 1 -K 30 -s 0  &
echo dada_diskdb -k $DADA_KEY -f $infile -z
dada_diskdb -k $DADA_KEY -f $infile -z &
dada_dbmonitor -k $DADA_KEY
wait

