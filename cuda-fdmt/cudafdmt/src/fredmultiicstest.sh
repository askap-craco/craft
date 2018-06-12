#!/bin/bash

trap 'kill $(jobs -p)' EXIT
infile=`ls *.dada`
hdr=`ls co*.hdr`
echo Infile $infile hdr=$hdr
#cat $hdr
export DADA=$HOME/psrdada-install
#export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdata
#export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
#cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt

#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/Debugtest_cuda/cudafdmt
#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/src/cudafdmt

ls -l $cudafdmt
# 84 samples/block
#block_size=8128512

# 512 samples/block
let block_size=72*336*4*512

all_keys=""
DADA_KEY="1000"

for indir in $@ ; do
    let DADA_KEY="$DADA_KEY + 100"
    all_keys="$all_keys $DADA_KEY"
    dada_db -d  -k $DADA_KEY > /dev/null 2>&1
    dada_db -a 32768 -b $block_size -n 8 -k $DADA_KEY 
    #$DADA/bin/dada_install_header -k $DADA_KEY -H $hdr
    #echo Header installed $hdr $DADA_KEY
    echo dada_diskdb -k $DADA_KEY  -z -f $indir/*.dada
    dada_diskdb -k $DADA_KEY -z -f $indir/*.dada &
    dada_dbmonitor -k $DADA_KEY &
done

rm -f *.dat

cudafdmt -N 70 -t 512 -d 512 -r 1  -s 0 -o fredda.multi.cand -p  -M 0.1 -T 0.1 -K 30 $all_keys &

#$cudafdmt -t 512 -d 512 -r 1  -s 0  -M 0.2 -T 0.2 -C 6.0  -o fredda.$1.cand *.fil
cudapid=$!
#cuda-gdb --args $cudafdmt -t 512 -d 512 $DADA_KEY -p -r 1 -D -r 1 -K 30 -s 0

#wait $cudapid
wait
#dada_db -d -k $DADA_KEY
kill $(jobs -p)
for key in all_keys; do
    dada_db -d $key
done

