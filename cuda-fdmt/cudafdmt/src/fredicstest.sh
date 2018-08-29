#!/bin/bash


function remove_db {
    for key in $all_keys; do
	dada_db -d -k $key
    done
}

trap 'kill $(jobs -p) ; remove_db' EXIT
infile=`ls *.dada`
hdr=`ls co*.hdr`
echo Infile $infile hdr=$hdr
#cat $hdr
export DADA=$HOME/psrdada-install
#export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdata
export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
#cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt

#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/Debugtest_cuda/cudafdmt
cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt

ls -l $cudafdmt
# 84 samples/block
#block_size=8128512

# 512 samples/block
let block_size=72*336*4*512

all_keys=""

for keyoff in $@ ; do
    let DADA_KEY="1000 + 100*${keyoff}"
    all_keys="$all_keys $DADA_KEY"

    $DADA/bin/dada_db -d  -k $DADA_KEY > /dev/null 2>&1
    $DADA/bin/dada_db -a 32768 -b $block_size -n 8 -k $DADA_KEY 
    #$DADA/bin/dada_install_header -k $DADA_KEY -H $hdr
    #echo Header installed $hdr $DADA_KEY
    echo dada_diskdb -k $DADA_KEY -f $infile -z
    dada_diskdb -k $DADA_KEY -f $infile -z &
    dada_dbmonitor -k $DADA_KEY &

done

rm -f *.dat
$cudafdmt -t 512 -d 512 -r 4 -o fredda.$1.cand -p -x 10 -M 0.1 -T 0.1 -C 20  $all_keys
#$cudafdmt -t 512 -d 512 -r 1  -s 0  -M 0.2 -T 0.2 -C 6.0  -o fredda.$1.cand *.fil
cudapid=$!
#cuda-gdb --args $cudafdmt -t 512 -d 512 $DADA_KEY -p -r 1 -D -r 1 -K 30 -s 0

#wait $cudapid
wait
#dada_db -d -k $DADA_KEY


