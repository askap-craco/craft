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
#export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
export PATH=$DADA/bin:$PATH
#export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
#cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt

#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/Debugtest_cuda/cudafdmt
cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt

ls -l $cudafdmt
# 84 samples/block
#block_size=8128512

# 512 samples/block
nt=512
let block_size="72*336*4*$nt"

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

let out_block_size="36*336*4*$nt"
#OUT_KEY=8001
#$DADA/bin/dada_db -d -k $OUT_KEY > /dev/null 2>&1
#$DADA/bin/dada_db -a 32768 -b $out_block_size -n 4 -k $OUT_KEY
#mkdir -p  icsout/
#rm -f icsout/*.dada
#$DADA/bin/dada_dbdisk -z -k $OUT_KEY -D icsout/ &
#all_keys="$all_keys $OUT_KEY"

#$cudafdmt -N 5 -t $nt -d 512 -r 1  -s 0 -o fredda.$1.cand -D -p $all_keys
#cmd="valgrind --leak-check=full -v  cudafdmt -t 128 -d 128 -p -N 2 $all_keys"
#Runnding $cmd
#$cmd
wait

#dada_db -d -k $DADA_KEY


