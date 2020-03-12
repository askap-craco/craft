#!/bin/bash

function remove_db {
    for key in $all_keys; do
	dada_db -d -k $key
    done
}

trap 'kill $(jobs -p) ; remove_db' EXIT
export MKL_DIR=/data/FRIGG_2/ban115/install/mkl/mkl/

export DADA=$HOME/psrdada-install
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$MKL_DIR/lib/intel64:$LD_LIBRARY_PATH

#export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdata
export DADA=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
export PATH=$DADA/bin:$PATH
export LD_LIBRARY_PATH=$DADA/lib:$LD_LIBRARY_PATH
#cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/src/cudafdmt
#cudafdmt=$HOME/craftdev/craft/cuda-fdmt/cudafdmt/Debugtest_cuda/cudafdmt
export SVD=/data/FRIGG_2/ban115/SpectralSVD
export PATH=$HOME/git/craft/cuda-fdmt/cudafdmt/Debugtest_cuda:$SVD:$PATH
#cudafdmt=$HOME/git/craft/cuda-fdmt/cudafdmt/src/cudafdmt
cudafdmt=`which cudafdmt`

ls -l $cudafdmt
# 84 samples/block
#block_size=8128512

# 512 samples/block
BLOCK_CYCLES=256
let block_size=72*336*4*${BLOCK_CYCLES}
let polsum_block_size=36*336*4*${BLOCK_CYCLES}

all_keys=""
DADA_KEY="1000"
ncount=1

for i in $(seq $ncount) ; do
    echo "COUNT $i"
    for indir in $@ ; do
	let DADA_KEY="$DADA_KEY + 100"
	all_keys="$all_keys $DADA_KEY"
	dada_db -d  -k $DADA_KEY > /dev/null 2>&1
	dada_db -a 32768 -b $block_size -n 4 -k $DADA_KEY  -l -p
	#$DADA/bin/dada_install_header -k $DADA_KEY -H $hdr
	#echo Header installed $hdr $DADA_KEY
	#echo dada_diskdb -k $DADA_KEY  -z -f $indir/*.dada
	dada_diskdb -k $DADA_KEY -z -f $indir/*00000.dada &

	if [[ $DADA_KEY == "1100" ]] ; then
	    echo dada_dbmonitor -k $DADA_KEY &
	fi
    done
done

rm -f *.dat
#cudafdmt -N 30 -t 512 -d 4096 -r 1  -s 0 -o fredda.multi.cand -p  -M 0.1 -T 0.1 -K 30 -B 16 -x 999999 $all_keys
#cudafdmt -D -N 30 -t 512 -d 512 -r 1  -s 0 -o /tmp/fredda.multi.cand -p  -B 16 -x 30 $all_keys

#$cudafdmt -t 512 -d 512 -r 1  -s 0  -M 0.2 -T 0.2 -C 6.0  -o fredda.$1.cand *.fil
#$cudafdmt -t 512 -d 512 -r 1  -s 0  -o fredda.$1.cand *.fil
CUDA_DEVICE=1
fredda_cand_file=fredda.cand
DADA_KEYS=$all_keys
icsdir=/tmp/

echo "DADA KEYS $DADA_KEYS"

which cudafdmt
ls -l `which cudafdmt`
#cuda-memcheck cudafdmt -N 3 -s 0.1 -S 0.0 -t $BLOCK_CYCLES -d 3072 -p -U localhost:$UDP_PORT  -g $CUDA_DEVICE -n 16384 -r 1 -T 2  -M 2 -K 2 -C 6 -x 10 -R -o $fredda_cand_file $fredda_extra $DADA_KEYS 

cudapid=$!
#cuda-gdb --args $cudafdmt -t 512 -d 512 $DADA_KEY -p -r 1 -D -r 1 -K 30 -s 0

input_keys=$all_keys
# Setup export buffer
ICS_EXPORT_KEY=1230
SVD_EXPORT_KEY=1240
FINAL_EXPORT_KEY=1250

all_keys="$input_keys $SVD_EXPORT_KEY $ICS_EXPORT_KEY $FINAL_EXPORT_KEY"
echo "BLOCK SIZE IS $block_size and polsum block size is $polsum_block_size"

echo "ALL KEYS $all_keys"
echo "INPUT KEYS" $input_keys
#dada_db -d -k $ICS_EXPORT_KEY
#dada_db -d -k $SVD_EXPORT_KEY
#dada_db -d -k $FINAL_EXPORT_KEY
dada_db -a 32768 -b $polsum_block_size -n 2 -k $ICS_EXPORT_KEY  -l -p -r 2
dada_db -a 32768 -b $polsum_block_size -n 2 -k $SVD_EXPORT_KEY  -l -p -r 2 
dada_db -a 32768 -b $polsum_block_size -n 2 -k $FINAL_EXPORT_KEY  -l -p
#rm -f 2*.dada

rm -f icsout/*.dada svdout/*.dada finalout/*.dada
mkdir icsout svdout finalout

dada_dbdisk -k $FINAL_EXPORT_KEY -z -D finalout &
dada_dbdisk -k $SVD_EXPORT_KEY -z -D svdout &
dada_dbdisk -k $ICS_EXPORT_KEY -z -D icsout &

mkdir fredda_ics
pushd fredda_ics
echo '************* STARTING FIRST FREDDA'
$cudafdmt -R  -t $BLOCK_CYCLES -d 2048 -p -r 1 -K 3 -I 2047 -W 3 -z 10 -C 10 -T 0.25 -M 0.007 -s 10 -o fredda_ics.cand -X $ICS_EXPORT_KEY $input_keys &
popd


dbsvddb  $ICS_EXPORT_KEY $SVD_EXPORT_KEY &
#dada_dbcopydb $ICS_EXPORT_KEY $SVD_EXPORT_KEY &

mkdir fredda_svd
pushd fredda_svd
flag_flags=-K 3 -M 0.1 -P -0.5 -Q 1.5 -A 0.7 -V 1.3 
echo $cudafdmt -t $BLOCK_CYCLES -d 2048 -r 1 -x 10 $flag_flags -X $FINAL_EXPORT_KEY -R -o fredda_svd.cand $SVD_EXPORT_KEY &
popd

#dada_dbmonitor -k $ICS_EXPORT_KEY &
dada_dbmonitor -k 1100 & 

#wait $cudapid
wait
#dada_db -d -k $DADA_KEY
kill $(jobs -p)
remove_db
