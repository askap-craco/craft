#!/bin/bash

trap 'kill $(jobs -p)' EXIT

fno=0
for f in `find . -name '*.fil'` ; do 
    
    echo "`date` Processing $f"

    dout=$f.heim

    if [[ -d $dout ]] ; then
	echo "Output directory $dout exists. Ignorring"
	continue
    fi
    mkdir $dout
    gpuno=$(($fno % 2))
    cmd="heimdall -output_dir $dout -zap_chans 184 184 -dm 0 4000 -f $f -nsamps_gulp  16384 -v"
	
    echo $cmd
	# 14 MB file of around 30s takes 7 s with a single heimdall
	# Using default gulp = 262144
	# All 36 beams fails if you  run it in the background
	# at blocksize=2 it takes 6s per file
	# at blocksz=6 it takes 10s per hfile
	# at blocksz=12 it takes 18s per file
	# at blocksz=18 it breaks with a memory allocation failure
    
	# With nsamps_gulp=16385
	# with blocksz=12 it takes 18s per file
	# with blocksz=18 it takes 21s per file
	# with blocksz=24it takes 21s per file
	# 36 = mem allocation failure

	
    blocksz=$((1))
    endblock=$(($blocksz - 1))
    blockno=$(($fno % $blocksz))
    
    $cmd &

    if [[ $blockno ==  $endblock ]] ; then
	echo "Waiting for heidall to finish"
	wait
	echo "`date` Finished $blocksz beams"
    fi
    fno=$(($fno + 1))
    echo "`date` Finished processing $dname"

done

wait
