#!/bin/bash

for d in 512 1024 2048 4096 ; do
    t=512
    while [[ $t -le $d ]] ; do
	nbeams=`ls *.fil | wc -l`
	for b in {1..11} ; do
	    echo $d $t $b 
	    files=`ls *.fil | head -n $b`
	    logf=cudalog_t${t}_d${d}_b${b}.log
	    echo logf $logf
	    cudafdmt -t $t -d $d $files > $logf 
	done
    done
done
