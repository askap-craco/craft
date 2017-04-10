#!/bin/bash

for f in `find . -name '*.hdr'` ; do
    d=`dirname $f`
    fname="${d}/fredda_galaxy3.cand"
    
    if [[ ! -e $fname ]] ; then
	echo $d
    fi
done
