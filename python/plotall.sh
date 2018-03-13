#!/bin/bash

fname=$1

for f in $@ ; do
    pushd $f
    
    
    popd
done
