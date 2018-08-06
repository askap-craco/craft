#!/bin/bash

for f in $@ ; do
    cat $f > /dev/null &
done
wait
