#!/bin/bash

sbid=$1
scan=$2
ant=$3

localdir=/Users/ban115/bolton/craft/data/
remotedir=galaxy:/scratch2/askap/askapops/craft/co/

d=$sbid/$scan/$ant
mkdir -p $localdir/$d
cd $localdir/$d
pwd
cmd="rsync -avz --progress $remotedir/$d/ ."
echo $cmd
$cmd
pwd $localdir/$d
