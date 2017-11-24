#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

export CRAFT=/home/ban115/craft/craft/
export PATH=$CRAFT/cuda-fdmt/cudafdmt/src:$CRAFT/python:$PATH:/home/ban115/bin:$CRAFT/jobs/
export PYTHONPATH=$CRAFT/python:$PYTHONPATH
export OMP_NUM_THREADS=24

module unload PrgEnv-cray > /dev/null
module load python/2.7.10 > /dev/null
module load astropy > /dev/null

thedir=$1
thebase=`basename $thedir`

echo `date` starting fixheaders

CRAFTDATA=/scratch2/askap/askapops/craft/co/

ls -d $CRAFTDATA/SB00* | parallel -j 24 fixhdr.sh

echo `date` finished fixheaders
