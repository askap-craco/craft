#!/bin/bash

#    frb_idm=$1
#    t0=$2 # start time in samples
#    frb_amp=$3
#    frb_sn=$4
#    frb_relpos=$5

craco_sim_single.sh 0 0 1 inf 0,0
# Canonical:  unit FRB at center at 0 DM
craco_sim_single.sh 0 1 1 inf 0,0 # canonical + toff=1 samlpes
craco_sim_single.sh 0 2 1 inf 0,0 # canonical + toff=2 samlpes
craco_sim_single.sh 0 3 1 inf 0,0 # canonical + toff=3 samlpes
craco_sim_single.sh 0 4 1 inf 0,0 # canonical + toff=4 samlpes
craco_sim_single.sh 0 9 1 inf 0,0 # canonical + toff=9 samlpes
craco_sim_single.sh 2 0 1 inf 0,0 # canonical + dm=2
craco_sim_single.sh 0 0 2 inf 0,0 # canonical + amplitude=2
craco_sim_single.sh 0 0 1 inf 100,200 # canonical + position offset = 100" in ra and 200" in dec
craco_sim_single.sh 0 0 1 10 0,0 # Canonical + S/N = 12

craco_sim_single.sh 6 2 1 inf 0,0 #
craco_sim_single.sh 6 2 1 inf 300,400 # 
craco_sim_single.sh 6 2 1 10  300,400 # 
