#!/bin/sh -l 
# FILENAME:  myqsub2.sh
#PBS -q standby
#PBS -l nodes=1:ppn=1 
#PBS -l walltime=4:00:00

# Change to the directory 
cd $PBS_O_WORKDIR
cd .. 
./FiberMpp 2
