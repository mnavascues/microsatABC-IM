#!/bin/bash
# définir le type de shell ?????
#$ -S /bin/bash

# Set name of job:
#$ -N Agrarius

# Chose queue: 
#$ -q mp.q
#$ -pe SMP 8

# change working directory (to current directory, 
# otherwise it must be specified): 
#$ -cwd

set -x 

# get current working directory to variable workingDir
workingDir=$PWD

# copy files to temporal directory on the working machine
cp -rp *.R $TMPDIR/
cp -rp *.str $TMPDIR/
cp -rp ms $TMPDIR/

# change direcory
cd $TMPDIR 

chmod +x ms
chmod +x main.R

# run script in R
./main.R > Agrarius.log

# copy files (and results) back to original folder
cp -rp $TMPDIR/* $workingDir/











