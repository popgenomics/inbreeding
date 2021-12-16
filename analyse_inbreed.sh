#!/bin/bash
binpath=/home/croux/Programmes/fasta2dos
python3 ${binpath}/fasta2dos.py ${1} > ${1}.dos 
python3 ${binpath}/hierf2inbreed.py ${1}.dos > ${1}.inbreedR
Rscript ${binpath}/inbred_stats.R input_hierfstat=${1}.dos input_inbreedR=${1}.inbreedR

