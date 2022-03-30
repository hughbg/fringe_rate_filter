#!/bin/bash
#PBS -q hera
#PBS -j oe
#PBS -N imaging
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -l mem=6GB
#PBS -M h.garsden@qmul.ac.uk

infile="cleaned_sim_small.uvh5"
outfile="xxxxxxsim_middle.png"

echo "$infile -> $outfile"

cd /lustre/aoc/projects/hera/hgarsden/fringe_rate_filter

~/anaconda3_all/bin/python extract_times.py $infile
rm -rf x.* output* 
cp middle.uvfits x.uvfits
casa --nogui --nologger < mk_image.py
cp output.png $outfile
