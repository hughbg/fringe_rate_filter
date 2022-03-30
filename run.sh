#!/bin/sh
#SBATCH --job-name='frf'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=232GB
#SBATCH --output=frf_sh_%j.log
#SBATCH --time=72:00:00
#SBATCH -p hera


if [ -z "$SLURM_JOBID" ]
then SLURM_JOBID=shell
fi

cd /lustre/aoc/projects/hera/hgarsden/fringe_rate_filter

date 

files=`cat files.txt`
beam="/lustre/aoc/projects/hera/H4C/beams/NF_HERA_Vivaldi_efield_beam_healpix.fits"
spw="180 265"

input=sim_1000_256_noise.uvh5
output=cleaned_sim_1000_256_noise.uvh5

echo "$input -> $output" > frf_${SLURM_JOBID}.log

~/anaconda3_all/bin/python hera_cal/scripts/tophat_frfilter_run.py $input --tol 1e-9 \
      --CLEAN_outfilename $output \
      --uvbeam $beam --percentile_low 5 --percentile_high 95\
      --clobber --verbose --mode dft_leastsq --frate_standoff 0.05 --min_frate_half_width 0.15 >>  frf_${SLURM_JOBID}.log 2>&1

date
