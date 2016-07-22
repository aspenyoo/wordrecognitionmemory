#!/bin/bash
#PBS -l nodes=1:ppn=13
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=26GB
#PBS -m abe
#PBS -N paramfit_wrapper

index=${PBS_ARRAYID}
job=${PBS_JOBID}
walltime_lim=${PBS_WALLTIME}
script_name=${PBS_JOBNAME}
module purge
module load matlab

export MATLABPATH=/home/ay963/matlab-scripts
cat<<EOF | matlab -nodisplay
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/wordrecognitionmemory'))

modelname = 'FP';
binningfn = 2;
subjids = 1:14;

cluster_wrapper(modelname, binningfn, subjids, $index+1)

EOF


