#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=20GB
#PBS -m abe
#PBS -N paramfit

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

modelname = 'REM';
nStartVals = 10;

fixparam = $index;

parfor isubj = 8:14;
	fitdata_cluster(isubj,modelname,'patternbayes', [1; fixparam],[],nStartVals); exit;
end

EOF


