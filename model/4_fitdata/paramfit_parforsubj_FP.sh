#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=30GB
#PBS -m abe
#PBS -N paramfit_parforsubj_FP

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
nStartVals = 10;

parfor isubj = 1:14;
fixparam = $index;
	nStartVal = max([10-countnum('FP',isubj,fixparam) 0]);
	fitdata_cluster(isubj,'FP',2,'patternbayes', [1 6; fixparam 0],[],[],nStartVal); exit;
end

EOF


