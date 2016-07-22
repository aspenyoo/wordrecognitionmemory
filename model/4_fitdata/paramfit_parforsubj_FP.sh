#!/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=10GB
#PBS -m abe
#PBS -N paramfit_parforsubj

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
nStartVals = 10;

fixparam = $index;

parfor isubj = 1:14;
	nStartVal = max([nStartVals-countnum(modelname,isubj,fixparam) 0]);
	try
	fitdata_cluster(isubj,modelname,binningfn,'patternbayes', [1; fixparam],[],[],nStartVal); exit;
	end
end

EOF


