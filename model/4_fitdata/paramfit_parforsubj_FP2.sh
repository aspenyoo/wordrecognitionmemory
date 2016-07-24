#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=30GB
#PBS -m abe
#PBS -N paramfit_parforsubj_FP2

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
fixparams = [1; $index];
optimMethod = 'patternbayes';

clusterwrapper_parforsubj2(modelname,binningfn,nStartVals,fixparams,optimMethod)

EOF

