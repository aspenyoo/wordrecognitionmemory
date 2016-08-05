#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -l walltime=36:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=28GB
#PBS -m abe
#PBS -N paramfit_wrapper_REM

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
binningfn = 1;

if $index > 100;
    blah = num2str($index);
    isubj = str2double(blah(1:end-2));
    joblistnum = str2double(blah(end-1:end)); 
else
    isubj = 1:14;
    joblistnum = $index;
end

cluster_wrapper(modelname, binningfn, isubj, joblistnum)

EOF


