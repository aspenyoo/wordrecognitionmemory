#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=4GB
#PBS -m abe
#PBS -N REM30

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
binningfn = 3;
memstrengthvar = 0;

if $index > 100;
    blah = num2str($index);
    isubj = str2double(blah(1:end-2));
    joblistnum = str2double(blah(end-1:end)); 
end
joblistfile = 'joblist_08112016.txt'; % will do a different job per person. (if you want same for all, write joblist name)

cluster_wrap2(modelname, binningfn, memstrengthvar, isubj, joblistnum,joblistfile,[1 7; nan 0])

EOF


