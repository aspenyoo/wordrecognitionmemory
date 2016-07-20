#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:01:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=1GB
#PBS -m abe
#PBS -N testing

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

fixparam = $index;
testmodelname = 'REM';
isubj = 1;

filepath = 'BPSfits/';
filename = [filepath 'paramfit_patternbayes_' testmodelname '_subj' num2str(isubj) '.txt'];

fileID = fopen(filename,permission);
formatSpec = '%4.4f testing \r\n';
fprintf(fileID, formatSpec, index); % save stuff in txt file
fclose(fileID);

EOF


