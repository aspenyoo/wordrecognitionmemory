#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=38:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=4GB
#PBS -m abe
#PBS -N paramfit_wrap2

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

if $index > 100;
    blah = num2str($index);
    model = str2double(blah(1));
    binningfn = str2doubld(blah(2));
    isubj = str2double(blah(3:4));
    joblistnum = str2double(blah(5:6)); 
end

switch model
case 1
	modelname = 'FP';
	d0idx = 4;
case 2
	modelname = 'FPheurs';
case 3
	modelname = 'REM';
	d0idx = 7;
end

if any(binningfn == [2 3 4]); d0idx = d0idx + 1; end

cluster_wrap2(modelname, binningfn, isubj, joblistnum,'joblist_07252016.txt',[1 d0idx; nan 0])

EOF


