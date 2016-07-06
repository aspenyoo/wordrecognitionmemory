#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -M aspen.yoo@nyu.edu
#PBS -l mem=2GB
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
addpath('/home/ay963/job-scripts')
addpath(genpath('/home/ay963/matlab-scripts'))
addpath(genpath('/home/ay963/wordrecognitionmemory'))

modelname = 'FPheurs';
nStartVals = 10;

% used in cluster to indicate both subject and fixed M
if $index > 100;
    blah = num2str($index);
    isubj = str2double(blah(1:end-2));
    fixparam = str2double(blah(end-1:end)); % fixing M value
    nStartVals = max([nStartVals-countnum(modelname,isubj,fixparam) 0]); % nStartVals. used for resubmitting
end

% fitdata_cluster(isubj, 'uneqVar','patternbayes');
fitdata_cluster(isubj,modelname,'patternbayes', [1; fixparam],[],nStartVals); exit;

status = 0;
while true
	if aborted % continue same model, subject, chain
		try mat = dlmread
		status=system(sprintf('cd /home/ay963/job-scripts; qsub -t $index %s',script_name))
	end
	if status==0; break; else; pause(65); end
end
EOF


