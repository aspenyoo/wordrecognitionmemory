#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8GB
#PBS -j oe
#PBS -M ay963@nyu.edu
#PBS -m abe
#PBS -N paramfit2
#PBS -q s48

module purge
module load matlab/2014a

index=${PBS_ARRAYID}

export MATLABPATH=$HOME/matlab-scripts
cat<<EOF| matlab -nodisplay
addpath('/home/ay963/job-scripts')
addpath(genpath('/home/ay963/matlab-scripts'))
if $index > 14;
    blah = num2str($index);
    isubj = str2double(blah(1:end-2))
    ijob(2,1) = str2double(blah(end-1:end)) % which job number to do
end
joblist = dlmread('joblist_2102016.txt');
% fitdata_cluster(isubj, 'uneqVar','patternbayes');
fitdata_cluster(isubj,'FPheurs','patternbayes', joblist(ijob,:)); 
exit;
EOF