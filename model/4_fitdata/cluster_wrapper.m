function cluster_wrapper(modelname, jobnum)
% a wrapper to fit multiple parameter values iteratively, so that the
% length of one job will be 48 hours.
%
% this function is used to combat the 500 job limit in the cluster. 

nStartVals = 1;

% load table with M numbers for particular job lists

nSubj = 14;
parfor isubj = 1:nSubj;
    
    alldata = dlmread('joblist_03302016.txt');
    MVec = alldata(jobnum,:);
    MVec = MVec(MVec ~= 0);
    
    for iM = 1:length(MVec);
        fixparam = MVec(iM)
        fitdata_cluster(isubj,modelname,'patternbayes', [1; fixparam],[],[],nStartVals);
    end
    
end