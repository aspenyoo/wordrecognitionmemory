function cluster_wrapper(modelname, subjids, jobnum)
% a wrapper to fit multiple parameter values iteratively, so that the
% length of one job will be 48 hours.
%
% this function is used to combat the 500 job limit in the cluster. 

nStartVals = 1;

% load table with M numbers for particular job lists

nSubj = length(subjids);
samejoblist = 0;
parfor isubj = 1:nSubj;
    
    subjid = subjids(isubj);
    
    if (samejoblist)
        alldata = dlmread('joblist_03302016.txt');
    else
        alldata = dlmread(['joblist_patternbayes_subj' num2str(subjid) '_maxtime48.txt']);
    end
    MVec = alldata(jobnum,:);
    MVec = MVec(MVec ~= 0);
    
    for iM = 1:length(MVec);
        fixparam = MVec(iM)
        fitdata_cluster(subjid,modelname,'patternbayes', [1; fixparam],[],[],nStartVals);
    end
    
end