function cluster_wrapper(modelname, binningfn, subjids, jobnum)
% a wrapper to fit multiple parameter values iteratively, so that the
% length of one job will be 48 hours.
%
% this function is used to combat the 500 job limit in the cluster. 
% 
% ======== INPUT VARIABLES ========
% MODELNAME: 'FP','FPheurs','UVSD','REM'
% BINNINGFN: 0 (linear binning), 1 (logistic binning), 2 (log binning)
% SUBJIDS: vector of subject IDs
% JOBNUM: which row of the joblist you want to be evaluated

if isempty(binningfn); binningfn = 1; end

nSubj = length(subjids);
nStartVals = 1;

% load table with M numbers for particular job lists

samejoblist = 1;
parfor isubj = 1:nSubj;
    
    subjid = subjids(isubj);
    
    if (samejoblist)
        alldata = dlmread('joblist_07252016.txt');
    else
        alldata = dlmread(['joblist_patternbayes_subj' num2str(subjid) '_maxtime48.txt']);
    end
    MVec = alldata(jobnum,:);
    MVec = MVec(MVec ~= 0);
    
    for iM = 1:length(MVec);
        fixparam = MVec(iM)
        fitdata_cluster(subjid, modelname, binningfn, 'patternbayes', [1; fixparam],[],[],nStartVals);
    end
    
end
