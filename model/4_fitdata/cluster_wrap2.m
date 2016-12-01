function cluster_wrap2(modelname, binningfn, subjids, jobnum, joblistfile, fixparams)
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

if nargin < 5; fixparams = [1; nan]; end

nSubj = length(subjids);
nStartVals = 1;

% load table with M numbers for particular job lists

for isubj = 1:nSubj;
    
    subjid = subjids(isubj);
    
    if ~isempty(joblistfile)
        alldata = dlmread(joblistfile);
    else
        alldata = dlmread(['joblist_' modelname '_subj' num2str(subjid) '.txt']);
    end
    MVec = alldata(jobnum,:);
    MVec = MVec(MVec ~= 0);
    
    for iM = 1:length(MVec);
        fixparams(2,1) = MVec(iM)
        fitdata_cluster(subjid, modelname, binningfn, fixparams,[],[],nStartVals);
    end
    
end