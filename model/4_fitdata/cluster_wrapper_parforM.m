function cluster_wrapper_parforM(modelname, binningfn, jobnum)
% a wrapper to fit multiple parameter values iteratively, so that the
% length of one job will be 48 hours.
%
% this function is used to combat the 500 job limit in the cluster. 
if isempty(binningfn); binningfn = 1; end % logistic is default

nStartVals = 1;

% load table with M numbers for particular job lists
Mjobs = dlmread('joblist_allM1.txt');


allsubjs = dlmread('joblist_04052016.txt');
subjVec = allsubjs(jobnum,:);
subjVec = subjVec(subjVec ~= 0);
for isubj = 1:length(subjVec);
    subjnum = subjVec(isubj);
    
    parfor iMVec = 1:size(Mjobs,1);
        MVec = Mjobs(iMVec,:);
        MVec = MVec(MVec ~= 0);
        
        for iM = 1:length(MVec);
            M = MVec(iM);
            fitdata_cluster(subjnum,modelname, binningfn, 'patternbayes', [1; M],[],[],nStartVals);
        end
    end
    
end
