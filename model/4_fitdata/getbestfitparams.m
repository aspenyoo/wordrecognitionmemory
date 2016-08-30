function getbestfitparams(modelname, binningfn,memstrengthvar,subjids,paramrange, filepath)
% gets the best fitting parameters from txt file across subjects and
% compiles it into a .mat file
% 
% MODELNAME: name of the model. FP, FPheurs, or uneqVar
% NSUBJ: number of subjects (e.g., 14, 25)
% 
% PARAMRANGE: a 3 x (number of parameters you want to enforce a range)
% matrix, where the first row is the parameter number, second row is the
% lower bound of that parameter, and third row is the upper bound of that parameter. 
if nargin < 5; paramrange = []; end
if nargin < 6; filepath = ['model' filesep '4_fitdata' filesep 'BPSfits' filesep]; end

switch modelname
    case 'UVSD'
        nParams = 1;
    case {'FP','FPheurs'}
        nParams = 2;
    case 'REM'
        nParams = 5;
end
switch binningfn 
    case {0,1}
        nParams = nParams + 3; 
    case 2
        nParams = nParams + 4;
    case 3
        nParams = nParams + 5;
end

nLLcol = nParams + 1; % column corresponding to nLLs for FP, FPheurs, and uneqVar models
nSubj = length(subjids);

bestdata = nan(nSubj,nLLcol*2);
for isubj = 1:nSubj;
    subjid = subjids(isubj);
    filename = [filepath 'paramfit_patternbayes_' modelname num2str(binningfn) ...
        num2str(memstrengthvar) '_subj' num2str(subjid) '.txt'];
    alldata = dlmread(filename);
    
    datasorted = sortrows(alldata,nLLcol);
    for iparam = 1:size(paramrange,2); % how many ranges you are imposing
        idx = datasorted(:,iparam) < paramrange(2,iparam); % deleting things to small
        datasorted(idx,:) = [];
        
        idx = datasorted(:,iparam) > paramrange(3,iparam); % deleting things to large
        datasorted(idx,:) = [];
    end
    bestdata(isubj,:) = datasorted(1,:);
    
end

bestFitParam = bestdata(:,1:nParams);
nLL_est = bestdata(:,nLLcol);
% startTheta = bestdata(:,6:9);
nLL_SD = bestdata(:,end);

matfilename = [ filepath 'paramfit_patternbayes_' modelname num2str(binningfn) num2str(memstrengthvar) '.mat'];
save(matfilename,'subjids','modelname','binningfn','bestFitParam','nLL_est','nLL_SD');