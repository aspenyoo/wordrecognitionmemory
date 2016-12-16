function getbestfitparams(modelname,binningfn,subjids,nStartVals,paramrange,filepath)
% gets the best fitting parameters from txt file across subjects and
% compiles it into a .mat file
% 
% MODELNAME: name of the model. FP, FPheurs, or uneqVar
% NSUBJ: number of subjects (e.g., 14, 25)
% 
% PARAMRANGE: a 3 x (number of parameters you want to enforce a range)
% matrix, where the first row is the parameter number, second row is the
% lower bound of that parameter, and third row is the upper bound of that parameter. 
if nargin < 4; nStartVals = []; end
if nargin < 5; paramrange = []; end
if nargin < 6; filepath = ['model' filesep '4_fitdata' filesep 'BPSfits' filesep]; end

Mmax = 50;

switch modelname
    case 'UVSD'
        nParams = 2;
    case {'FP','FPheurs'}
        nParams = 2;
    case 'REM'
        nParams = 5;
end

switch binningfn
    case {0,1} % linear, logistic
        % slope, y-int, sigma_mc
        nParams = nParams + 3;
    case 2     % logarithmic
        % a, b, d0, sigma_mc
        nParams = nParams + 4;
    case 3      % power law
        % a, b, gamma, d0, sigma_mc
        nParams = nParams + 5;
    case 4 % weibull
        % a, b, shape, scale, d0, sigma_mc
        nParams = nParams + 6;
end

nLLcol = nParams + 1; % column corresponding to nLLs for FP, FPheurs, and uneqVar models
nSubj = length(subjids);

bestdata = nan(nSubj,nLLcol*2);
for isubj = 1:nSubj
    subjid = subjids(isubj);
    filename = [filepath 'paramfit_patternbayes_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt'];
    alldata = dlmread(filename);
    
    % deleting entries that exceed nStartVals
    if ~isempty(nStartVals)
        for iM = 1:Mmax
            idx = find(alldata(:,1) == iM);
            if length(idx) > nStartVals
                alldata(idx(nStartVals+1:end),:) = [];
            end
        end
    end
    
    datasorted = sortrows(alldata,nLLcol);
    for iparam = 1:size(paramrange,2) % how many ranges you are imposing
        idx = datasorted(:,iparam) < paramrange(2,iparam); % deleting things too small
        datasorted(idx,:) = [];
        
        idx = datasorted(:,iparam) > paramrange(3,iparam); % deleting things too large
        datasorted(idx,:) = [];
    end
    
    bestdata(isubj,:) = datasorted(1,:);
    
end

bestFitParam = bestdata(:,1:nParams);
nLL_est = bestdata(:,nLLcol);
% startTheta = bestdata(:,6:9);
nLL_SD = bestdata(:,end);

matfilename = [ filepath 'paramfit_patternbayes_' modelname num2str(binningfn) '.mat'];
save(matfilename,'subjids','modelname','bestFitParam','nLL_est','nLL_SD');