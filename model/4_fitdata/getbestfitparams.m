function getbestfitparams(modelname, binningfn, subjids,filepath)
% gets the best fitting parameters from txt file across subjects and
% compiles it into a .mat file
% 
% MODELNAME: name of the model. FP, FPheurs, or uneqVar
% NSUBJ: number of subjects (e.g., 14, 25)
if nargin < 4; filepath = ['model' filesep '4_fitdata' filesep 'BPSfits' filesep]; end

switch modelname
    case {'FP','FPheurs','uneqVar'}
        nParams = 4;
    case 'REM'
        nParams = 7;
end

if any(binningfn == [2 3 4]); nParams = nParams + 2; end 
nLLcol = nParams + 1; % column corresponding to nLLs for FP, FPheurs, and uneqVar models
nSubj = length(subjids);

bestdata = nan(nSubj,nLLcol*2);
for isubj = 1:nSubj;
    subjid = subjids(isubj);
    filename = [filepath 'paramfit_patternbayes_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt'];
    alldata = dlmread(filename);
    
    datasorted = sortrows(alldata,nLLcol);
    bestdata(isubj,:) = datasorted(1,:);
    
end

bestFitParam = bestdata(:,1:nParams);
nLL_est = bestdata(:,nLLcol);
% startTheta = bestdata(:,6:9);
nLL_SD = bestdata(:,end);

matfilename = [ filepath 'paramfit_patternbayes_' modelname num2str(binningfn) '.mat'];
save(matfilename,'subjids','modelname','binningfn','bestFitParam','nLL_est','nLL_SD');