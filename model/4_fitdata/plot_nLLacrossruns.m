function plot_nLLacrossruns(modelname, binningfn, subjids,filepath)
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

figure;
for isubj = 1:nSubj;
    subjid = subjids(isubj);
    filename = [filepath 'paramfit_patternbayes_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt'];
    alldata = dlmread(filename);
    
    subplot(4,4,subjid)
    scatter(alldata(:,1),alldata(:,nLLcol),'k.');
    defaultplot
end