function getbestfitparams(modelname, nSubj)
% gets the best fitting parameters from txt file across subjects and
% compiles it into a .mat file
% 
% MODELNAME: name of the model. FP, FPheurs, or uneqVar
% NSUBJ: number of subjects (e.g., 14, 25)

nLLcol = 5; % column corresponding to nLLs for FP, FPheurs, and uneqVar models

bestdata = nan(nSubj,10);
for isubj = 1:nSubj;
    isubj
    filename = ['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt'];
    alldata = dlmread(filename);
    
    datasorted = sortrows(alldata,nLLcol);
    bestdata(isubj,:) = datasorted(1,:);
    
end

bestFitParam = bestdata(:,1:4);
nLL_est = bestdata(:,5);
% startTheta = bestdata(:,6:9);
nLL_SD = bestdata(:,10);

matfilename = ['4_fitdata/paramfit_patternbayes_' modelname '.mat'];
save(matfilename,'bestFitParam','nLL_est','nLL_SD');