

%% SIMULATE DATA / CALCULATE LL OF DATA

clear all

modelname = 'FP';
binningfn = 3;

nSubjs = 14;

% simulate data 
for isubj = 1:nSubjs
    
    if (isubj)
        load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d_freed0.mat',modelname,binningfn))
        theta = bestFitParam(isubj,:);
    else
        switch modelname
            case 'UVSD'
            case 'FP'
                theta = [5 .1]; % M sigma
            case 'REM'
                
        end
        
        switch binningfn
            case 3 % power law
                theta = [theta 1 0 .1 0 .1]; % a b gamma d0 sigma_mc
            case 4 % weibull
                theta = [theta 1 0 .1 .1 0 .1]; % a b shape scale d0 sigma_mc
        end
    end
    
    % number of old and new ratings for participant (only relevant aspect when
    % simulating data is that sum(nnew_part) = sum(nold_part) = 150)
    [nnew_part, nold_part] = deal([150 zeros(1,19)]);
    
    
    [ pnewMat(isubj,:), pold(isubj,:) ] = calc_nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part,[],10,10);%, logflag, nConf )

end

new_modelM = mean(pnewMat);
new_modelSD = std(pnewMat)./sqrt(nSubjs);


% % MAKE FIGURE
% oldcolor = [.7 .7 .5];
% newcolor = [.5 .4 .7];

% figure;
% plot_summaryfit(1:20,partM,partSD,new_modelM,new_modelSD,newcolor,newcolor);
% hold on;
% plot_summaryfit(1:20,partM,partSD,old_modelM,old_modelSD,oldcolor,oldcolor);
% defaultplot

%% FIT DATA

clear all

nSubjs = 14;
modelname = 'UVSDd';
binningfn = 3;

nStartVals = 10;

switch modelname
    case 'UVSDd' % no sigma_mc
        nStartVals = 2;
        switch binningfn
            case 3
                fixparams = [7; 0];
            case 4
                fixparams = [8; 0];
        end
    case 'UVSDx'
        nStartVals = 2;
end

for isubj = 11:nSubjs
    fitdata_cluster(isubj, modelname, binningfn, fixparams,[], 20, nStartVals)
end

% [bestFitParam, nLL_est, startTheta, outputt] = paramfit_patternbayes(testmodelname, binningfn, nnew_part, nold_part, fixparam ,1);

%% get best fit params

clear all
modelname = 'UVSDd';
binningfn = 4;
subjids = 1:14;

getbestfitparams(modelname,binningfn,subjids);

%%  PLOT MODEL VALIDATION PLOTS
% group confidence ratings, indvl confidence ratings, indvl d distributions

clear all
modelname = 'UVSDd';
binningfn = 4;

if strcmp(modelname,'UVSDx') || strcmp(modelname,'UVSDd')
    additionalstr = '';
else
    additionalstr = '_freed0';
end

load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d%s.mat',modelname,binningfn,additionalstr));


plotparamfits(modelname, bestFitParam, binningfn);%, nConf, fakedata, subjnum, selectiveplot)












