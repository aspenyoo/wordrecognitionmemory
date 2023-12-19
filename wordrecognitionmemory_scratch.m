

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

nStartVals = 3;

switch modelname
    case 'UVSDd' % no sigma_mc
        switch binningfn
            case 3
                fixparams = [7; 0];
            case 4
                fixparams = [8; 0];
        end
    case 'UVSDx'
        fixparams = [];
end

for isubj = 1:nSubjs
    fitdata_cluster(isubj, modelname, binningfn, fixparams,[], 20, nStartVals)
end

% [bestFitParam, nLL_est, startTheta, outputt] = paramfit_patternbayes(testmodelname, binningfn, nnew_part, nold_part, fixparam ,1);

%% get best fit params

clear all
modelname = 'UVSDx';
binningfn = 3;
subjids = 1:14;

getbestfitparams(modelname,binningfn,subjids);

%%  PLOT MODEL VALIDATION PLOTS
% group confidence ratings, indvl confidence ratings, indvl d distributions

clear all
modelname = 'FP';
binningfn = 3;

% if strcmp(modelname,'UVSDx') || strcmp(modelname,'UVSDd')
    additionalstr = '';
% else
%     additionalstr = '_freed0';
% end

load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d%s.mat',modelname,binningfn,additionalstr));

plotparamfits(modelname, bestFitParam, binningfn);%, nConf, fakedata, subjnum, selectiveplot)

%% PARAMETER VALUES

clear all

modelname = 'UVSDx';
binningfn = 3;

load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d.mat',modelname,binningfn))
% load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d_freed0.mat',modelname,binningfn))

switch modelname
    case 'FP'
        logflag = [1 1];
    case 'REM'
        logflag = [1 0 0 0 0];
    case {'UVSDx','UVSDd'}
        logflag = [0 1];
end
switch binningfn
    case 3 % power law
        logflag = [logflag 0 0 0];
    case 4 % cumulative weibull
        logflag = [logflag 0 0 0 0];
end
logflag = [logflag 0 1];
logflag = logical(logflag);

bestFitParam(:,logflag) = exp(bestFitParam(:,logflag));
median(bestFitParam)
prctile(bestFitParam,[25 75])


%% MODEL COMPARISON

clear all

modelnameVec = {'UVSDx','FP','REM'};
binningfnVec = [3 4];
nSubjs = 14;
n = 300;
imodelref = 5;

nModels = length(modelnameVec);
nbinningfns = length(binningfnVec);

for ibinningfn = 1:nbinningfns
    binningfn = binningfnVec(ibinningfn);
    
    nLLMat{ibinningfn} = nan(nModels,nSubjs);
    for imodel = 1:nModels
        modelname = modelnameVec{imodel};
        
        load(sprintf('fitdata/BPSfits/paramfit_patternbayes_%s%d.mat',modelname,binningfn))
        
        nParamsVec{ibinningfn}(imodel) = size(bestFitParam,2);
        nLLMat{ibinningfn}(imodel,:) = nLL_est;
    end
end

LLMat = [nLLMat{1}; nLLMat{2}];
nParamsVec = [nParamsVec{1} nParamsVec{2}];

AIC = 2*bsxfun(@plus,LLMat,nParamsVec');
BIC = bsxfun(@plus,2*LLMat,log(n)*nParamsVec');
AICc = bsxfun(@plus,AIC,2.*bsxfun(@rdivide,nParamsVec'.*(nParamsVec'+1),(bsxfun(@minus,n,nParamsVec')-1)));

Delta_AICc = bsxfun(@minus,AICc,AICc(imodelref,:));
Delta_BIC = bsxfun(@minus,BIC,BIC(imodelref,:));

med_AICc = median(Delta_AICc,2);
med_BIC = median(Delta_BIC,2);

% if all models same number of parameters, do LL comp
if (sum(nParamsVec == min(nParamsVec)) == length(nParamsVec))
    Delta_LL = bsxfun(@minus,LLMat,LLMat(imodelref,:));
    med_LL= median(Delta_LL,2);
    
    figure;
    for imodel = 1:nModels
        
        LLVec = Delta_LL(imodel,:);
        
        CI_LL = sort(median(LLVec(randi(nSubjs.(exptype),nSubjs.(exptype),1000))));
        CI_LL = CI_LL([25 975]);
        fill([0.55 1.45 1.45 0.55]+imodel-1,CI_LL([1 1 2 2]),0.7*ones(1,3)); hold on;
        plot([0.55 1.45]+imodel-1,[med_LL(imodel) med_LL(imodel)],'k-')
    end
    
    hold on;
    violinplot(Delta_LL');
    plot([0 nModels+0.5],[0 0],'k-')
    set(gca,'Xtick',1:nModels,'XTickLabel',modelVec)
    ylabel(sprintf('AICc(model) - AICc(%s)',modelVec{imodelref}))
    title('LL')
    defaultplot
    
else
    AIC = 2*bsxfun(@plus,LLMat,nParamsVec');
    BIC = bsxfun(@plus,2*LLMat,log(n)*nParamsVec');
    AICc = bsxfun(@plus,AIC,2.*bsxfun(@rdivide,nParamsVec'.*(nParamsVec'+1),(bsxfun(@minus,n,nParamsVec')-1)));
    
    Delta_AICc = bsxfun(@minus,AICc,AICc(imodelref,:));
    Delta_BIC = bsxfun(@minus,BIC,BIC(imodelref,:));
    
    med_AICc = median(Delta_AICc,2);
    med_BIC = median(Delta_BIC,2);
    
    figure;
    for imodel = 1:nModels
        
        daicc = Delta_AICc(imodel,:);
        dbic = Delta_BIC(imodel,:);
        
        subplot(1,2,1)
        CI_AICc = sort(median(daicc(randi(nSubjs,nSubjs,1000))));
        CI_AICc = CI_AICc([25 975]);
        fill([0.55 1.45 1.45 0.55]+imodel-1,CI_AICc([1 1 2 2]),0.7*ones(1,3)); hold on;
        plot([0.55 1.45]+imodel-1,[med_AICc(imodel) med_AICc(imodel)],'k-')
        
        subplot(1,2,2)
        CI_BIC = sort(median(dbic(randi(nSubjs,nSubjs,1000))));
        CI_BIC = CI_BIC([25 975]);
        fill([0.55 1.45 1.45 0.55]+imodel-1,CI_BIC([1 1 2 2]),0.7*ones(1,3)); hold on;
        plot([0.55 1.45]+imodel-1,[med_BIC(imodel) med_BIC(imodel)],'k-')
    end
    
    subplot(1,2,1); hold on;
    violinplot(Delta_AICc');
    plot([0 2*nModels+0.5],[0 0],'k-')
    set(gca,'Xtick',1:(2*nModels),'XTickLabel',[modelnameVec modelnameVec])
    ylabel(sprintf('AICc(model) - AICc(%s)',modelnameVec{imodelref}))
    title('AICc')
    defaultplot
    
    subplot(1,2,2); hold on;
    violinplot(Delta_BIC');
    plot([0 2*nModels+0.5],[0 0],'k-')
    set(gca,'Xtick',1:(2*nModels),'XTickLabel',[modelnameVec modelnameVec])
    ylabel(sprintf('BIC(model) - BIC(%s)',modelnameVec{imodelref}))
    title('BIC')
    defaultplot

end







