



%% plot data and model ROCs
clear all; close all

modelnameVec = {'UVSD'};
binningfnVec = [3];
YN_indvlplots = 1;
YN_aveplot = 0;
subjVec = 1:14;

nModels = length(modelnameVec);
nBinningfns = length(binningfnVec);
nSubj = length(subjVec);

% get human data
load('subjdata.mat')
data_pNew = bsxfun(@rdivide,nNew_part,sum(nNew_part,2));
data_pOld = bsxfun(@rdivide,nOld_part,sum(nOld_part,2));


%     figure(1); hold on;
    
    data_cumHit = norminv(cumsum(flipud(data_pOld'))');
    data_cumFA = norminv(cumsum(flipud(data_pNew'))');
    
    mean_data_cumHit = mean(data_cumHit);
    mean_data_cumFA = mean(data_cumFA);
    sem_data_cumHit = std(data_cumHit)./nSubj;
    sem_data_cumFA = std(data_cumFA)./nSubj;
    
    errorbarxy(mean_data_cumFA,mean_data_cumHit,...
        sem_data_cumFA,sem_data_cumHit,{'k','k','k'})


figureidx = 1;
for imodel = 1:nModels
    modelname = modelnameVec{imodel};
    
    for ibinningfn = 1:nBinningfns
        binningfn = binningfnVec(ibinningfn);
        
        figureidx = figureidx + 1;
        
        % load in all best fit parameters
        load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'])
        
        % for each subject
        [model_pNew, model_pOld, data_pNew, data_pOld] = deal(nan(nSubj,20));
        for isubj = 1:nSubj
            subjid = subjVec(isubj)
            
            theta = bestFitParam(subjid,:);

            % calculate pnew, pold
            [ model_pnew, model_pold ] = nLL_approx_vectorized( modelname, theta, binningfn, [150 zeros(1,19)], [150 zeros(1,19)]);
            model_pNew(isubj,:) = model_pnew;
            model_pOld(isubj,:) = model_pold;
            
            % calculate ROC
            model_cumfa = norminv(cumsum(flipud(model_pnew(:))));
            model_cumhit = norminv(cumsum(flipud(model_pold(:))));
            
            if (YN_indvlplots)
                % plot individual ROC plots
                figure(figureidx); hold on;
                subplot(4,4,isubj)
                plot(model_cumfa,model_cumhit,'-k'), hold on
                plot(data_cumFA(isubj,:),data_cumHit(isubj,:),'ok')
%                 axis([0 1 0 1])
                defaultplot;
                axis square
            else
                
            end
        end
        
        if (YN_aveplot)
            figure(1); hold on;
            
            % calculate ROC for all data
            model_cumHit = cumsum(flipud(model_pOld'))';
            model_cumFA = cumsum(flipud(model_pNew'))';
            
            % get mean and sem of model
            mean_model_cumHit = mean(model_cumHit);
            mean_model_cumFA = mean(model_cumFA);
            sem_model_cumHit = std(model_cumHit)./nSubj;
            sem_model_cumFA = std(model_cumFA)./nSubj;
            
            % plots ROC
            errorbarxy(mean_model_cumFA,mean_model_cumHit,...
                sem_model_cumFA,sem_model_cumHit)
            defaultplot
            axis([0 1 0 1])
        end
        
    end
end

%% plot group fits for all models
clear all

% fullmodelnameVec = {'FP3','FP4','UVSD3','UVSD4','REM3','REM4'};
fullmodelnameVec = {'FP3','UVSD3','REM3','FP4','UVSD4','REM4'};
% fullmodelnameVec = {'UVSD3','UVSD4'};
subjids = [1:14];
nModels = length(fullmodelnameVec);

for imodel = 1:nModels
    fullmodelname = fullmodelnameVec{imodel};
    
    % load in information
    load(['paramfit_patternbayes_' fullmodelname '.mat'])
    
    subplot(2,3,imodel)
    plotparamfits(modelname,bestFitParam(subjids,:),str2double(fullmodelname(end)), 20, 0, subjids, [ 0 1 0 0])
    title(fullmodelname)
end

%% indvl and group fit for all models
clear all; close all

fullmodelnameVec = {'FP3','FP4','UVSD3','UVSD4','REM3','REM4'};
modelnameVec = {'MDG-P','MDG-W','UVSD-P','UVSD-W','REM-P','REM-W'};
subjids = [1:14];
nModels = length(fullmodelnameVec);

figure
for imodel = 1:nModels
    fullmodelname = fullmodelnameVec{imodel};
    
    % load in information
    load(['paramfit_patternbayes_' fullmodelname '.mat'])
    
    subplot(2,6,6+imodel)
    plotparamfits(modelname,bestFitParam(subjids,:),str2double(fullmodelname(end)), 20, 0, subjids, [ 0 1 0 0])
end


subjid = 4;
for imodel = 1:nModels;
    fullmodelname = fullmodelnameVec{imodel};
    
    % load in information
    load(['paramfit_patternbayes_' fullmodelname '.mat'])
    
    % plot
    subplot(2,6,imodel);
    plotparamfits(modelname,bestFitParam(subjid,:),str2double(fullmodelname(end)),20,[],subjid,[1 0 0 0])
    subplot(2,6,imodel);
    ylim([0 0.4])
    title(modelnameVec{imodel})
end


%% plot fig plots for each subject
subjids = [2 5 7:14];
for isubj = 1:length(subjids);
    subjid = subjids(isubj);
    
    fig_qualtrends(subjid)
    savefig(['qualtrends_subj' num2str(subjid) '.fig'])
end