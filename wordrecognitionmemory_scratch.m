% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  RANDOM STUFF
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% october 11, 2016
% weibull
clear all

a = 1;
lambda = 1;
k = .5;

xx = linspace(0,2.5,50);
y = a.*(1-exp(-(xx./lambda).^k));
plot(xx,y);

%%
for i= 1:6;
    subplot(2,6,i)
    xlim([1 20])
    set(gca,'XTick',[1 10 20])
end

%%

for i = 7:12;
    subplot(2,6,i)
    set(gca,'YTick',0:0.05:0.2,'XTick',[1 10 20])
    xlim([1 20])
end

%% checking to see SD of nLLs for different model
% 6/29/2016

modelname = 'REM';
nSamples = 100;

% simulate data
[resp, tp] = simulate_resp(modelname,1,1,[])

paramtry = tp;
nLLVec = nan(1,nSamples);
for isample = 1:nSamples
    isample
    nLLVec(isample) = nLL_approx_vectorized(modelname,paramtry,1,resp.new,resp.old);
end

hist(nLLVec,15); defaultplot
std(nLLVec)

%% looking at diff parameters effect on log function
% july 7, 2016

xx = linspace(-10,10,100);
a = 1;
b = 1;
% c = .05;
nConf = 20;

% yy = a.*log(xx) + b;
signxx = sign(xx);
% yy = a.*log(xx) + b + nConf/2 + 1;
% yy2 = a.*log(xx.*exp(b)./a) + nConf/2 + 1;
% plot(xx,yy,'k',xx,yy2,'r--'); defaultplot
yy = min(max(round(a.*log(abs(xx)) + b + nConf/2 + 1),nConf/2 +1),nConf); % log function shifted up to confidence region
yy2 = min(max(round(a.*log(abs(xx) + b) + nConf/2 + 1),nConf/2 +1),nConf); % log function shifted up to confidence region
subplot(1,2,1);
plot(xx,yy,'k',xx,yy2,'r'); defaultplot
ylim([11 20])

yy(signxx == -1) = nConf+1 - yy(signxx == -1);
yy2(signxx == -1) = nConf+1 - yy2(signxx == -1);
subplot(1,2,2);
plot(xx,yy,'k',xx,yy2,'r'); defaultplot
ylim([ 1 20])
% blahhist = hist(yy,1:nConf)

% =======================================
%     MODEL FITTING RELATED
% =======================================


%% remove txt spacing for models and subjects

modelVec = {'REM','FP'};
binningfnVec = [3 4];

nModels = length(modelVec);
nBinningfns = length(binningfnVec);
subjidVec = [15:36];
nSubj = length(subjidVec);


for imodel = 1:nModels
    model = modelVec{imodel};
    
    for ibinningfn = 1:nBinningfns
        binningfn = binningfnVec(ibinningfn);
        
        disp([model num2str(binningfn)])
        
        for isubj = 1:nSubj
            subjid = subjidVec(isubj)
            removetxtspaces(model,binningfn,subjid);
        end
    end
end


%% separate jobs for each person
clear

modelname = 'REM';
binningfn = 3;
subjidVec = [1:14];
nSubj = length(subjidVec);
MVec = [1:50 55:5:65 70:10:120 200];
% MVec = [1:65 70:5:90 100 110 120];
nMs = length(MVec);
filepath = 'model/4_fitdata/';
approxTime = MVec;%linspace(.22*1000/3600,4.61*1000/3600,nMs);
maxTime = 500;
nJobs = [];
nStartVals = 10;

for isubj = 1:nSubj
    
    subjid = subjidVec(isubj)
    
    jobnumVec = []; estTimeVec = [];
    for iM = 1:nMs
        M = MVec(iM);
        counts = max([nStartVals - countnum2(modelname, binningfn,subjid, [1;M]) 0]);
        jobnumVec = [jobnumVec repmat(M,1,counts)];
        estTimeVec = [estTimeVec repmat(approxTime(iM),1,counts)];
    end
    
    jobfilename = [filepath 'joblist_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt' ];
    create_joblist(jobfilename, jobnumVec, estTimeVec, maxTime, nJobs)
end



%% create joblist

nStartVals = 10;
esttimeVec = linspace(0.08,8,120);
jobnumVec = repmat([55:5:70 80:10:120],1,nStartVals);
esttimeVec = esttimeVec(jobnumVec);
maxTime = 36;
filepath = 'model/4_fitdata/';

% for subjid = 1:14;
%     jobfilename = [filepath 'joblist_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt' ];
%     create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);
% end

jobfilename = [filepath 'joblist_55_120.txt'];
create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%           CHECKING THINGS
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%% checking nLLs are consistent (debugging)
% 08.15.2016

clear

modelname = 'FP';
optimMethod = 'patternbayes';
binningfn = 0;
memstrengthvar = 1;
load(['paramfit_' optimMethod '_' modelname num2str(binningfn) num2str(memstrengthvar) '.mat'])
subjids = 1;

load('subjdata.mat')
nLL = nan(1,length(subjids));
for isubj = subjids
    isubj
    nLL(isubj) = nLL_approx_vectorized( modelname, bestFitParam(isubj,:), binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
    %     bfp = [bestFitParam(isubj,[1 2 end 5 3 4]) 1e-3 1]; % for binningfn = 2, memstrengthvar = 1;
    %     bfp = [bestFitParam(isubj,[1 2 5 4]) 1 0 2 bestFitParam(isubj,3)]; % for binningfn = 1; memstrength = 0
    %     bfp = [bestFitParam(isubj,[1 2 5 4 3]) 0 2 1 ];
    bfp = [bestFitParam(isubj,[1 2 5 4 3]) 0 2 1]; % 01: pcorr --> lin
    memstrengthfn = 1;
    nLL2(isubj) = nLL_approx_vectorized2( modelname, bfp, memstrengthfn, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
end

[nLL_est(subjids) nLL' nLL2'] % nLL3']

%% debugging
isubj = 2;
[pnew pold] = nLL_approx_vectorized( modelname, bestFitParam(isubj,:), binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
figure; plot(1:20,pnew,1:20,pold)
ylim([0 0.3])

%% checking that nLLs are consistnent between old and new nLL functions

clear

modelname = 'REM';
optimMethod = 'patternbayes';
binningfn1 = 1;
memstrengthvar = 0;
binningfn2 = 1;
load(['paramfit_' optimMethod '_' modelname num2str(binningfn2) '.mat'])
load('subjdata.mat')
isubj = 14;

nSamples = 100;
newnLL = nan(1,nSamples);
oldnLL = nan(1,nSamples);
for isample = 1:nSamples;
    isample
    if binningfn1 == 1;
        bfp = [bestFitParam(isubj,:) 0];
        bfp(end-1) = -bfp(end-1);
    end
    newnLL(isample) = nLL_approx_vectorized( modelname, bfp, binningfn1, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
    oldnLL(isample) = nLL_approx_vectorized_old( modelname, bestFitParam(isubj,:), binningfn2, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
end

if nSamples > 10
    figure
    [counts,centers] = hist(newnLL);
    [counts2,centers2] = hist(oldnLL);
    plot(centers,counts,centers2,counts2)
    defaultplot
else
    [newnLL; oldnLL]
end

%% changing file names to be consistent with new model naming scheme
clear all

modelname = 'FPheurs';
binningfnnew = 1;
binningfnold = 1;
memstrengthvar = 0;
nParamsold = 4; % how many parameters does the old way have
filepath = 'model/4_fitdata/BPSfits/';

nSubj = 14;
for isubj = 1:nSubj
    % load old file
    filename = ['paramfit_patternbayes_' modelname num2str(binningfnold) ...
        '_subj' num2str(isubj) '.txt'];
    alldata = dlmread(filename);
    
    % flip d0
    alldata(:,nParamsold) = -alldata(:,nParamsold);
    alldata(:,2*nParamsold+1) = -alldata(:,2*nParamsold+1);
    
    % add 0 sigma_mc
    nrows = size(alldata,1);
    alldata = [alldata(:,1:nParamsold) zeros(nrows,1) ...
        alldata(:,nParamsold+1:2*nParamsold+1) zeros(nrows,1) alldata(:,2*nParamsold+2)];
    
    % write to new file
    newfilename = [filepath 'paramfit_patternbayes_' modelname num2str(binningfnold) ...
        num2str(memstrengthvar) '_subj' num2str(isubj) '.txt'];
    
    fid = fopen(newfilename, 'w');
    formatSpec = repmat('%4.4f \t ',1,2*nParamsold+4);
    formatSpec = [formatSpec(1:end-3) '\r\n'];
    
    for irows = 1:nrows;
        fprintf(fid, formatSpec, alldata(irows,:));
    end
    fclose(fid);
    
end

%% checking nLL of fit parameters with nLL given

bestFitParam(isubj,:)
nLL_est(isubj,:)

[nnew_part, nold_part] = loadsubjdata(isubj,modelname,20);
nLL_approx_vectorized(modelname,bestFitParam(isubj,:),binningfn,nnew_part,nold_part)

%% saving data and model for each subject
clear all

modelname = 'FP';
binningfn = 3;
subjids = 1:14;
nSubj = length(subjids);
optimMethod = 'patternbayes';

load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])

% do this part while on debugger
plotparamfits(modelname,binningfn,optimMethod,bestFitParam(subjids,:),20, 0, 0, subjids, [1 1 0 0])
model = [];
model.new = round(pNew_est*150);
model.old = round(pOld_est*150);
data = [];
data.new = nNew_part;
data.old = nOld_part;
save('datamodel_forRonald.mat','data','model');

%% looking at subject stuff
% debugging for REM subject 4

nSubj = 14;
nnew_part = nan(nSubj,20); nold_part = nnew_part;
for isubj = 1:nSubj
    [nnew_part(isubj,:), nold_part(isubj,:)] = loadsubjdata(isubj);
end

%% Figure 2B from Mickes et al., 2007 paper
% (percent correct for each confidence)

clear all

modelname = 'FP';
binningfn = 3;
optimMethod = 'patternbayes';
nX = 100;
nS = 50;

subjids = [1:10 12:14];
nSubj = length(subjids);

% get data for subjects
load('subjdata.mat')
nNew_part = sum(nNew_part(subjids,:));
nOld_part = sum(nOld_part(subjids,:));

% get number of responses for model
load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'],'bestFitParam') % load bestFitParam
nNew_mod = nan(14,20); nOld_mod = nNew_mod;
for isubj = 1:14;
    subjid = subjids(isubj)
    theta = bestFitParam(subjid,:);
    
    [nNew_mod(isubj,:), nOld_mod(isubj,:)] = simulate_data(modelname, theta, binningfn, nX, nS);
end

% sum of model predictions
nNew_mod_all = sum(nNew_mod(subjids,:));
nOld_mod_all = sum(nOld_mod(subjids,:));

% plot of responses of frequencies
figure; bar([nNew_part; nNew_mod_all; nOld_part; nOld_mod_all]');
legend('data new','model new', 'data old', 'model old')
axis([0.5 20.5 0 250]); defaultplot; defaultplot;
xlabel('Rating')
ylabel('Frequency')

% pc for confidence
PC_part = [nNew_part(1:10)./(nOld_part(1:10) + nNew_part(1:10)) nOld_part(11:20)./(nNew_part(11:20) + nOld_part(11:20))];
PC_mod = [nNew_mod_all(1:10)./(nOld_mod_all(1:10) + nNew_mod_all(1:10)) nOld_mod_all(11:20)./(nNew_mod_all(11:20) + nOld_mod_all(11:20))];
figure; bar([PC_part; PC_mod]'); axis([0.5 20.5 0.4 1]); legend('data','model'); defaultplot;
xlabel('Rating')
ylabel('Proportion Correct')

%% check variance of LL estimates
% clear

nSamples = 10;
modelname = 'REM';
binningfn = 1;
memstrengthvar = 0;
optimMethod = 'patternbayes';
subjid = 2;

% get best fit parameters and subject data
[nnew_part, nold_part] = loadsubjdata(subjid,modelname);
% load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])

nLLVec = nan(1,nSamples); timeVec = nan(1,nSamples);
for isamp = 1:nSamples
    isamp
    t0 = GetSecs;
    nLLVec(isamp) = nLL_approx_vectorized(modelname,bestFitParam(subjid,:),binningfn,memstrengthvar,nnew_part,nold_part,[],30,50);
    timeVec(isamp) = GetSecs - t0;
end

nLLVec

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% model comparison
% 08182016
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all

modelnameVec = {'FP3','FP4','UVSD3','UVSD4','REM3','REM4'};
optimMethod = 'patternbayes';
nModels = length(modelnameVec);
nSubj = 14;
numObs = 300;

[BICMat, AICcMat, nLLMat] = deal(nan(nSubj,nModels));
for imodel = 1:nModels
    modelname = modelnameVec{imodel};
    load(['paramfit_' optimMethod '_' modelname '.mat']) % load bestFitParam
    k = size(bestFitParam,2) - 1;
    
    nLLMat(:,imodel) = nLL_est;
    AICcMat(:,imodel) = 2*(k + nLL_est) + (2*(k+1)*(k+2))/(numObs-k-2);
    BICMat(:,imodel) = 2*nLL_est + k.*log(numObs);
end

% % get means and SEMs of difference

refmodelidx = 2; % which column is the reference model
modelidxVec = 1:nModels;
modelidxVec(refmodelidx) = [];

% 
% [blah, idx] = sortrows(bsxfun(@minus,AICcMat,AICcMat(:,refmodelidx)));
blah = bsxfun(@minus,AICcMat,AICcMat(:,refmodelidx));
idx = 1:14;
mean_AICcMat = mean(blah);
sem_AICcMat = std(blah)./sqrt(size(AICcMat,1));
median_AICcMat = median(blah);
xx = sort(blah);
IQR_AICcMat = [mean(xx(3:4,:)); mean(xx(10:11,:))];

bleh = bsxfun(@minus,BICMat,BICMat(:,refmodelidx));
bleh = bleh(idx,:);
mean_BICMat = mean(bleh);
sem_BICMat = std(bleh)./sqrt(size(BICMat,1));
median_BICMat = median(bleh);
xx = sort(bleh);
IQR_BICMat = [mean(xx(3:4,:)); mean(xx(10:11,:))];

%% plot model comparison

fillamnt = 0.45;

figure, hold on;
for imodel = 1:nModels-1
    currmodidx = modelidxVec(imodel);
    
    fill([imodel-fillamnt imodel+fillamnt imodel+fillamnt imodel-fillamnt],...
        [mean_AICcMat(currmodidx)-sem_AICcMat(currmodidx) mean_AICcMat(currmodidx)-sem_AICcMat(currmodidx) ...
        mean_AICcMat(currmodidx)+sem_AICcMat(currmodidx) mean_AICcMat(currmodidx)+sem_AICcMat(currmodidx)],...
        0.8*ones(1,3),'EdgeColor','none')
end
bar(blah(:,modelidxVec)')
% errorbar(1:size(AICcMat,2),mean_AICcMat,sem_AICcMat);
title('\Delta AICc')
defaultplot
set(gca,'XTick',1:nModels-1,'XTickLabel',modelnameVec(modelidxVec))

figure, hold on;
for imodel = 1:nModels-1
    currmodidx = modelidxVec(imodel);
    
    fill([imodel-fillamnt imodel+fillamnt imodel+fillamnt imodel-fillamnt],...
        [mean_BICMat(currmodidx)-sem_BICMat(currmodidx) mean_BICMat(currmodidx)-sem_BICMat(currmodidx) ...
        mean_BICMat(currmodidx)+sem_BICMat(currmodidx) mean_BICMat(currmodidx)+sem_BICMat(currmodidx)],...
        0.8*ones(1,3),'EdgeColor','none')
end
bar(bleh(:,modelidxVec)')
title('\Delta BIC')
defaultplot
set(gca,'XTick',1:nModels-1,'XTickLabel',modelnameVec(modelidxVec))

%% just averaged bars for AICc BIC

fillamnt = 0.45;
sideamnt = 0.01;
modelnameVec = {'MDG-P','MDG-W','UVSD-P','UVSD-W','REM-P','REM-W'}; 

figure, hold on;
for imodel = 1:nModels-1
    currmodidx = modelidxVec(imodel);
    
    % AICc average
    fill([imodel-fillamnt imodel-sideamnt imodel-sideamnt imodel-fillamnt],...
        [IQR_AICcMat(1,currmodidx) IQR_AICcMat(1,currmodidx)...
        IQR_AICcMat(2,currmodidx) IQR_AICcMat(2,currmodidx)],...
        aspencolors('salmon'),'EdgeColor','none')
    
    % BIC average
    fill([imodel+sideamnt imodel+fillamnt imodel+fillamnt imodel+sideamnt],...
        [IQR_BICMat(1,currmodidx) IQR_BICMat(1,currmodidx)...
        IQR_BICMat(2,currmodidx) IQR_BICMat(2,currmodidx)],...
        aspencolors('booger'),'EdgeColor','none')
    
    % AICc indvl
    xx = imodel - fillamnt/2 - sideamnt;
    plot(xx,blah(:,currmodidx),'.','Color',aspencolors('berry'),'MarkerSize',10)
    
    % BIC indvl
    xx = imodel + fillamnt/2 - sideamnt;
    plot(xx,bleh(:,currmodidx),'.','Color',aspencolors('green'),'MarkerSize',10)
    
    % median lines
    plot([imodel-fillamnt imodel-sideamnt],median_AICcMat(:,currmodidx)*ones(1,2),'Color',aspencolors('berry'))
    plot([imodel+sideamnt imodel+fillamnt],median_BICMat(:,currmodidx)*ones(1,2),'Color',aspencolors('green'))
    
end
plot([0.5 5.5],[0 0],'k--','LineWidth',1)
axis([0.5 5.5 -150 100])
set(gcf, 'Position', [100, 100, 500, 400]);
legend('AICc','BIC')
defaultplot
set(gca,'XTick',1:nModels-1,'XTickLabel',modelnameVec(modelidxVec))
xlabel('model')
ylabel(['model - ' modelnameVec(refmodelidx)]);

%% get pairwise counts of winning models

[counts_AICc, counts_BIC] = deal(nan(nModels));
for imodel = 1:nModels
    Delta_AICc = bsxfun(@minus,AICcMat,AICcMat(:,imodel));
    Delta_BIC = bsxfun(@minus,BICMat,BICMat(:,imodel));
    
    counts_AICc(imodel,:) = sum(sign(Delta_AICc) == -1);
    counts_BIC(imodel,:) = sum(sign(Delta_BIC) == -1);
    
end

%% % % % % % % % % % % % % % %
% UVSD STUFF
% % % % % % % % % % % % % % % %

%% get UVSD binning to work
% 08.26.2016
clear all
theta = [1 1.25 1 0];
binningfn = 0;

[pnew, pold, confBounds] = responses_uneqVar(theta, binningfn);
centers_new = linspace(-3,3,50);
centers_old = linspace(theta(1)-3*theta(2),theta(1)+3*theta(2),50);
counts_new = normpdf(centers_new);
counts_old = normpdf(centers_old,theta(1),theta(2));

figure;
plot(centers_new, counts_new,'Color',aspencolors('greyblue')); hold on;
plot(centers_old, counts_old,'Color',aspencolors('dustygold'))
plot(repmat(confBounds(:),1,2)', repmat([0 max([counts_new(:); counts_old(:)])],length(confBounds),1)','Color',0.7*ones(1,3));

%% fit UVSD model
% 08.26.2016
clear all

modelname = 'UVSD';
binningfn = 4;
switch binningfn
    case 3
        fixparams = [6; 0];
    case 4
        fixparams = [7; 0];
end
nStartVals = 4;

blah = GetSecs;
for isubj = 13:14
    isubj
    fitdata_cluster(isubj, modelname, binningfn, fixparams,[],[],nStartVals);
end
% timee = GetSecs - blah

%% looking at how confhist is calculated. debugging
% 08.29.2016

load('paramfit_patternbayes_FP21.mat')
binningfn = 2;
memstrengthvar = 1;
isubj = 7;

theta = bestFitParam(isubj,:);
[nnew_part,nold_part] = loadsubjdata(isubj);
nLL_approx_vectorized( modelname, theta, binningfn, memstrengthvar, nnew_part, nold_part)

%% =====================================================
%       DOING STUFF WITH FIT PARAMETERS
% ======================================================
%  clear all

modelname = 'REM';
binningfn = 4;
optimMethod = 'patternbayes';
subjids = [1:14];

%% remove txt spacing
for isubj = subjids
    removetxtspaces(modelname,binningfn,isubj);
end

%% get MLE parameter estimates
nStartVals = 10;
getbestfitparams(modelname,binningfn,subjids,nStartVals)

%% load MLE parameter estimates
load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])

%% plot best fit parameters
subjids = [1:14];
plotparamfits(modelname,bestFitParam(subjids,:),binningfn, 20, 0, subjids, [ 1 0 0 0])

%% calculate pnew and pold and save in file for ronald
load('subjdata.mat')

nSubj = 14;
nNew_mod = nan(nSubj,20); nOld_mod = nan(nSubj,20);
for isubj = 1:nSubj
    isubj
    [nNew_mod(isubj,:),nOld_mod(isubj,:)] = nLL_approx_vectorized(modelname,bestFitParam(isubj,:),binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 100, 50);
end

nNew_mod = bsxfun(@times,nNew_mod,sum(nNew_part,2));
nOld_mod = bsxfun(@times,nOld_mod,sum(nOld_part,2));
save(['model/4_fitdata/BPSfits/' modelname num2str(binningfn) num2str(memstrengthvar) '_datamodel.mat'],'nNew_part','nOld_part','nNew_mod','nOld_mod')

%% getting the nLL of MLE estimate for each M


%% plotting binning function

xx = linspace(0,10,100);
nSubj = size(bestFitParam,1);
nConf = 20;

figure;
for isubj = 1:nSubj
    theta = bestFitParam(isubj,:);
    
    switch modelname
        case {'FP','FPheurs'}
            M = theta(1);
            sigma = theta(2);
        case 'REM'
            M = theta(1);                   % number of features
            g = theta(2);                   % probability of success (for geometric distribution)
            ustar = theta(3);               % probability of encoding something
            c = theta(4);                   % probability of encoding correct feature value
            m = theta(5);                   % number of storage attempts
            
            p0 = (1-ustar)^m;               % probability of x_ij= 0
            pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
            pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)
    end
    
    binvalues = 1.5:(nConf/2 - 0.5);
    switch binningfn
        case 0
            slope = theta(end-2);
        case 1
            k = theta(end-2);
        case 2 % logarithmic
            a = theta(end-3);
            b = theta(end-2);
        case 3 % power law
            a = theta(end-4);
            b = theta(end-3);
            gamma = theta(end-2);
            yy = a.*((xx.^gamma - 1)./gamma) + b;
            tempp = gamma.*(binvalues -b)./a + 1;
            tempp(tempp < 0) = nan;
            confbounds = tempp.^(1/gamma);
        case 4 % weibull
            scale = theta(end-5);
            shift = theta(end-4);
            a = theta(end-3);
            b = theta(end-2);
            yy = a.*(1-exp(-(xx./scale).^shift)) + b;
            tempp = 1 - (binvalues -b)./a;
            tempp(tempp < 0) = nan;
            tempp(tempp > 1) = nan;
            confbounds = scale.*(-log(tempp)).^(1/shift);
    end
    d0 = theta(end-1);
    
    subplot(4,4,isubj);
    plot(xx,yy,'k-')
    hold on;
    plot(confbounds+d0,binvalues,'ok')
    defaultplot
    ylim([1 10])
    title(['subj ' num2str(isubj)])
    
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  PARAMETER RECOVERY
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%% SIMULATING DATA WITH COVARIANCE STRUCTURE
clear

modelname = 'REM';
binningfn = 3;
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'])
switch modelname
    case 'FP'
        switch binningfn
            case 3
                d0idx = 6;
            case 4 
                d0idx = 7;
        end
    case 'REM'
        switch binningfn
            case 3
                d0idx = 9;
            case 4 
                d0idx = 10;
        end
end

MU = mean(bestFitParam);
SIGMA = cov(bestFitParam);
nParams = length(MU);

nSubj = 50;
trueparams = [];
while size(trueparams,1)<nSubj
    tempparams = mvnrnd(MU,SIGMA,nSubj-size(trueparams,1));
    
    tempparams(:,1) = round(tempparams(:,1)); % M
    tempparams(tempparams(:,1)>75,:) = []; % deleting M larger than 75
    tempparams(tempparams(:,1)<0,:) = []; % negative M
    tempparams(tempparams(:,end)<0,:) = []; % negative MC noise
    
    switch modelname
        case 'FP'
            tempparams(tempparams(:,2)<0,:) = []; % negative sigma
        case 'REM'
            tempparams(tempparams(:,2)<0,:) = []; % negative g: probability of success (for geometric distribution)
            tempparams(tempparams(:,2)>1,:) = [];
            tempparams(tempparams(:,3)<0,:) = []; % negative ustar: probability of encoding something
            tempparams(tempparams(:,3)>1,:) = [];
            tempparams(tempparams(:,4)<0,:) = []; % negative c: probability of encoding correct feature value
            tempparams(tempparams(:,4)>1,:) = [];
            tempparams(tempparams(:,5)<0,:) = []; % number of storage attempts
            tempparams(:,5) = round(tempparams(:,5)); 
    end
    
    switch binningfn
        case 3
            
        case 4
            tempparams(tempparams(:,end-5)<0,:) = []; % positive scale
    end
    
    trueparams = [trueparams; tempparams];
end

trueparams(:,d0idx) = 0; % d0

% generate sim data
nnew_part = [150 zeros(1,19)];
nold_part = nnew_part;
nNew = nan(nSubj,20); nOld = nNew;
for isubj = 1:nSubj
    isubj
    theta = trueparams(isubj,:);
    [nNew(isubj,:), nOld(isubj,:)] = nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part);
end

% save sim data 
load('model/subjdata.mat')
modname = [modelname num2str(binningfn)];
simdata.(modname).nnew = [nan(14,20); round(150.*nNew)];
simdata.(modname).nold = [nan(14,20); round(150.*nOld)];
simdata.(modname).trueparam = [nan(14,nParams); trueparams];
save('model/subjdata.mat','nNew_part','nOld_part','simdata')


%% plot some simulated subjects

clear all
modelname = 'REM';
binningfn = 3;
nSubj = 25;
modname = [modelname num2str(binningfn)];
load('model/subjdata.mat')

figure;
for isubj = 1:nSubj;
    subplot(5,5,isubj)
    plot(1:20,simdata.(modname).nnew(isubj+14,:)); hold on
    plot(1:20,simdata.(modname).nold(isubj+14,:))
end

%% qualtrends: correlation for one subject

clear all; close all

operationalization = 'median';
nSamples = 10000;

% get theta
modelname = 'FP';
binningfn = 4;
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'],'bestFitParam')

for isubj = 5:14;
% isubj = 4;
isubj
theta = bestFitParam(isubj,:);

[rho, pvalue, slope] = deal(nan(1,nSamples));
for isample = 1:nSamples;
    [~, ~, rho(isample), pvalue(isample), b] = qualtrend_correlation(modelname, theta, binningfn, operationalization);
    slope(isample) = b(2);
end

% figure;
% hist(rho);

rho = sort(rho); 
pvalue = sort(pvalue);

CI = [rho(.05*nSamples) rho(.95*nSamples)]
end
% CI2 = [pvalue(.05*nSamples) pvalue(.95*nSamples)]
% 
% sum(pvalue > .05)/nSamples
% sum(pvalue > .01)/nSamples
% sum(pvalue > .001)/nSamples


% SUBJECTS AND CIs
% median
% 1: [-0.3646 -0.0428]
% 2: [-0.3302 -0.0695]
% 3: [-0.469 -0.208]
% 4: [-0.4841 -0.2278]
% 5: [-0.3675 -0.1095]
% 6: [-0.4751 -0.1993]
% 7: [-0.2656 0.0120]
% 8: [-0.2904 -0.0258]
% 9: [-0.3314 0.0248]
% 10: [-0.4889 -0.2140]
% 11: 
% 12: 
% 13:
% 14: 

%% quualtrends: across subjects

clear all; close all

operationalization = 'min';
modelname = 'FP';
binningfn = 4;
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'],'bestFitParam')

nSubj = 14;
nX = 30;
nS = 20;

[rho, pvalue, slope] = deal(nan(1,nSubj));
for isubj = 1:nSubj;
    theta = bestFitParam(isubj,:);
    [~, ~, rho(isubj), pvalue(isubj), b] = qualtrend_correlation(modelname, theta, binningfn, operationalization, nX, nS);
    slope(isubj) = b(2);
end

mean(rho)
std(rho)/sqrt(nSubj)

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% POOR MAN'S HEIRARCHICAL MODELING: seeing which M fits the best
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear all

mVec = [1:50];
nM = length(mVec);
subjVec = 1:14;
nSubj = length(subjVec);
modelVec = {'FP3','FP4','REM};
nModels = length(modelVec);

bestnLL = struct;
parameterestimates = struct;
for imodel = 1:nModels;
    model = modelVec{imodel};
    
    bestnLL.(model) = nan(nSubj,nM);
    parameterestimates.(model) = cell(1,nM);
    for iM = 1:nM;
        M = mVec(iM)
        
        [bestFitParam, nLL_est] = getbestfitparams(model(1:end-1),str2double(model(end)),subjVec,[],[1; M; M]);
        parameterestimates.(model){iM} = bestFitParam;
        bestnLL.(model)(:,iM) = nLL_est;
    end
   
end


for imodel = 1:nModels;
    model = modelVec{imodel};
    
    figure
    % get posteriors
    postt = -sum(bestnLL.(model));%./(14*300);
%     postt = exp(postt);
%     postt = postt./sum(postt);
    posterior.(model) = postt;
    plot(postt,'k-'); 
    title(model);
    defaultplot
end

%%
subjid = 1;
model = 'FP3';

paramMat = nan(nM,2);
for iM = 1:nM;
    paramMat(iM,:) = parameterestimates.(model){iM}(subjid,1:2);
end
plot(paramMat(:,1),paramMat(:,2),'o')

%% sample and calculate LL of best fit params 

modelname = 'FP';
binningfn = 3;
load('allmodels_fits.mat')
nSamples = 1;
nSubj = 14;
nM = 50;

nLLVec = cell(1,nM);
for isubj = 1:nSubj;
    isubj
    [nnew_part, nold_part] = loadsubjdata(isubj);
    
    for iM = 1:nM;
        iM
        theta = parameterestimates.([modelname num2str(binningfn)]){iM}(isubj,:);
        
%         nLLVec{iM} = nan(nSubjs,nSamples);
        for isample = 1:nSamples
            [ nLLVec{iM}(isubj,isample) ] = nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part);
        end
    end
end

%% doing stuff with the above calculation (previoulsy calculated)
% clear all

% load('FP3_nLLsamples.mat')
nLL = cellfun(@(x) sum(mean(x,2)),nLLVec,'UniformOutput',false);
nLL = cell2mat(nLL);

%% calculate IQRs of parameters for each model
clear

model = 'REM4';
load(['paramfit_patternbayes_' model '.mat'],'bestFitParam')
bestFitParam = sort(bestFitParam);

percentile25 = mean(bestFitParam(3:4,:))
percentile75 = mean(bestFitParam(10:11,:))
median = mean(bestFitParam(7:8,:))
% plot?

%% K-L divergence
clear all

modelVec = [1 5 7];
subjVec = 1:14;

nModels = length(modelVec);
nSubj = length(subjVec);

gamma = 0.577215; % eulers constant
idx = repmat(1:155,2,1);
A = 2./(1:2:309);
G = -gamma - log(2) + cumsum(A);
G = [0 -gamma-log(2) G(idx(:))];

KLMat = nan(nModels,nSubj);
for isubj = 1:nSubj;
    subjid = subjVec{isubj};
    
    % load data
    [data] = concatdata(subjid,'Detection');
    [Xdet] = conditionSeparator(data);
    nCond = length(Xdet);
    
    nBins = 10;
    entr = nan(1,nCond);
    % estimating entropy p*log(p)
    for icond = 1:nCond;
        blah = Xdet{icond};
        
        % grassberger unbiased approximation
        nresp1 = sum(blah(blah(:,1)==0,2));
        ntrials = size(blah(blah(:,1)==0,:),1);
        entr(icond) = ntrials*(G(ntrials) - 1/ntrials*(nresp1*G(nresp1+1) + (ntrials-nresp1)*G(ntrials-nresp1+1)));
        
        %         % biased calculation
        %         presp_nochange = mean(blah(blah(:,1)==0,2));
        %         entr(icond) = presp_nochange*log(presp_nochange) + (1-presp_nochange)*log(1-presp_nochange);
        %         entr(icond) = entr(icond).*size(blah(blah(:,1) == 0,1),1);
        
        changetrials = sortrows(abs(blah(blah(:,1)~=0,:)));
        nchangeTrials = size(changetrials,1);
        nTrialsperBin = floor(nchangeTrials/nBins);
        leftover = nchangeTrials - nTrialsperBin*nBins;
        
        %         % biased
        %         presp_quantile = [mean(changetrials(1:nTrialsperBin+leftover,2)) nan(1,nBins-1)];
        %         for ibin = 2:nBins;
        %             presp_quantile(ibin) = mean(changetrials((ibin-1)*nTrialsperBin+leftover+1:ibin*nTrialsperBin+leftover,2));
        %         end
        %         nTrialsbinVec = [nTrialsperBin+leftover nTrialsperBin*ones(1,nBins-1)];
        %         nTrialsbinVec((presp_quantile == 0) | (presp_quantile == 1)) = [];
        %         presp_quantile((presp_quantile == 0) | (presp_quantile == 1)) = [];
        %
        %         entr(icond) = entr(icond) + nTrialsbinVec*(presp_quantile.*log(presp_quantile) + (1-presp_quantile).*log(1-presp_quantile))';
        
        % using grassberger estimator
        nTrialsbinVec = [nTrialsperBin+leftover nTrialsperBin*ones(1,nBins-1)];
        trialscumsum = [0 cumsum(nTrialsbinVec)];
        for ibin = 1:nBins;
            ntrials = nTrialsbinVec(ibin);
            nresp1 = sum(changetrials((trialscumsum(ibin)+1):trialscumsum(ibin+1),2));
            entrr(ibin) = ntrials*(G(ntrials) - 1/ntrials*(nresp1*G(nresp1+1) + (ntrials-nresp1)*G(ntrials-nresp1+1)));
        end
        entr(icond) = - entr(icond) - sum(entrr);
        
        %     plot(presp_quantile,'ko'); defaultplot; pause
    end
    entr = sum(entr) % mean entropy over conditions
    
    for imodel = 1:nModels;
        model = modelVec(imodel);
        
        % load model parameters
        load(['paramfit_model' num2str(model) '.mat']) % load model
        bleh = cellfun(@(x) strcmp(x,subjid),subjids, 'UniformOutput',false);
        subjidx = find(cell2mat(bleh) == 1);
        theta = bestFitParam(subjidx,:);
        bias = 0;
        
        % calculating log likelihood q*log(p)
        theta(1:nCond) = log(theta(1:nCond));
        if model == 5; theta(nCond+1) = log(theta(nCond+1)); end
        [ll] = loglike(Xdet, model, theta);
        
        KLMat(imodel,isubj) = entr-ll;
    end
    
end
%% seeing which model wins as a function of nTrials (NN data)

modelVec = [1 5 7];
subjVec = 1:5;

nModels = length(modelVec);
nSubjs = length(subjVec);

AICMat = nan(nModels,nSubjs,10);
nParamsVec = nan(1,nModels);
count = zeros(10,3);
for isubj = 1:nSubjs;
    subjnum = subjVec(isubj);
    
    subjids = {'1','2','3','4','5','6','7','8','9','10'};
    subjids = cellfun(@(x) [ num2str(subjnum) '_' x '0000'],subjids,'UniformOutput',false);
    
    for imodel = 1:nModels;
        model = modelVec(imodel);
        [bestFitParam, nLL_est] = combineparamfit(model,subjids);
        AICMat(imodel,isubj,:) = nLL_est + 2*size(bestFitParam,2);
    end
    
    for i = 1:10;
        count(i,:) = count(i,:) + (AICMat(:,isubj,i) == min(AICMat(:,isubj,i)))';
    end
    
    nParamsVec(imodel) = size(bestFitParam,2);
end
%% % % % % % % % % % % % % % % % % % % % % % % 
%             P L O T S ! 
% % % % % % % % % % % % % % % % % % % % % % % 
%% plot LL as a function of M for REM
clear all

load('allmodels_fits.mat')
model = 'FP3';

for isubj = 1:14;
    subplot(4,4,isubj);
    plot(1:50,bestnLL.(model)(isubj,:),'ko')
    defaultplot
end

%% plot 2D likelihood landscape of all parameters for one subject

clear all;

modelname = 'REM';
binningfn = 4;
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'])

switch modelname
    case 'FP'
        %  M, sigma
        plb = [1 1e-3 ];
        pub = [50 3 ];
        paramlabels = {'M','sigma'};
    case 'UVSD'
        % mu, sigma
        plb = [0 1e-3];
        pub = [2 5];
        paramlabels = {'mu','sigma'};
    case 'REM'
        %  M g ustar c m
        plb = [1 1e-3 1e-3 1e-3 1];
        pub = [50 1 1 1 15];
        paramlabels = {'M','g','ustar','c','m'};
end

switch binningfn
    case 3      % power law
        % a, b, gamma
        plb = [plb 0 -10 -5 1e-6];
        pub = [pub 10 10 5 3];
        paramlabels = [paramlabels,'a','b','gamma','sigma_mc'];
    case 4 % weibull binning
        % scale, shift, a, b
        plb = [plb 0 0 0 -3];
        pub = [pub 10 10 3 3];
        paramlabels = [paramlabels,'scale','shift','gamma','sigma_mc'];
end


isubj = 7;
theta = bestFitParam(isubj,:);
theta(end-1) = []; % deleting d0. 
[nnew_part, nold_part] = loadsubjdata(isubj);
nParams = length(plb);
assert(nParams == length(theta));

compMat = nchoosek(1:nParams,2);
nComps = size(compMat,1); % number of pairwise comparisons

nn = 10;
nLLMat = cell(1,nComps);
for icomp = 1:nComps
    icomp
    idx_param1 = compMat(icomp,1);
    idx_param2 = compMat(icomp,2);
    temptheta = theta;
    
    param1Vec = linspace(plb(idx_param1),pub(idx_param1),nn);
    param2Vec = linspace(plb(idx_param2),pub(idx_param2),nn);
    
    if sum(strcmp(modelname,{'FP','REM'})) && idx_param1 == 1; 
        param1Vec = [1 10:5:50];
    end
    
    for ip1 = 1:nn;
        temptheta(idx_param1) = param1Vec(ip1);
        for ip2 = 1:nn;
            temptheta(idx_param2) = param2Vec(ip2);
            
            nLLMat{icomp}(ip1,ip2) = nLL_approx_vectorized( modelname, [temptheta(1:end-1) 0 temptheta(end)], binningfn, nnew_part, nold_part);
        end
    end
    
end

% plot the nLLMats
figure;
clims = [min(min(cell2mat(nLLMat))) max(max(cell2mat(nLLMat)))];
for icomp = 1:nComps
    idx_param1 = compMat(icomp,1);
    idx_param2 = compMat(icomp,2);
    temptheta = theta;
    
    param1Vec = linspace(plb(idx_param1),pub(idx_param1),nn);
    param2Vec = linspace(plb(idx_param2),pub(idx_param2),nn);
    
    if sum(strcmp(modelname,{'FP','REM'})) && idx_param1 == 1; 
        param1Vec = [1 10:10:50];
    end
    subplot(nParams,nParams,nParams*(idx_param1-1)+idx_param2-1); hold off
    imagesc(param1Vec,param2Vec,nLLMat{icomp},clims); hold on
    plot(theta(idx_param1),theta(idx_param2),'r*');
    defaultplot
end
%% plot 2D likelihood landscape of one parameter pair for all subjects

