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

%% debugging paramfit_parforsubj.sh
fixparam = 15;
modelname = 'REM';
nStartVals = 1;

parfor isubj = 8:10;
    isubj
    fitdata_cluster(isubj,modelname,'patternbayes', [1; fixparam],[],[],nStartVals);
end
exit;

%% checking which jobs need to be resubmitted

modelname = 'FP';
binningfn = 5;
optimMethod = 'patternbayes';
nSubj = 14;
z = nan(50,2);

for isubj = 1:nSubj;
    removetxtspaces(modelname,binningfn,isubj,optimMethod)
    isubj
    for i=1:50;
        z(i,2) = countnum2(modelname,binningfn,isubj,[1; i]);
        z(i,1)=i;
    end
    z'
end

%% separate jobs for each person
clear 

modelname = 'FP';
binningfn = 2;
memstrengthvar = 1;
optimMethod = 'patternbayes';
subjids = 14;
Mmax = 50;
filepath = 'model/4_fitdata/';
approxTime = linspace(.22*1000/3600,4.61*1000/3600,50);
maxTime = 8;
nJobs = [];
nStartVals = 10;

for isubj = subjids;
    subjid = isubj;%subjids(isubj);
    
    jobnumVec = []; estTimeVec = [];
    for iM = 1:Mmax;
        counts = max([nStartVals - countnum2(modelname, binningfn,memstrengthvar,subjid, [1;iM]) 0]);
        jobnumVec = [jobnumVec repmat(iM,1,counts)];
        estTimeVec = [estTimeVec repmat(approxTime(iM),1,counts)];
    end
    
    jobfilename = [filepath 'joblist_' modelname num2str(binningfn) num2str(memstrengthvar) '_subj' num2str(subjid) '.txt' ];
    create_joblist(jobfilename, jobnumVec, estTimeVec, maxTime, nJobs)
end



%% create joblist

nStartVals = 10;
esttimeVec = repmat([0.08 0.164 0.25 0.33 .42 0.58 0.84 1.26 1.68 2.10 2.52 3.36 4.2],1,nStartVals);
jobnumVec = repmat([1:5 7 10 15 20 25 30 40 50],1,nStartVals);
maxTime = 12;
filepath = 'model/4_fitdata/';

% for subjid = 1:14;
%     jobfilename = [filepath 'joblist_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt' ];
%     create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);
% end

jobfilename = [filepath 'joblist_08302016.txt'];
create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);

%% get best parameter fits
clear all

modelname = 'FP';
binningfn = 1;
memstrengthvar = 0;
optimMethod = 'patternbayes';
subjids = [1:14];

for isubj = subjids;
    removetxtspaces(modelname,binningfn,memstrengthvar,isubj,optimMethod);
end

getbestfitparams(modelname,binningfn,memstrengthvar,subjids)

%%
load(['paramfit_' optimMethod '_' modelname num2str(binningfn) num2str(memstrengthvar) '.mat'])
subjids = [1:14];
plotparamfits(modelname,binningfn,memstrengthvar,optimMethod,bestFitParam(subjids,:),20, 0, 0, subjids, [0 0 1 0])

%% checking nLLs are consistent (debugging)
% 08.15.2016

clear

modelname = 'REM';
optimMethod = 'patternbayes';
binningfn = 1;
memstrengthvar = 0;
load(['paramfit_' optimMethod '_' modelname num2str(binningfn) num2str(memstrengthvar) '.mat'])
binningfn = 2;
memstrengthvar = 1;
subjids = 1:14;

load('subjdata.mat')
nLL = nan(1,length(subjids));
for isubj = subjids;
    isubj
    nLL(isubj) = nLL_approx_vectorized( modelname, bestFitParam(isubj,:), binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
%     nLL2(isubj) = nLL_approx_vectorized( modelname, [bestFitParam(isubj,1:end-1) 0], binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
%     nLL3(isubj) = nLL_approx_vectorized( modelname, [bestFitParam(isubj,1:end-1) 1e-3], binningfn, memstrengthvar, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
%     nLL2(isubj) = nLL_approx_vectorized_old( modelname, bestFitParam(isubj,1:end-1), binningfn, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
end

[nLL_est(subjids) nLL']% nLL2' nLL3']

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

if nSamples > 10;
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
for isubj = 1:nSubj;
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
for isubj = 1:nSubj;
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
for isamp = 1:nSamples;
    isamp
    t0 = GetSecs;
    nLLVec(isamp) = nLL_approx_vectorized(modelname,bestFitParam(subjid,:),binningfn,memstrengthvar,nnew_part,nold_part,[],30,50);
    timeVec(isamp) = GetSecs - t0;
end

nLLVec

%% model comparison (AIC)
% 08182016

clear all

modelnames = {'FP','FP'};
binningfns = [3 5];
optimMethod = 'patternbayes';
nModels = length(modelnames);
nSubj = 14;

AICVec = nan(nSubj,nModels);
nLLMat = nan(nSubj,nModels);
for imodel = 1:nModels;
    modelname = modelnames{imodel};
    binningfn = binningfns(imodel);
    load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat']) % load bestFitParam
    k = size(bestFitParam,2) -1
    
    nLLMat(:,imodel) = nLL_est;
    AICVec(:,imodel) = 2*(k + nLL_est);
end

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
binningfn = 1;
optimMethod = 'patternbayes';
fixparams = [4; 0];
memstrengthvar = [];

blah = GetSecs;
for isubj = 1:14;
    isubj
    fitdata_cluster(isubj, modelname, binningfn, memstrengthvar, optimMethod, fixparams,[],[],3);
end
timee = GetSecs - blah

%% looking at how confhist is calculated. debugging
% 08.29.2016

load('paramfit_patternbayes_FP21.mat')
binningfn = 2;
memstrengthvar = 1;
isubj = 7;

theta = bestFitParam(isubj,:);
[nnew_part,nold_part] = loadsubjdata(isubj);
nLL_approx_vectorized( modelname, theta, binningfn, memstrengthvar, nnew_part, nold_part)