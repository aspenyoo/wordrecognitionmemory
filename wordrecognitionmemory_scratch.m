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
binningfn = 3;
optimMethod = 'patternbayes';
nSubj = 14;
z = nan(50,2);

for isubj = 1:14;
    removetxtspaces(modelname,binningfn,isubj,optimMethod)
    isubj
    for i=1:50;
        z(i,2) = countnum2(modelname,binningfn,isubj,[1; i]);
        z(i,1)=i;
    end
    z'
end

%% separate jobs for each person

modelname = 'FP';
binningfn = 3;
optimMethod = 'patternbayes';
subjids = 1:14;
nSubj = length(subjids);
Mmax = 50;
filepath = 'model/4_fitdata/';
approxTime = linspace(.1,3,50);
maxTime = 36;
nJobs = [];
nStartVals = 10;

for isubj = 1:nSubj;
    subjid = subjids(isubj);
    
    jobnumVec = []; estTimeVec = [];
    for iM = 1:Mmax;
        counts = nStartVals - countnum2(modelname, binningfn, subjid, [1;iM]);
        jobnumVec = [jobnumVec repmat(iM,1,counts)];
        estTimeVec = [estTimeVec repmat(approxTime(iM),1,counts)];
    end
    
    jobfilename = [filepath 'joblist_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt' ];
    create_joblist(jobfilename, jobnumVec, estTimeVec, maxTime, nJobs)
end

%% checking nLL of fit parameters with nLL given

modelname
isubj = 1;
binningfn = 2;

bestFitParam(isubj,:)
nLL_est(isubj,:)

[nnew_part, nold_part] = loadsubjdata(isubj,modelname,20);
nLL_approx_vectorized(modelname,bestFitParam(isubj,:),binningfn,nnew_part,nold_part)

%% create joblist

nStartVals = 10;
esttimeVec = repmat(linspace(0.08,4.2,50),1,nStartVals);
jobnumVec = repmat(1:50,1,nStartVals);
maxTime = 36;
filepath = 'model/4_fitdata/';

% for subjid = 1:14;
%     jobfilename = [filepath 'joblist_' modelname num2str(binningfn) '_subj' num2str(subjid) '.txt' ];
%     create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);
% end

jobfilename = [filepath 'joblist_08112016.txt'];
create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);

%% get best parameter fits
clear all

modelname = 'FP';
binningfn = 3;
optimMethod = 'patternbayes';
subjids = [1:14];

% for isubj = 1:10
%     removetxtspaces(modelname,binningfn,isubj,optimMethod);
% end

getbestfitparams(modelname,binningfn,1:14)

load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])
subjids = [1:14];
plotparamfits(modelname,binningfn,optimMethod,bestFitParam(subjids,:),20, 0, 0, subjids, [1 1 0 0])

%% checking nLLs are consistent (debugging)

load('subjdata.mat')
for isubj = 1:14;
    isubj
    [ nLL(isubj) ] = nLL_approx_vectorized( modelname, bestFitParam(isubj,:), binningfn, nNew_part(isubj,:), nOld_part(isubj,:), [], 50, 30 );
end


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
optimMethod = 'patternbayes';
subjid = 2;

% get best fit parameters and subject data
[nnew_part, nold_part] = loadsubjdata(subjid,modelname);
% load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])

nLLVec = nan(1,nSamples); timeVec = nan(1,nSamples);
for isamp = 1:nSamples;
    isamp
    t0 = GetSecs;
    nLLVec(isamp) = nLL_approx_vectorized(modelname,bestFitParam(subjid,:),binningfn,nnew_part,nold_part,[],30,50);
    timeVec(isamp) = GetSecs - t0;
end

nLLVec