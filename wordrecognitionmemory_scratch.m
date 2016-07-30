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

modelname = 'REM';
binningfn = 1;
optimMethod = 'patternbayes';
nSubj = 14;
z = nan(50,2);

for isubj = 5:nSubj;
    removetxtspaces(modelname,binningfn,isubj,optimMethod)
    isubj
    for i=1:50;
        z(i,2) = countnum2(modelname,binningfn,isubj,[1; i]);
        z(i,1)=i;
    end
    z'
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
esttimeVec = repmat(linspace(.3,4,50),1,nStartVals);
jobnumVec = repmat(1:50,1,nStartVals);
maxTime = 36;

filepath = 'model/4_fitdata/';
jobfilename = [filepath 'joblist_07252016.txt'];
create_joblist(jobfilename, jobnumVec, esttimeVec, maxTime);

%% get best parameter fits

modelname = 'FP';
binningfn = 2;
subjids = [1:14];
optimMethod = 'patternbayes';

nSubj = length(subjids);
for isubj = 1:nSubj;
    subjid = subjids(isubj);
    removetxtspaces(modelname,binningfn,subjid,optimMethod);
end

getbestfitparams(modelname,binningfn,subjids)

load(['paramfit_' optimMethod '_' modelname num2str(binningfn) '.mat'])
plotparamfits(modelname,binningfn,optimMethod,bestFitParam,20, 0, 0, subjids, [1 1 0 0])

%% looking at subject stuff
% debugging for REM subject 4

nSubj = 14;
nnew_part = nan(nSubj,20); nold_part = nnew_part;
for isubj = 1:nSubj;
   [nnew_part(isubj,:), nold_part(isubj,:)] = loadsubjdata(isubj); 
end

