% this code checks if the mirror, variance, and zROC length effect hold for
% a particular model and parameter combination
% 
% encoding/memory strength can be defined either as a decreased sigma or
% increased M and are defined in THETAMAT


modelname = 'FPheurs';
thetaMat = [10 1 4 -5;     % lower encoding strength
            15 1 4 -5];      % higher encoding strength
islogbinning = 1; 
nX = 10; 
nS = 10;
nConf = 20;
N = [150 150];

% nSimilarities = 5;
% similarityVec = linspace(1e-3,2,nSimilarities); 
similarityVec = [0 0.5 1 1.5];
nCond = size(thetaMat,1);

nIter = 1;
nNew = nan(nIter,nConf);
nOld = nan(nIter,nConf);
mean_nnew = nan(nCond,nConf);
mean_nold = nan(nCond,nConf);
for icond = 1:nCond;
    theta = thetaMat(icond,:);
    
    for iiter = 1:nIter;
        [nNew(iiter,:), nOld(iiter,:)] = simulate_data(modelname, theta, islogbinning, nX, nS, nConf, N);
    end
    
    mean_nnew(icond,:) = squeeze(mean(nNew,1))./sum(squeeze(mean(nNew,1)));     % averaging across iterations and normalizing
    mean_nold(icond,:) = squeeze(mean(nOld,1))./sum(squeeze(mean(nOld,1)));     % averaging across iterations and normalizing
end

figure;
clf
newsimcolors = aspencolors(nCond+1,'blue');
oldsimcolors = [0.9 0.9 0.7; aspencolors('dustygold')];
for icond = 1:nCond;
    plot(1:20,mean_nnew(icond,:),'Color',newsimcolors(icond+1,:)); hold on;
    plot(1:20,mean_nold(icond,:),'Color',oldsimcolors(icond,:));
end
hold off
defaultplot
legend('weak new','weak old','strong new','strong old')
xlabel('confidence report')
ylabel('probability density')

% MIRROR EFFECT: FA_SN < FA_WN < H_WO < H_SO
FA_SN = sum(mean_nnew(2,11:20));
FA_WN = sum(mean_nnew(1,11:20));
H_WO = sum(mean_nold(1,11:20));
H_SO = sum(mean_nold(2,11:20));
if (FA_SN < FA_WN) && (FA_WN < H_WO) && (H_WO < H_SO);
    display('the mirror effect holds!')
else
    display('the mirror effect DOES NOT hold')
end

% VARIANCE EFFECT: SO > WO, SN > WN
var_SO = var(1:20,mean_nold(2,:)); % weighting 1-20 by their probability density
var_WO = var(1:20,mean_nold(1,:));
var_SN = var(1:20,mean_nnew(2,:));
var_WN = var(1:20,mean_nnew(1,:));
if (var_SO > var_WO) && (var_SN > var_WN)
    display('variance effect holds!')
else
    display('variance effect DOES NOT hold')
end

% ZROC LENGTH EFFECT: zROCs from the strong condition are shorter than
% zROCs from the weak condition
zSO = norminv(roundn(cumsum(fliplr(mean_nold(2,:))),-5));
zSN = norminv(roundn(cumsum(fliplr(mean_nnew(2,:))),-5));
del_strong = logical(isnan(zSN) + isnan(zSO) + isinf(zSN) + isinf(zSO)); % deleting nans and infs
zSO(del_strong) = []; zSN(del_strong) = [];
vector_strong = [(zSN(end) - zSN(1)), (zSO(end) - zSO(1))];
zROClength_strong = sqrt(sum(vector_strong.^2));

zWO = norminv(roundn(cumsum(fliplr(mean_nold(1,:))),-5));
zWN = norminv(roundn(cumsum(fliplr(mean_nnew(1,:))),-5));
del_weak = logical(isnan(zWN) + isnan(zWO) + isinf(zWN) + isinf(zWO));
zWO(del_weak) = []; zWN(del_weak) = [];
vector_weak = [(zWN(end) - zWN(1)), (zWO(end) - zWO(1))];
zROClength_weak = sqrt(sum(vector_weak.^2));

if zROClength_strong < zROClength_weak;
    display('zROC length effect holds!')
else
    display('zROC length effect DOES NOT hold!')
end

figure;
plot(zSN,zSO,'Color',zeros(1,3));hold on
plot(zWN,zWO,'Color',0.7*ones(1,3)); % plot both zROCs
defaultplot; xlabel('zFA'); ylabel('zHit')
legend('strong condition','weak condition')
title('zROC')
