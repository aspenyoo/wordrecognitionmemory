% load and plot parameter recovery data 
% 
% Aspen Yoo -- February 2, 2016

modelname = 'FPheurs';

switch modelname
    case 'FP'; paramnames = {'M','\sigma','k','d_0'};
    case 'FPheurs'; paramnames = {'M','\sigma','k','d_0'};
    case 'uneqVar'; paramnames = {'\mu_{old}', '\sigma_{old}','k','d_0'};
end

% load ML parameter estimates
load(['paramfit_patternbayes_' modelname '.mat'])
estParam = bestFitParam(15:end,:);
nSubj = size(estParam,1);

% load parameter values for simulated data
load('subjdata.mat')
trueParam = simdata.(modelname).trueparam(15:14+nSubj,:);

% plot parameter recovery! 
bigtitle = [modelname ' parameter recovery'];
paramrecovplot(paramnames, trueParam, estParam, bigtitle)

% plot other plots
nrealSubj = 14;
nSubj = 25;

fakedata.new = simdata.(modelname).nnew(15:nSubj,:)./150;
fakedata.old = simdata.(modelname).nold(15:nSubj,:)./150;
% nnew_part = nan(nSubj-nrealSubj,20); nold_part = nnew_part;
% for isubj = 15:nSubj;
%     [nnew_part(isubj-14,:), nold_part(isubj-14,:)] = simulate_data(modelname, trueParam(isubj-14,:), 1, 30, 20);
% end
% fakedata.new = nnew_part./150;
% fakedata.old = nold_part./150;
pause(0.1);

plotparamfits(modelname, 1, 'patternbayes', estParam, 20, 0, fakedata,[],[1 0 0 0])
% plotparamfits(modelname, islogbinning, optimizationMethod, bestFitParam, nConf, isdatasave, fakedata, subjnum, selectiveplot)