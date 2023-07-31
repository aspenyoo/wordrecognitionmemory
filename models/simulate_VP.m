function [ pnew, pold ] = simulate_VP( theta, islogbinning, cplot, dplot )
% gives proportion of new and old words (pnew and pold) answered in 
% particular d given theta values
% 
% also paramdistlogVP.m
% Aspen Yoo
% march 6, 2015
% 
% updated June 26, 2015 to include linear binning

% theta = [0.8026    1.4474   10.0000    1.4737];
% pplot = 1;
% rplot = 1;

% --------------------------------------------------------------------
if nargin < 3;
    cplot = 0;
    dplot = 0;
end

M = roundn(theta(1),-4);
Jbar = theta(2);
tau = theta(3);
if (islogbinning); k = theta(4); else c1 = theta(4); end
if length(theta) < 6;
    L = 10;
    if length(theta) < 5;
        if (islogbinning); d0 = 0; else c2 = 0; end
    else
        if (islogbinning); d0 = theta(5); else c2 = theta(5); end
    end
else
    if (islogbinning); d0 = theta(5); else c2 = theta(5); end
    L = theta(6);
end

nExp = 15;           % number of experiments (different Ss)
nRuns = 15;          % different X's per S
N = 150;            % number of words

% NEW WORDS p(r_i|theta,new)
% dimensions: N(nWords), M(nFeatures), nRuns (nX per S), nExp(nS)
S = randn(N, M, 1, nExp);
t1 = S;
sigma = sqrt(1./(gamrnd(Jbar/tau,tau,[N,M,nRuns,nExp])));
sigma(isinf(sigma)) = 1e3;

X = bsxfun(@plus, S, sigma.*randn(N, M, nRuns, nExp));
t0 = randn(N, M, 1, nExp);
sigfrac = 1 + 1./sigma.^2;
logsigfrac = log(sigfrac);

% for looping over trial words (t)  nTrials (N) times.
d_new = nan(N,nRuns*nExp);
d_old = nan(N,nRuns*nExp);

% if 0
% t0 = randn(N, M, 1, nExp);
% for i = 1:N;
%     tempd0 = log(squeeze(mean(exp(sum(0.5*logsigfrac-1/2.*bsxfun(@minus,sigfrac.* ...
%         bsxfun(@minus, t0(i,:,:,:),(X./(1+sigma.^2))).^2, t0(i,:,:,:).^2),2)),1)));
%     tempd1 = log(squeeze(mean(exp(sum(0.5*logsigfrac-1/2.*bsxfun(@minus,sigfrac.* ...
%         bsxfun(@minus, t1(i,:,:,:),(X./(1+sigma.^2))).^2, t1(i,:,:,:).^2),2)),1)));
%     memDist_new(i,:) = tempd0(:);
%     memDist_old(i,:) = tempd1(:);
% end
% else
Xtemp = X./(1+sigma.^2);
sumt0sq = 0.5*sum(t0.^2,2);
sumt1sq = 0.5*sum(t1.^2,2);

logsigfracsum = sum(0.5*logsigfrac,2);
for i = 1:N;
    tempd0 = log(squeeze(mean(exp(logsigfracsum - bsxfun(@minus, 0.5*sum(sigfrac.* ...
        bsxfun(@minus, t0(i,:,:,:), Xtemp).^2,2), sumt0sq(i,:,:,:))),1)));
    tempd1 = log(squeeze(mean(exp(logsigfracsum - bsxfun(@minus, 0.5*sum(sigfrac.* ...
        bsxfun(@minus, t1(i,:,:,:), Xtemp).^2,2), sumt1sq(i,:,:,:))),1)));
    
    %tempd0 = log(squeeze(mean(exp(logsigfracsum - 0.5*sum(bsxfun(@minus,sigfrac.* ...
    %    bsxfun(@minus, t0(i,:,:,:), Xtemp).^2, t0sq(i,:,:,:)),2)),1)));
    %tempd1 = log(squeeze(mean(exp(logsigfracsum - 0.5*sum(bsxfun(@minus,sigfrac.* ...
    %    bsxfun(@minus, t1(i,:,:,:), Xtemp).^2, t1sq(i,:,:,:)),2)),1)));
    d_new(i,:) = tempd0(:);
    d_old(i,:) = tempd1(:);
end
% end

d_new = sort(d_new(:));
d_old = sort(d_old(:));

% memDist --> confdist
if (islogbinning) 
confDist_new = round(L.*(2./(1+exp(-(d_new-d0)/k)) - 1)+10.5);
confDist_old = round(L.*(2./((1+exp(-(d_old-d0)/k))) - 1)+10.5);
confDist_new(confDist_new == 21) = 20;
confDist_old(confDist_old == 21) = 20;
else
    m = 18/(c2-c1);
    b = 1.5-m*c1;

    confDist_new = round(m.*d_new + b);
    confDist_old = round(m.*d_old + b);
    
    confDist_new(confDist_new < 1) = 1;
    confDist_new(confDist_new > 20) = 20;
    confDist_old(confDist_old < 1) = 1;
    confDist_old(confDist_old > 20) = 20;
end

confDist_new = histc(confDist_new,1:20);
confDist_old = histc(confDist_old,1:20);

% changing any 0 freq to 1, (prevents LL from going to -Inf)
confDist_new(confDist_new==0) = 1;
confDist_old(confDist_old==0) = 1;

pnew = confDist_new/sum(confDist_new);
pold = confDist_old/sum(confDist_old);

% ===== PLOTS =====

if dplot == 1;
    if cplot ==1;
        subplot(2,1,1)
    end
    binningparameters = theta(4:end);
    memdistplot(d_old, d_new, binningparameters, islogbinning);
end

if cplot == 1;
    if dplot == 1;
        subplot(2,1,2)
    else
        figure;
    end
    confdistplot(pnew, pold)
end