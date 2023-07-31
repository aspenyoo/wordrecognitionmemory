function [ pnew, pold ] = simulate_VPheurs( theta, islogbinning, cplot, dplot )
% gives proportion of new and old words (pnew and pold) answered in 
% particular d given theta values
% 
% Aspen Yoo -- march 6, 2015
% also paramdistlogheursVP.m

% close all; clear all;
% theta = [1.8 .001 30 10 8];
% pplot = 0;
% rplot = 1;

% --------------------------------------------------------------------
if nargin < 3; cplot = 0; end
if nargin < 4; dplot = 0; end

% random number generator
% rng('shuffle');

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
X = bsxfun(@plus,S,sigma.*randn(N, M, nRuns, nExp));
t0 = randn(N, M, 1, nExp);

% calculating distance between t's and X
d_old = nan(N, nRuns, nExp);
d_new = nan(N, nRuns, nExp);
for i = 1:N;
    dist1 = sqrt(sum(bsxfun(@minus,t1(i,:,:,:),X).^2,2));
    dist0 = sqrt(sum(bsxfun(@minus,t0(i,:,:,:),X).^2,2));
    
    d_old(i,:,:,:) = (min(dist1));
    d_new(i,:,:,:) = (min(dist0));
end
d_old = -d_old(:);
d_new = -d_new(:);

% memDist --> confdist
if (islogbinning) 
newHist = round(L.*(2./(1+exp(-(d_new-d0)/k)) - 1)+10.5);
oldHist = round(L.*(2./((1+exp(-(d_old-d0)/k))) - 1)+10.5);
newHist(newHist == 21) = 20;
oldHist(oldHist == 21) = 20;
else
    m = 18/(c2-c1);
    b = 1.5-m*c1;

    newHist = round(m.*d_new + b);
    oldHist = round(m.*d_old + b);
    
    newHist(newHist < 1) = 1;
    newHist(newHist > 20) = 20;
    oldHist(oldHist < 1) = 1;
    oldHist(oldHist > 20) = 20;
end
newHist = histc(newHist,1:20);
oldHist = histc(oldHist,1:20);

% changing any 0 freq to 1, (prevents LL from going to -Inf)
newHist(newHist==0) = 1;
oldHist(oldHist==0) = 1;

% pnew = newHist/sum(newHist);
% pold = oldHist/sum(oldHist);
pnew = newHist/sum(newHist);
pold = oldHist/sum(oldHist);

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