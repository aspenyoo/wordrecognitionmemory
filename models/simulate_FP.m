function [ pnew, pold ] = responses_FP( theta, islogbinning, cplot, dplot, nS, nX )
% gives proportion of new and old words (pnew and pold) answered in
% particular d given theta values
% [ pnew, pold ] = paramdist( THETA, ISLOGBINNING) where theta is a 4-dim
% vector of parameters. ISLOGBINNING = 1 if log binning and lin if
% ISLOGBINNING = 0
% 
% CPLOT: confidence distribution
% DPLOT: d distribution (memory distribution)
%
% Aspen Yoo -- June 26, 2015


% --------------------------------------------------------------------

if nargin < 4; dplot = 0; end
if nargin < 3; cplot = 0; end

N = 150;            % number of words
M = theta(1);
sigma = theta(2);
if (islogbinning); k = theta(3); else c1 = theta(3); end
if length(theta) < 5;
    L = 10;
    if length(theta) < 4;
        if (islogbinning); d0 = 0; else c2 = 0; end
    else
        if (islogbinning); d0 = theta(4); else c2 = theta(4); end
    end
else
    if (islogbinning); d0 = theta(4); else c2 = theta(4); end
    L = theta(5);
end

if M ~= floor(M)
    assert('M must be a whole number')
end
if nargin < 5;
    nS = 15;           % number of experiments (different Ss)
    nX = 15;          % different X's per S
end


% NEW WORDS p(r_i|theta,new)
% dimensions: N(nWords), M(nFeatures), nRuns (nX per S), nExp(nS)
S = randn(N, M, 1, nS);
t1 = S;
X = bsxfun(@plus,S, sigma.*randn(N, M, nX, nS));
t0 = randn(N, M, 1, nS);

% for looping over trial words (t)  nTrials (N) times.
d_new = nan(N,nX*nS);
d_old = nan(N,nX*nS);

for i = 1:N;
    tempd0 = log(squeeze((1+1/sigma^2)^(M/2).*mean(exp(-.5.*sum(bsxfun(@minus,(1+1/sigma^2)...
        .*bsxfun(@minus, t0(i,:,:,:),X./(1+sigma^2)).^2,t0(i,:,:,:).^2),2)),1)));
    tempd1 = log(squeeze((1+1/sigma^2)^(M/2).*mean(exp(-.5.*sum(bsxfun(@minus,(1+1/sigma^2)...
        .*bsxfun(@minus, t1(i,:,:,:),X./(1+sigma^2)).^2,t1(i,:,:,:).^2),2)),1)));
%     sigfrac = 1 + 1./sigma.^2;
%     tempd0 = log(squeeze(mean(prod(sigfrac.^(1/2).*exp(-1/2.*bsxfun(@minus,sigfrac.* ...
%         bsxfun(@minus, t0(i,:,:,:),(X./(1+sigma.^2))).^2, t0(i,:,:,:).^2)),2),1)));
%     tempd1 = log(squeeze(mean(prod(sigfrac.^(1/2).*exp(-1/2.*bsxfun(@minus,sigfrac.* ...
%         bsxfun(@minus, t1(i,:,:,:),X./(1+sigma.^2)).^2, t1(i,:,:,:).^2)),2),1)));
    d_new(i,:) = tempd0(:);
    d_old(i,:) = tempd1(:);
end
d_new = sort(d_new(:));
d_old = sort(d_old(:));

% memDist --> confdist
if (islogbinning) % log binning
    newHist = round(L.*(2./(1+exp(-(d_new-d0)/k)) - 1)+10.5);
    oldHist = round(L.*(2./((1+exp(-(d_old-d0)/k))) - 1)+10.5);
    newHist(newHist == 21) = 20;
    oldHist(oldHist == 21) = 20;
else % lin binning
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
newHist(newHist==0) = 1e-3;
oldHist(oldHist==0) = 1e-3;

pnew = newHist/sum(newHist);
pold = oldHist/sum(oldHist);

% ===== PLOTS =====

if dplot == 1;
    if cplot ==1;
        subplot(2,1,1)
    end
    binningparameters = theta(3:end);
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