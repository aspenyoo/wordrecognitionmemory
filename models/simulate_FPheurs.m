function [ pnew, pold ] = simulate_FPheurs( theta, islogbinning, cplot, dplot )
% gives proportion of new and old words (pnew and pold) answered in
% particular d given theta values
% [ pnew, pold ] = paramdist( theta ), where theta is the 3-dimensional
% vector [ M sigma k ].
%
% Aspen Yoo -- 12/8/14


% --------------------------------------------------------------------

if nargin < 4; dplot = 0; end
if nargin < 3; cplot = 0; end


M = theta(1);
sigma = theta(2);
if (islogbinning); k = theta(3); else c1 = theta(3); end
if length(theta) < 5
    L = 10;
    if length(theta) < 4
        if (islogbinning); d0 = 0; else c2 = 0; end
    else
        if (islogbinning); d0 = theta(4); else c2 = theta(4); end
    end
else
    if (islogbinning); d0 = theta(4); else c2 = theta(4); end
    L = theta(5);
end

assert(M == floor(M),'error: M must be a whole number')

nExp = 15;           % number of experiments (different Ss)
nRuns = 15;          % different X's per S
N = 150;            % number of words

% NEW WORDS p(r_i|theta,new)
% dimensions: N(nWords), M(nFeatures), nRuns (nX per S), nExp(nS)
S = repmat(randn(N, M, 1, nExp),[1, 1, nRuns, 1]);
t1 = S;
X = S + sigma.*randn(N, M, nRuns, nExp);
t0 = repmat(randn(N, M, 1, nExp), [1 1 nRuns 1]);

% for looping over trial words (t)  nTrials (N) times.
d_new = nan(N,nRuns,nExp);
d_old = nan(N,nRuns,nExp);
for i = 1:N
    dist1 = sqrt(sum(bsxfun(@minus,t1(i,:,:,:),X).^2,2));
    dist0 = sqrt(sum(bsxfun(@minus,t0(i,:,:,:),X).^2,2));
    
    d_old(i,:,:) = squeeze(min(dist1));
    d_new(i,:,:) = squeeze(min(dist0));
end
d_new = -d_new(:);
d_old = -d_old(:);

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
newHist(newHist==0) = 1;
oldHist(oldHist==0) = 1;

pnew = newHist/sum(newHist);
pold = oldHist/sum(oldHist);

% ===== PLOTS =====

if dplot == 1
    if cplot ==1
        subplot(2,1,1)
    end
    binningparameters = theta(3:end);
    memdistplot(d_old, d_new, binningparameters, islogbinning);
end

if cplot == 1
    if dplot == 1
        subplot(2,1,2)
    else
        figure;
    end
    confdistplot(pnew, pold)
end