function [ pnew, pold ] = responses_FP_sampling( theta, islogbinning, cplot, dplot, nX, nS )
% gives proportion of new and old words (pnew and pold) answered in
% particular d given theta values
% [ pnew, pold ] = paramdist( THETA, ISLOGBINNING) where theta is a 4-dim
% vector of parameters. ISLOGBINNING = 1 if log binning and lin if
% ISLOGBINNING = 0
%
% CPLOT: confidence distribution
% DPLOT: d distribution (memory distribution)
%
% Aspen Yoo -- Sept 7, 2015
% BAS'S EDIT!


% --------------------------------------------------------------------

if nargin < 4; dplot = 0; end
if nargin < 3; cplot = 0; end
if nargin < 6; nS = 10000; end
if nargin < 5; nX = 1; end

N = 150;
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

d_new = nan(N,nX,nS);
d_old = nan(N,nX,nS);
for iX = 1:nX;
    %     iruns
    X = randn(N, M)*sqrt(sigma^2+1);
    
    for iS = 1:nS;
        sNew = randn(N,M);
        sOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(N,M)*(1/sqrt(1/sigma^2 + 1));
        for i = 1:N;
            % new trials
            snew = sNew(i,:);
            d_new(i,iX,iS) = M/2*log(1+1/sigma^2) + 0.5*snew*snew' + log(mean(exp(-0.5*((1/sigma.^2+1)*sum((repmat(snew,[N,1])-X/(1+sigma^2)).^2,2)))));
            
            % old trials
            sold = sOld(i,:);
            d_old(i,iX,iS) = M/2*log(1+1/sigma^2) + 0.5*sold*sold' + log(mean(exp(-0.5*((1/sigma.^2+1)*sum((repmat(sold,[N,1])-X/(1+sigma^2)).^2,2)))));
        end
    end
end

% BINNING. p(conf|X0,C)
if (islogbinning) % log binning
    newHist = round(L.*(2./(1+exp(-(d_new-d0)./k)) - 1)+10.5);
    oldHist = round(L.*(2./((1+exp(-(d_old-d0)./k))) - 1)+10.5);
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

newHist = reshape(newHist,[20, nX*nS]);
oldHist = reshape(oldHist,[20, nX*nS]);

pnew = bsxfun(@rdivide,newHist,sum(newHist,1));
pold = bsxfun(@rdivide,oldHist,sum(oldHist,1));

% ===== PLOTS =====

if dplot == 1;
    if cplot ==1;
        subplot(2,1,1)
    end
    binningparameters = theta(3:end);
    
    d_old = reshape(d_old,[150,nX*nS]);
    d_new = reshape(d_new,[150,nX*nS]);
    d_old = d_old'; d_new = d_new';
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