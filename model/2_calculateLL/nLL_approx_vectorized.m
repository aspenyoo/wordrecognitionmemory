function [ nLL ] = nLL_approx_vectorized( modelname, theta, islogbinning, nnew_part, nold_part, fixparams, nX, nS, nConf )
% nLL_approx calculates the negative log likelihood using an approximation
% method
%
% [nLL] = nLL_approx_vectorized(MODELNAME, THETA, ISLOGBINNING, NNEW_PART,
% NOLD_PART) calculates the negative log p(nnew_part,nold_part|theta,model)
%
% ===== INPUT VARIABLES =====
% MODELNAME: 'FP','FPheurs','VP','VPheurs','uneqVar'
% THETA: parameter values
% ISLOGBINNING: 1: logistic. 0: linear
% NNEW_PART: 1x20 vector of responses for new distribution (total 150)
% NOLD_PART: 1x20 vector of responses for old distribution (total 150)
% FIXPARAMS: a 2xn matrix in which first row corresponds to index of which
% parameter to fix and second row corresponds to the value of that
% parameter.
% NX: how many Xs (noisy memories) are drawn from a single S (word matrix)
% Ns: how many Ss (word matrices) are drawn
% NCONF: how many confidence values to have
%
% ===== OUTPUT VARIABLES =====
% NLL: negative log likelihood
%
% Aspen Yoo - January 28, 2016

if nargin < 6; fixparams = []; end
if nargin < 7; nX = 30; end
if nargin < 8; nS = 20; end
if nargin < 9; nConf = 20; end

rng('shuffle')

% if theta has fixed parameters, adjust accordingly
if ~isempty(fixparams);
    nParams = length(theta) + size(fixparams,2);
    tempTheta = nan(1,nParams);
    unfixedparams = 1:nParams;
    unfixedparams(fixparams(1,:)) = [];
    tempTheta(fixparams(1,:)) = fixparams(2,:);
    tempTheta(unfixedparams) = theta;
    theta = tempTheta;
end


if strcmp('uneqVar',modelname)
    [pnew, pold] = responses_uneqVar(theta, islogbinning);
    nLL = -sum(log(pnew).*nnew_part) - sum(log(pold).*nold_part);
else
    M = theta(1);
    sigma = theta(2);
    if (islogbinning); k = theta(3); else c1 = theta(3); end
    if length(theta) < 4;
        if (islogbinning); d0 = 0; else c2 = 0; end
    else
        if (islogbinning); d0 = theta(4); else c2 = theta(4); end
    end
    Nold = sum(nold_part); Nnew = sum(nnew_part);
    L = nConf/2;
    J = 1/sigma.^2+1;
    lambda = 0.01;             % lapse rate
    
    if M ~= floor(M)
        M = floor(M);
        %     assert('M must be a whole number')
    end
    
    if ~(islogbinning)
        m =(nConf-2)/(c2-c1);       % slope
        b = 1.5-m*c1;               % y-intercept
    end
    
    d_old = nan(Nold*nS,nX);
    for iX = 1:nX;
        
        X = randn(Nold, M)*sqrt(sigma^2+1);
        Xrep = repmat(X,[nS 1]);
        
        SNew = randn(Nnew*nS,M);
        SOld = (Xrep/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
        
        switch modelname
            case 'FP'
                d_new = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                    sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2)))));
                d_old(:,iX) = M/2*log(1+1/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                    sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold*nS])).^2,2)))));
            case 'FPheurs'
                d_new = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew*nS])).^2,2),[],1));
                d_old(:,iX) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold*nS])).^2,2),[],1));
        end
        
        % binning new words.
        if (islogbinning)
            newHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf); % bounds: [1 20]
        else
            newHist = min(max(round(m.*d_new(:) + b),1),20); % bounds: [1 20]
        end
        newHist = histc(newHist,1:nConf);
        pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
        LL_new(iX) = nnew_part*log(pnew);
        
    end
    
    % BINNING OLD WORDS. p(conf|X0,C)
    if (islogbinning) % log binning
        oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
    else % lin binning
        oldHist = max(min(round(m.*d_old(:) + b),nConf),1);
    end
    oldHist = histc(oldHist,1:nConf); % histogram
    % oldHist(oldHist==0) = 1e-3; % changing any 0 freq to 1, (prevents LL from going to -Inf)
    pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
    
    % calculating nLL
    LL_old = nold_part*log(pold);
    LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
    nLL = -LL_new-LL_old;
    
end


