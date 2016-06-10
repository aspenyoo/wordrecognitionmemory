function [nnew_part, nold_part] = simulatedata_varySnew(modelname, theta, islogbinning, nX, nS, nConf, N)
% simulatedata_varySnew allows you to vary the distribution that new words
% (Snew) are drawn from 
%
% ======= INPUT VARIABLES ========
% MODELNAME: 'FP','FPheurs'
% THETA: a vector of length 6. [theta mu_Snew sigma_Snew]. if length 4,
% will procede where new words are drawn from the same distribution as old
% words.
% ISLOGBINNING: 1-logistic binning, 0-linear binning
% NX: number of X samples
% NS: number of S samples per X
% NCONF: number of confidence ratings (default: 20)
% N: vector of length 2 indicating number of new and old words. if scalar,
% will assume N = Nnew = Nold. (default = 150 = [150 150])

if nargin < 3; islogbinning = 1; end
if nargin < 4; nX = 1; end
if nargin < 5; nS = 1; end
if nargin < 6; nConf = 20; end
if nargin < 7; N = [150 150]; end

if isempty(N); N = [150 150]; end
if length(N) == 1; N = [N N]; end
if isempty(nConf); nConf = 20; end
Nnew = N(1); Nold = N(2);       % number of words

if strcmp(modelname,'uneqVar') % if UVSD model
    [pnew, pold] = responses_uneqVar(theta, islogbinning);
    nnew_part = round(Nnew*pnew);
    nold_part = round(Nold*pold);
    
else % if a optimal or heuristic model
    M = theta(1);                   % number of features
    sigma = theta(2);               % memory noise
    k = theta(3);                   % slope of logistic binning function
    d0 = theta(4);                  % shift of logistic binning function
    if length(theta) > 4; 
        mu_new = theta(5);          % mean of the new word distributino
    else
        mu_new = 0;
    end
    if length(theta) > 5;
        sigma_new = theta(6);       % mean of the old word distribution
    else
        sigma_new = 1;
    end
    
    L = nConf/2;
    
    sigs = 1;                       % width of word feature distribution
    J = 1/sigs^2 + 1/sigma^2;
    
    newHist = nan(nX,nConf);
    d_old = nan(Nold,nX,nS);
    d_new = nan(Nnew,nS);
    dNew = nan(Nold,nX,nS);
    for iX = 1:nX;
        X = randn(Nold, M)*sqrt(sigma^2+1);
        
        for iS = 1:nS;
            SNew = mu_new + randn(Nnew,M).*sigma_new; % drawing new words
            SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
            
            switch modelname
                case 'FP'
                    d_new(:,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
                    d_old(:,iX,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
                case 'FPheurs'
                    d_new(:,iS) = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1));
                    d_old(:,iX,iS) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2),[],1));
            end
        end
        dNew(:,iX,:) = d_new;
        
        % binning new words.
        if (islogbinning)
            newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
        else
            newHisttemp = min(max(round(m.*d_new(:) + b),1),nConf);
        end
        newHist(iX,:) = histc(newHisttemp,1:nConf);
    end    
    
    % binning old words
    if (islogbinning) % log binning
        oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
    else % lin binning
        oldHist = max(min(round(m.*d_old(:) + b),nConf),1);
    end
    oldHist = histc(oldHist,1:nConf); % histogram
    
    nnew_part = round(mean(newHist,1)/(nS));
    nold_part = round(oldHist'/(nX*nS));
end
