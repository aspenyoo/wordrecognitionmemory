function [d_mean, lengths, nnew_part, nold_part] = simulatedata_snewonline(modelname, theta, islogbinning, nX, nS, nConf, N)
% simulatefake simulates fake data using FP model
%
% ======== INPUT VARIABLES =======
% MODELNAME: FP, FPheurs, uneqVar
% THETA: model parameters
% ISLOGBINNING: 1: logistic binning. 0: linear binning
% NX: number of X samples
% NS: number of S samples
% NCONF: number of confidence values (default 20)
% N: vector of length 2 the the number of new and old words, respectively.

dissimilarityMax = 3; % max distance from center
if nargin < 3; islogbinning = 1; end
if nargin < 4; nX = 1; end
if nargin < 5; nS = 1; end
if nargin < 6; nConf = 20; end
if isempty(nConf); nConf = 20; end
if nargin < 7; N = [150 150]; end
if isempty(N); N = [150 150]; end
if length(N) == 1; N = [N N]; end

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
    L = nConf/2;
    
    sigs = 1;                       % width of word feature distribution
    J = 1/sigs^2 + 1/sigma^2;
    
    newHist = nan(nX,nConf);
    dOld = nan(Nold,nX,nS);
    dNew = nan(Nnew,nX,nS);
    d_new = nan(Nnew,nS);
    
    lengths = linspace(0,dissimilarityMax,Nnew);
    SNew = repmat(lengths',[1 M])./sqrt(M); % drawing new words that increase in dissimilarity (that are equally long in each feature dimension)
    for iX = 1:nX;
        X = randn(Nold, M)*sqrt(sigma^2+1);
        
        for iS = 1:nS;
            %             SNew = randn(Nnew,M); % drawing new words
            SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
            
            switch modelname
                case 'FP'
                    d_new(:,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
                    dOld(:,iX,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
                case 'FPheurs'
                    d_new(:,iS) = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1));
                    dOld(:,iX,iS) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2),[],1));
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
        oldHist = min(round(L.*(2./((1+exp(-(dOld(:)-d0)./k))) - 1)+10.5),nConf);
    else % lin binning
        oldHist = max(min(round(m.*dOld(:) + b),nConf),1);
    end
    oldHist = histc(oldHist,1:nConf); % histogram
    
    nnew_part = round(mean(newHist,1)/(nS));
    nold_part = round(oldHist'/(nX*nS));
    
    % average d as a with increasing similarity
    d_mean = squeeze(mean(mean(dNew,3),2));
end

plot(lengths,d_mean)
defaultplot
