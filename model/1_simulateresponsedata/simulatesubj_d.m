function [dnew, dold] = simulatesubj_d(modelname, theta, N, nX, nrep)
if nargin < 3; N = [150 150]; end
if nargin < 4; nX = 5; end
if nargin < 5; nrep = 30; end
% if nargin < 2; N = 150; end

% ======= INPUT VARIABLES =========
% MODELNAME: 
% THETA: 
% N: 
% NX: 
% NREP: 

Nnew = N(1); Nold = N(2);
% clear all
% N = 150;
% nX = 1;
% nrep = 100;
% theta = [11 0.8848    1.2506   -0.0599];

M = theta(1);
sigma = theta(2);
% nConf = 20;

sigs = 1;
J = (1/sigs^2 + 1/sigma^2);

% nX = 1;
% nsamp = 10;
nsamp_tnew = nrep*Nnew;
nsamp_told = nrep*Nold;
% Simulate fake data
% S = randn(N,M) * sigs;
% X = S + randn(N,M) * sigma;

dnew = nan(nX,Nnew*nrep); dold = nan(nX,Nold*nrep);
for iX = 1:nX;
    
    S = randn(Nold,M);
    Tnew = randn(nsamp_tnew,M);
    X = S + randn(Nold,M)*sigma;
    % X = randn(N,M) * sqrt(sigs^2 + sigma^2);
    
    allX = repmat(X,[1 1 nsamp_tnew]);
    % ========== NEW TRIALS ==========
    % t = randn(nsamp_t,M) * sigs;
    allTnew = permute(repmat(Tnew,[1 1 Nold]),[3 2 1]);
    
    temp = sum(bsxfun(@minus, allTnew, (1/sigs^2 + X/sigma^2)/J).^2,2);
    switch modelname
        case 'FP'
            dnewtemp = M/2 * log(1+sigs^2/sigma^2) + 0.5*sum(Tnew.^2,2) + log(squeeze(mean(exp(-.5.*temp*J))));
            dnew(iX,:) = dnewtemp';
        case 'FPheurs'
            dnew(iX,:) = -min((sum((allTnew-allX).^2,2)).^0.5); % 1 x nsamp_t
    end
    %     r = min(round(nConf/2.*(2./((1+exp(-(dnew-d0)/k))) - 1)+nConf/2+0.5),nConf);
    %     pnew = hist(r,1:nConf);
    %     pnew = pnew + 1;
    %     pnew = pnew/sum(pnew);
    
    % ========== OLD TRIALS ==========
    Xrep = repmat(X,[nrep 1]);
    allX = repmat(X,[1 1 nsamp_told]);
    Told = Xrep/sigma^2/J + randn(nsamp_told,M)/sqrt(J);
    allTold = permute(repmat(Told,[1 1 Nold]),[3 2 1]);
    temp = sum((allTold - (1/sigs^2 + allX/sigma^2)/J).^2,2);
    switch modelname
        case 'FP'
            doldtemp = M/2 * log(1+sigs^2/sigma^2) + 0.5*sum(Told.^2,2)+  log(squeeze(mean(exp(-.5.*temp*J))));
            dold(iX,:) = doldtemp';
        case 'FPheurs'
            dold(iX,:) = -min((sum((allTold-allX).^2,2)).^0.5); % 1 x nsamp_t
    end
    %     % confidence values etc
    %     r = min(round(nConf/2.*(2./((1+exp(-(dold-d0)/k))) - 1)+nConf/2 + 0.5),nConf);
    %     r = reshape(r,N,nrep);  % N x nrep
    %     pold = hist(r',1:nConf); % nConf x N matrix of histogram for each conf value
    %     pold = pold + 1; % correcting for potential zeros
    %     pold = pold/sum(pold); % normalizing
    
end


% plot(1:nConf,pnew,1:nConf,pold)