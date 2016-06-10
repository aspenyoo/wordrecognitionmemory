function [rho_new, pvalue_new] = correlation_LPR_distanceNewX(modelname, M, sigma)
% calculates pearson r correlation coefficient between 1) the distance of 
% new words and the closest memory and 2) the LPR
% for a given MODELNAME, M, SIGMA

% modelname = 'FP';
% M = 20;
% sigma = 1;

nX = 1;
nS = 1;
nConf = 20;
N = [150 150];
Nold = N(1); Nnew = N(2);

sigma_s = 1;                       % width of word feature distribution
J = 1/sigma_s^2 + 1/sigma^2;

%     newHist = nan(nX,nConf);
dOld = nan(Nold,nX,nS);
d_new = nan(Nnew,nS);
dNew = nan(Nold,nX,nS);
dist_new = nan(Nnew,nS);
distNew = nan(Nold,nX,nS);
distOld = nan(Nold,nS,nS);
for iX = 1:nX;
    X = randn(Nold, M)*sqrt(sigma^2+1);
    
    for iS = 1:nS;
        SNew = randn(Nnew,M); % drawing new words
        SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
        
        switch modelname
            case 'FP'
                d_new(:,iS) = M/2*log(1+sigma_s^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                    sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
                dOld(:,iX,iS) = M/2*log(1+sigma_s^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                    sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
                dist_new(:,iS) = sqrt(squeeze(min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1))); % distance between new words and the closest noisy memory
                distOld(:,iX,iS) = sqrt(squeeze(min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1))); % distance between old words and the closest noisy memory
                %                     dist_new(:,iS) = sqrt(squeeze(min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(SOld,[1,1,Nnew])).^2,2),[],1))); % distance from the closest old word
                %                     distOld(:,iX,iS) = sqrt(squeeze(min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(SOld,[1,1,Nnew])).^2,2),[],1))); % distance from the closest old word
            case 'FPheurs'
                d_new(:,iS) = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1));
                dOld(:,iX,iS) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2),[],1));
        end
    end
    dNew(:,iX,:) = d_new;
    distNew(:,iX,:) = dist_new;
end

% statistics: correlation
[rho_new, pvalue_new] = corr([distNew dNew],'type','Spearman'); % new word: distance by LPR
% using spearman rho because it doesn't assume gaussian distributed data
% and works for ordinal data. 
rho_new = rho_new(2);
pvalue_new = pvalue_new(2);

