function [confidence, distancee, rho, pvalue, b, statss] = qualtrend_correlation(modelname, theta, binningfn, operationalization,nX,nS)
% simulate_data simulates fake data

if nargin < 6; nX = 1; end
if nargin < 7; nS = 1; end
% nConf = 20;
N = [150 150];
Nold = N(1); Nnew = N(2);

switch binningfn
    case 0
        slope = theta(end-2);
        nParams = 5;
    case 1
        k = theta(end-2);
        nParams = 5;
    case 2 % logarithmic
        a = theta(end-3);
        b = theta(end-2);
        nParams = 6;
    case 3 % power law
        a = theta(end-4);
        b = theta(end-3);
        gamma = theta(end-2);
        nParams = 7;
    case 4 % weibull
        scale = theta(end-5);
        shift = theta(end-4);
        a = theta(end-3);
        b = theta(end-2);
        nParams = 8;
end
d0 = theta(end-1);
sigma_mc = theta(end);
M = theta(1);
sigma = theta(2);

% check to make sure the length of theta is correct
assert(nParams == length(theta),'length of theta is not correct')

%     newHist = nan(nX,nConf);
[dOld, dNew, distancee] = deal(nan(Nold,nX,nS));
[d_new, distt] = deal(nan(Nnew,nS));
for iX = 1:nX;
    X = randn(Nold, M)*sqrt(sigma^2+1);
    
    for iS = 1:nS;
        SNew = randn(Nnew,M); % drawing new words
        SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
        
        [d_new(:,iS), dOld(:,iX,iS)] = calculate_d_FP(M, sigma, 1, Nnew, Nold, SNew, SOld, X);
        
        switch operationalization
            case 'min'  % distance bewteen new words and the closest old word
                distt(:,iS) = sqrt(squeeze(min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(SOld,[1,1,Nnew])).^2,2),[],1)));
            case 'mean'
                distt(:,iS) = sqrt(sum(bsxfun(@minus,SNew,mean(SOld)).^2,2));
            case 'median'
                distt(:,iS) = sqrt(sum(bsxfun(@minus,SNew,median(SOld)).^2,2));
        end
        dNew(:,iX,:) = d_new;
        distancee(:,iX,:) = distt;
    end
end
distancee = distancee(:);

% using absolute value for mapping
dNew_sign = sign(dNew(:)+d0); % -1 for respond new, +1 for respond old
q = abs(dNew(:)+d0);

% non-rounded confidence values
switch binningfn
    case 0 % linear
        confidence = slope.*q + d0;
    case 2 % logarithmic
        confidence = a.*log(q) + b;
    case 3 % power law
        confidence = a.*((q.^gamma - 1)./gamma) + b;
    case 4 % weibull
        confidence = a.*(1-exp(-(q./scale).^shift)) + b;
end

% add metacognitive noise
confidence = round(confidence + randn(Nold*nX*nS,1)*sigma_mc);
confidence(confidence < 1) = 1;
confidence(confidence > 20) = 20;
confidence(dNew_sign < 0,:) = 11 - confidence(dNew_sign < 0,:);
confidence(dNew_sign > 0,:) = 10 + confidence(dNew_sign > 0,:);

% statistics: correlation
[rho, pvalue] = corr([distancee confidence]); % correlation between each new word and operationalization of similarity
rho = rho(2); pvalue = pvalue(2);

% statistics: regression
[b,~,~,~,statss] = regress(confidence,[ones(length(confidence),1) distancee]);
% bestfitline = [1 min(distancee); 1 max(distancee)]*b;

