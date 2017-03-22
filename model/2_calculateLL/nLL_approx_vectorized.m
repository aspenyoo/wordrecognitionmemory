function [ varargout ] = nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part, logflag, fixparams, nX, nS, nConf )
% nLL_approx calculates the negative log likelihood using an approximation
% method
%
% [nLL] = nLL_approx_vectorized(MODELNAME, THETA, ISLOGBINNING, NNEW_PART,
% NOLD_PART) calculates the negative log p(nnew_part,nold_part|theta,model)
%
% ===== INPUT VARIABLES =====
% MODELNAME: 'FP','FPheurs','VP','VPheurs','uneqVar', 'REM'
% THETA: parameter values
% BINNINGFN: mapping from MEMSTRENGTHVAR to confidence reports. 0: linear,
% 1: logistic, 2: log, 3: power law mapping, 4: weibull
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
% Aspen Yoo - Nov 30, 2016
if nargin < 6; 
    switch modelname
        case 'FP'
            logflag = [1 1];
        case 'REM'
            logflag = [1 0 0 0 0];
        case 'UVSD'
            logflag = [0 1];
    end
    switch binningfn
        case 3 % power law
            logflag = [logflag 0 0 0];
        case 4 % cumulative weibull
            logflag = [logflag 0 0 0 0];
    end
    logflag = [logflag 0 1];
    logflag = logical(logflag);
end
if nargin < 7; fixparams = []; end
if nargin < 8; nX = 30; end
if nargin < 9; nS = 100; end
if nargin < 10; nConf = 20; end

if ~isempty(fixparams) && (nargin < 6); logflag(fixparams(1,:)) = [];end
theta(logflag) = exp(theta(logflag)); % exponentiating the appropriate free paraemters

% if theta has fixed parameters, adjust accordingly
if ~isempty(fixparams)
    nParams = length(theta) + size(fixparams,2);
    tempTheta = nan(1,nParams);
    unfixedparams = 1:nParams;
    unfixedparams(fixparams(1,:)) = [];
    tempTheta(fixparams(1,:)) = fixparams(2,:);
    tempTheta(unfixedparams) = theta;
    theta = tempTheta;
end

switch modelname
    case {'FP','FPheurs','UVSD'}
        nParams = 2;
    case 'REM'
        nParams = 5;
end
switch binningfn
    case 0
        slope = theta(end-2);
        nParams = nParams + 3;
    case 1
        k = theta(end-2);
        nParams = nParams + 3;
    case 2 % logarithmic
        a = theta(end-3);
        b = theta(end-2);
        nParams = nParams + 4;
    case 3 % power law
        a = theta(end-4);
        b = theta(end-3);
        gamma = theta(end-2);
        nParams = nParams + 5;
    case 4 % weibull
        scale = theta(end-5);
        shape = theta(end-4);
        a = theta(end-3);
        b = theta(end-2);
        nParams = nParams + 6;
end
d0 = theta(end-1);
sigma_mc = theta(end);

if strcmp(modelname,'UVSD') % if uneqVar
    [pnew, pold, confbounds] = responses_uneqVar(theta, binningfn);
    nLL = -sum(log(pnew).*nnew_part) - sum(log(pold).*nold_part);
else % if FP, FPheurs, or REM
    
    % stuff that doesn't depend on model, mapping, or memstrengthvar
    Nold = sum(nold_part); Nnew = sum(nnew_part);
    lapse = 0.01;             % lapse rate
    
    % parameter names
    switch modelname
        case {'FP','FPheurs'}
            M = round(theta(1));
            sigma = theta(2);
        case 'REM'
            M = round(theta(1));                   % number of features
            g = theta(2);                   % probability of success (for geometric distribution)
            ustar = theta(3);               % probability of encoding something
            c = theta(4);                   % probability of encoding correct feature value
            m = theta(5);                   % number of storage attempts
            
            p0 = (1-ustar)^m;               % probability of x_ij= 0
            pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
            pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)
    end
    
    % old parameterization of REM and FP fits
    %     sigma_mc = theta(end-5);
    %     d0 = theta(end-4);
    %     a = theta(end-3);
    %     b = theta(end-2);
    %     scale = theta(end-1);
    %     shift = theta(end);
    %     nParams = nParams +6;
    
    % check to make sure the length of theta is correct
    assert(nParams == length(theta),'length of theta is not correct')
    
    % calculate LL
    LL_new = nan(1,nX);
    d_old = nan(Nold*nS,nX);
    d_newtotal = nan(Nnew*nS,nX);
    newHisttotal = nan(nX,nConf);
    for iX = 1:nX
        
        % calculate d
        switch modelname
            case 'FP'
                X = randn(Nold, M)*sqrt(sigma^2+1);
                SNew = randn(Nnew*nS,M);
                SOld = (repmat(X,[nS 1])/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
                
                [d_new, d_old(:,iX)] = calculate_d_FP(M, sigma, nS, Nnew, Nold, SNew, SOld, X);
                d_newtotal(:,iX) = d_new;
            case 'FPheurs'
                X = randn(Nold, M)*sqrt(sigma^2+1);
                SNew = randn(Nnew*nS,M);
                SOld = (repmat(X,[nS 1])/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
                
                d_new = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew*nS])).^2,2),[],1));
                d_old(:,iX) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold*nS])).^2,2),[],1));
                d_newtotal(:,iX) = d_new;
            case 'REM'
                X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
                SNew = geornd(g,[Nnew*nS M])+1; % new words
                idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
                SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
                
                maxidx = max(X(:));
                matchoddsVec = (c+(1-c).*g.*(1-g).^(0:maxidx))./(g.*(1-g).^(0:maxidx));
                [d_new, d_old(:,iX)] = calculate_d_REM(M, g, c, nS, Nnew, Nold, SNew, SOld, X, matchoddsVec);
                d_newtotal(:,iX) = d_new;
        end
        
        % using absolute value for mapping
        d_new_sign = sign(d_new(:)+d0); % -1 for respond new, +1 for respond old
        q = abs(d_new(:)+d0);
        
        % non-rounded confidence values
        switch binningfn
            case 0 % linear
                conf = slope.*q + d0;
            case 1 % logistic
            case 2 % logarithmic
                conf = a.*log(q) + b;
            case 3 % power law
                conf = a.*((q.^gamma - 1)./gamma) + b;
            case 4 % weibull
                conf = a.*(1-exp(-(q./scale).^shape)) + b;
        end
        
        % binning with or without metacognitive noise
        if (sigma_mc)
            newHist = [zeros(length(conf),nConf/2) 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf/2-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf/2-0.5)],conf)./(sigma_mc*sqrt(2))))];
            newHist(d_new_sign < 0,:) = fliplr(newHist(d_new_sign < 0,:));
            newHist = (sum(newHist)./(nS*nX))';
        else
            conf = min(max(round(conf),1),nConf/2)+nConf/2;
            conf(d_new_sign<0) = nConf+1 - conf(d_new_sign<0);
            newHist = histc(conf,1:nConf); % histogram
        end
        newHisttotal(iX,:) = newHist'; % in case you want nnew_part
        
        pnew = lapse/nConf + (1-lapse)*(newHist/sum(newHist));
        LL_new(iX) = nnew_part*log(pnew);
    end
    
    % BINNING OLD WORDS. p(conf|X0,C)
    
    % using absolute value for mapping
    d_old_sign = sign(d_old(:)+d0); % -1 for respond new, +1 for respond old
    q = abs(d_old(:)+d0);
    
    % non-rounded confidence values
    switch binningfn
        case 0 % linear
            conf = slope.*q + d0;
        case 1 % logistic
        case 2 % logarithmic
            conf = a.*log(q) + b;
        case 3 % power law
            conf = a.*((q.^gamma - 1)./gamma) + b;
        case 4 % weibull
            conf = a.*(1-exp(-(q./scale).^shape)) + b;
    end
    
    % histograms of confidence
    if (sigma_mc) % if there is metacognitive noise
        oldHist = [zeros(length(conf),nConf/2) 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf/2-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
            (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf/2-0.5)],conf)./(sigma_mc*sqrt(2))));];
        oldHist(d_old_sign < 0,:) = fliplr(oldHist(d_old_sign < 0,:));
        oldHist = (sum(oldHist)./(nS*nX))';
    else
        conf = min(max(round(conf),1),nConf/2)+nConf/2;
        conf(d_old_sign<0) = nConf+1 - conf(d_old_sign<0);
        oldHist = histc(conf,1:nConf); % histogram
    end
    
    pold = lapse/nConf + (1-lapse)*(oldHist/sum(oldHist)); % normalizing
    
    % calculating nLL
    LL_old = nold_part*log(pold);
    LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
    nLL = -LL_new-LL_old;
    
end

switch nargout
    case 1
        varargout = {nLL};
    case 2
        if ~strcmp(modelname,'UVSD')
            pnew = sum(newHisttotal)/sum(newHisttotal(:));
        end
        varargout = {pnew, pold};
    case 5
        nBins = 50;
        
        if strcmp(modelname,'UVSD')
            nSamples = 10000;
            
            mu_old = theta(1);
            sigma_old = theta(2);
            
            old_samples = mu_old + randn(1,nSamples).*sigma_old;
            new_samples = randn(1,nSamples);
            d_newtotal = -log(sigma_old) - 1/2.*( ((new_samples-mu_old).^2)./sigma_old.^2 - new_samples.^2);
            d_old = -log(sigma_old) - 1/2.*( ((old_samples-mu_old).^2)./sigma_old.^2 - old_samples.^2);
            
        end
        %
        %             x = linspace(min([0 mu_old]-nSDs.*[1 sigma_old]),max([0 mu_old]+nSDs.*[1 sigma_old]),nSamples);
        %             pold_x = normpdf(x,mu_old,sigma_old);
        %             pnew_x = normpdf(x,0,1);
        %             d = -log(sigma_old) - 1/2.*( ((x-mu_old).^2)./sigma_old.^2 - x.^2);
        %
        %             counts_new = pnew_x.*d;
        %             centers_new = linspace(-3,3,50);
        %             centers_old = linspace(mu_old-3*sigma_old,mu_old+3*sigma_old,50);
        %             counts_new = normpdf(centers_new);
        %             counts_old = normpdf(centers_old,theta(1),theta(2));
        
        
        %             xx = linspace(0,10,100);
        %             figure;
        %             plot(xx,a.*(1-exp(-(xx./gamma).^k)) + b,'k-')
        %             hold on;
        %             plot(confbounds+d0,binvalues,'or')
        
        [counts_new,centers_new] = hist(d_newtotal(:),nBins);
        [counts_old,centers_old] = hist(d_old(:),nBins);
        
        
        switch binningfn
            case 0 % linear
                binvalues = 1.5:(nConf-0.5);
                confbounds = (binvalues-nConf/2-0.5)./shape - d0;
            case 1 % logistic
                binvalues = 1.5:(nConf-0.5);
                confbounds = -shape.*log((nConf+0.5)./(binvalues-0.5)-1);
            case 2 % logarithmic
                binvalues = 1.5:(nConf/2 -0.5);
                confbounds = exp((binvalues-b)./a);
            case 3 % power law
                binvalues = 1.5:(nConf/2 -0.5);
                tempp = gamma.*(binvalues -b)./a + 1;
                tempp(tempp < 0) = nan;
                confbounds = tempp.^(1/gamma);
            case 4 % weibull
                binvalues = 1.5:(nConf/2 -0.5);
                tempp = 1 - (binvalues -b)./a;
                tempp(tempp < 0) = nan;
                tempp(tempp > 1) = nan;
                confbounds = scale.*(-log(tempp)).^(1/shape);
        end
        
        varargout = {centers_new, counts_new, centers_old, counts_old, confbounds};
end