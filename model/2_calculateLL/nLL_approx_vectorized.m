function [ varargout ] = nLL_approx_vectorized( modelname, theta, binningfn, memstrengthvar, nnew_part, nold_part, fixparams, nX, nS, nConf )
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
% 1: logistic, 2: log, 3: power law mapping
% MEMSTRENGTHVAR: variable used for "memory strength). 0: LPR, 1: p(correct), 2: 1/(p(incorrect)))
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
% Aspen Yoo - Aug 19, 2016

if nargin < 7; fixparams = []; end
if nargin < 8; nX = 30; end
if nargin < 9; nS = 50; end
if nargin < 10; nConf = 20; end

rng('shuffle')

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

if strcmp(modelname,'UVSD') % if uneqVar
    [pnew, pold, confbounds] = responses_uneqVar(theta, binningfn);
    nLL = -sum(log(pnew).*nnew_part) - sum(log(pold).*nold_part);
else % if FP, FPheurs, or REM
    
    % stuff that doesn't depend on model, mapping, or memstrengthvar
    Nold = sum(nold_part); Nnew = sum(nnew_part);
    L = nConf/2;
    lapse = 0.01;             % lapse rate

    % parameter names
    switch modelname
        case {'FP','FPheurs'};
            M = theta(1);
            sigma = theta(2);
            nParams = 2;
        case 'REM';
            M = theta(1);                   % number of features
            g = theta(2);                   % probability of success (for geometric distribution
            ustar = theta(3);               % probability of encoding something
            c = theta(4);                   % probability of encoding correct feature value
            m = theta(5);                   % number of storage attempts
            nParams = 5;
            
            p0 = (1-ustar)^m;               % probability of x_ij= 0
            pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
            pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)
    end
    
    switch binningfn
        case {0,1}              % linear (0) and logistic (1) mapping
            k = theta(end-2);
            d0 = theta(end-1);
            sigma_mc = theta(end);
            nParams = nParams + 3;
        case 2        % logarithmic mapping
            a = theta(end-3);
            b = theta(end-2);
            d0 = theta(end-1);
            sigma_mc = theta(end);
            nParams = nParams + 4;
        case 3              % power law mapping from p(correct)
            a = theta(end-4);
            b = theta(end-3);
            d0 = theta(end-2);
            lambda = theta(end-1);
            sigma_mc = theta(end);
            nParams = nParams + 5;
    end
    
    % check to make sure the length of theta is correct
    assert(nParams == length(theta),'length of theta is not correct')
    
    % calculate LL
    LL_new = nan(1,nX);
    d_old = nan(Nold*nS,nX);
    newHisttotal = nan(nX,nConf);
    for iX = 1:nX;
        
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
        
        %             % plot so see how nS affects data
        %             LL_newtemp = nan(1,nS);
        %             for iS = 1:nS;
        %                 % binning new words.
        %                 switch binningfn
        %                     case 0 % linear
        %                         newHist = min(max(round(m.*d_new(1:iS*Nnew) + b),1),nConf);                            % bounds: [1 20]
        %                     case 1 % logistic
        %                         newHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(1:iS*Nnew)-d0)./k)) - 1)),nConf);   % bounds: [1 20]
        %                     case 2 % log
        %                         d_new_sign = sign(d_new(1:iS*Nnew)+d0);                                                % -1 for respond new, +1 for respond old
        %                         newHist = min(max(round(a.*log(abs(d_new(1:iS*Nnew)+d0))+b+randn(length(q),1).*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
        %                         newHist(d_new_sign < 0) = nConf+1 - newHist(d_new_sign < 0);                     % changing respond "new" words back to 1 to 10
        %                     case 3 % log mapping on p(correct|evidence) instead of LPR
        %                         d_new_sign = sign(d_new(1:iS*Nnew)+d0);                                                % -1 for respond new, +1 for respond old
        %                         q = 1./(1+exp(-abs(d_new(1:iS*Nnew)+d0)));
        %                         newHist = min(max(round(a.*log(q)+b+randn(length(q),1).*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
        %                         newHist(d_new_sign < 0) = nConf+1 - newHist(d_new_sign < 0);                     % changing respond "new" words back to 1 to 10
        %                     case 4 % log mapping on 1/(1-p(correct))
        %                         d_new_sign = sign(d_new(1:iS*Nnew)+d0);
        %                         q = 1./(1+exp(-abs(d_new(1:iS*Nnew)+d0)));
        %                         newHist = min(max(round(a.*log(1-q)+b+randn(length(q),1).*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
        %                         newHist(d_new_sign < 0) = nConf+1 - newHist(d_new_sign < 0);                     % changing respond "new" words back to 1 to 10
        %                 end
        %                 newHist = histc(newHist,1:nConf);
        %                 pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
        %                 LL_newtemp(iS) = nnew_part*log(pnew);
        %             end
        %             figure;
        %             plot(1:nS,LL_newtemp,'ok')
        %             defaultplot
        %             pause
        
        % using absolute value for logarithmic and power law mapping
        if any(binningfn == [2 3])
            d_new_sign = sign(d_new(:)+d0); % -1 for respond new, +1 for respond old
            q = abs(d_new(:)+d0);
        else
            q = d_new(:) + d0;
        end
        
        % adjust variable depending on what memstrengthvar is
        switch memstrengthvar
            case 1      % p(correct|evidence). [0.5 1]
                q = 1./(1+exp(-q));
            case 2      % 1/(1-p(correct|evidence)). [2, Inf)
                q = 1./(1-1./(1+exp(-q)));
        end
        
        % calculate LLnew
        switch binningfn
            case 0 % linear
                conf = k.*q;                            % bounds: [1 20]
            case 1 % logistic
                conf= 0.5+ 2*L./(1+exp(-(q)./k));   % bounds: [1 20]
            case 2                      % logarithmic mapping
                conf = a.*log(q)+b;
            case 3 % generalized power law mapping
                conf = a.*((q.^lambda - 1)./lambda)+b;
        end
        
        if (binningfn == 0) && (memstrengthvar ==0); % confidence is 10.5 at decision boundary (0)
            conf = conf + nConf/2+0.5;
        end
        
        % binning with or without metacognitive noise
        if (sigma_mc)
            switch binningfn
                case {0,1}
                    newHist = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                    (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
                case {2,3}
                    newHist = [zeros(length(conf),nConf/2) 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf/2-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                        (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf/2-0.5)],conf)./(sigma_mc*sqrt(2))))];
                    newHist(d_new_sign < 0,:) = fliplr(newHist(d_new_sign < 0,:));
            end
            newHist = (sum(newHist)./(nS*nX))';
        else
            switch binningfn
                case {0,1}
                    conf = min(max(round(conf),1),nConf);
                case {2,3}
                    conf = min(max(round(conf),1),nConf/2)+nConf/2;
                    conf(d_new_sign<0) = nConf+1 - conf(d_new_sign<0);
            end
            newHist = histc(conf,1:nConf); % histogram
        end
        newHisttotal(iX,:) = newHist'; % in case you want nnew_part
        
        pnew = lapse/nConf + (1-lapse)*(newHist/sum(newHist));
        LL_new(iX) = nnew_part*log(pnew);
    end
    
    % BINNING OLD WORDS. p(conf|X0,C)
    
    % using absolute value for logarithmic and power law mapping
    if any(binningfn == [2 3])
        d_old_sign = sign(d_old(:)+d0); % -1 for respond new, +1 for respond old
        q = abs(d_old(:)+d0);
    else
        q = d_old(:) + d0;
    end
    
    % adjust variable depending on operationalization of "memory strength"
    switch memstrengthvar
        case 1      % p(correct|evidence). [0.5 1]
            q = 1./(1+exp(-q));
        case 2      % 1/(1-p(correct|evidence)). [2, Inf)
            q = 1./(1-1./(1+exp(-q)));
    end
    
    % non-rounded confidence values
    switch binningfn
        case 0 % linear
            conf = k.*q;                            % bounds: [1 20]
        case 1 % logistic
            conf= 0.5+2*L./(1+exp(-(q)./k));   % bounds: [1 20]
        case 2 % log mapping on p(correct|evidence) instead of LPR
            conf = a.*log(q)+b;
        case 3 % power law mapping
            conf = a.*((q.^lambda - 1)./lambda)+b;
    end
    
    % histograms of confidence
    if (sigma_mc) % if there is metacognitive noise
        switch binningfn
            case {0,1}
                oldHist = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                    (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
            case {2,3}
                oldHist = [zeros(length(conf),nConf/2) 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf/2-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                    (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf/2-0.5)],conf)./(sigma_mc*sqrt(2))));];
                oldHist(d_old_sign < 0,:) = fliplr(oldHist(d_old_sign < 0,:));
        end
        oldHist = (sum(oldHist)./(nS*nX))';
    else
        switch binningfn
            case {0,1}
                conf = min(max(round(conf),1),nConf);
                
            case {2,3}
                conf = min(max(round(conf),1),nConf/2)+nConf/2;
                conf(d_old_sign<0) = nConf+1 - conf(d_old_sign<0);
        end
        oldHist = histc(conf,1:nConf); % histogram
    end
    
    %         % plotting for X samples
    %         for iX = 1:nX;
    %             % BINNING OLD WORDS. p(conf|X0,C)
    %             d_old_sign = sign(d_old(1:iX*nS*Nold)+d0);                                                % -1 for respond new, +1 for respond old
    %             q = 1./(1+exp(-abs(d_old(1:iX*nS*Nold)+d0)))';
    %             q(d_old_sign < 0) = nConf+1 - q(d_old_sign <0);
    %             oldHist = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],q)./(sigma_mc*sqrt(2))) - ...
    %                     (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],q)./(sigma_mc*sqrt(2))));
    %             oldHist = round(sum(oldHist)./nS)';
    %             pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
    %
    %             % calculating nLL
    %             LL_oldtemp(iX) = nold_part*log(pold);
    %             LL_newtemp(iX) = max(LL_new(1:iX)) + log(mean(exp(LL_new(1:iX)-max(LL_new(1:iX))))); % average over X
    %             maxnew(iX) = max(LL_new(1:iX));
    % %             blah(iX) = max(LL_new(iX))-LL_new(1:iX);
    %         end
    %         nLLtemp = -LL_newtemp-LL_oldtemp;
    % %         figure; plot(1:nX,nLLtemp,'ko');defaultplot
    % %         hold on;
    % %         plot(1:nX,-LL_oldtemp,'bo');
    %
    %         figure; hold on;
    %         plot(1:nX,LL_newtemp,'o','Color',aspencolors('gold'));
    %         plot(1:nX,LL_new,'*','Color',aspencolors('dustygold'));
    %         plot(1:nX,maxnew,'.')
    %         plot(1:nX,LL_oldtemp,'o','Color',aspencolors('babyblue'));
    %         defaultplot
    %         xlabel('nX')
    %         ylabel('LL for new words')
    %         legend('new LL','LL for each sample','max of samples','old LL')
    
    pold = lapse/nConf + (1-lapse)*(oldHist/sum(oldHist)); % normalizing
    
    % calculating nLL
    LL_old = nold_part*log(pold);
    LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
    nLL = -LL_new-LL_old;
    
end

switch nargout
    case 1;
        varargout = {nLL};
    case 2;
        if ~strcmp(modelname,'UVSD')
            pnew = sum(newHisttotal)/sum(newHisttotal(:));
        end
        varargout = {pnew, pold};
    case 5;
        if strcmp(modelname,'UVSD')
            centers_new = linspace(-3,3,50);
            centers_old = linspace(theta(1)-3*theta(2),theta(1)+3*theta(2),50);
            counts_new = normpdf(centers_new);
            counts_old = normpdf(centers_old,theta(1),theta(2));
        else
            [counts_old,centers_old] = hist(d_old(:),50);
            [counts_new,centers_new] = hist(d_newtotal(:),50);
            
            switch binningfn
                case 0 % linear
                    binvalues = 1.5:(nConf-0.5);
                    switch memstrengthvar
                        case 0 % LPR
                            confbounds = binvalues./k - d0;
                        case 1 % p(corr)
                            confbounds = -log(1./(binvalues./k - d0 + 0.5)-1);
                        case 2 % 1/p(incorr)
                    end
                case 1 % logistic
                    binvalues = 1.5:(nConf-0.5);
                    switch memstrengthvar
                        case 0 % LPR
                            confbounds = -k.*log((nConf+0.5)./(binvalues-0.5)-1);
                        case 1 % p(corr)
                            confbounds = -log(1./(-k.*log((nConf+0.5)./(binvalues-0.5)-1)+0.5)-1);
                        case 2 % 1/p(incorr)
                    end
                case 2 % logarithmic
                    switch memstrengthvar
                        case 0 % LPR
                        case 1 % p(corr)
                        case 2 % 1/p(incorr)
                    end
                case 3 % power law
                    switch memstrengthvar
                        case 0 % LPR
                        case 1 % p(corr)
                        case 2 % 1/p(incorr)
                    end
            end

        end
        varargout = {centers_new, counts_new, centers_old, counts_old, confbounds};
end