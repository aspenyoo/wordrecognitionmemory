function [nnew_part, nold_part] = simulate_data(modelname, theta, binningfn, nX, nS, nConf, N, plotstuff)
% simulatefake simulates fake data

if nargin < 3; binningfn = 1; end
if nargin < 4; nX = 1; end
if nargin < 5; nS = 1; end
if nargin < 6; nConf = 20; end
if nargin < 7; N = [150 150]; end
if nargin < 8; plotstuff = zeros(1,2); end

if isempty(N); N = [150 150]; end
if length(N) == 1; N = [N N]; end
if isempty(nConf); nConf = 20; end
Nnew = N(1); Nold = N(2);       % number of words

switch modelname
    
    case 'uneqVar' % UVSD model
        [pnew, pold] = responses_uneqVar(theta, binningfn);
        nnew_part = round(Nnew*pnew);
        nold_part = round(Nold*pold);
        
    case {'FP','FPheurs'} % optimal or heuristic model
        M = theta(1);                   % number of features
        sigma = theta(2);               % memory noise
        k = theta(3);                   % slope of logistic binning function
        d0 = theta(4);                  % shift of logistic binning function
        L = nConf/2;
        
        switch binningfn
            case {2,3,4};
                a = theta(3);
                b = theta(4);
                d0 = theta(5);
                sigma_mc = theta(6);
            case 5
                a = theta(3);
                b = theta(4);
                d0 = theta(5);
                lambda = theta(6);
                sigma_mc = theta(7);
        end
        
        sigs = 1;                       % width of word feature distribution
        J = 1/sigs^2 + 1/sigma^2;
        
        newHist = nan(nX,nConf);
        d_old = nan(Nold*nS,nX);
        d_new = nan(Nnew,nS);
        for iX = 1:nX;
            
            X = randn(Nold, M)*sqrt(sigma^2+1);
            Xrepp = repmat(X,[nS 1]);
            
            SNew = randn(Nnew*nS,M);
            SOld = (Xrepp/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
            
            switch modelname
                case 'FP'
                    [d_new, d_old(:,iX)] = calculate_d_FP(M, sigma, nS, Nnew, Nold, SNew, SOld, X);
%                     d_new = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
%                         sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2)))));
%                     d_old(:,iX) = M/2*log(1+1/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
%                         sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold*nS])).^2,2)))));
                case 'FPheurs'
                    d_new = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew*nS])).^2,2),[],1));
                    d_old(:,iX) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold*nS])).^2,2),[],1));
            end
            
            % binning new words.
            switch binningfn
                case 0
                    newHisttemp = min(max(round(m.*d_new(:) + b),1),nConf); % bounds: [1 20]
                    newHist(iX,:) = histc(newHisttemp,1:nConf);
                case 1
                    newHisttemp= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf); % bounds: [1 20]
                    newHist(iX,:) = histc(newHisttemp,1:nConf);
                case 2  % log mapping from LPR to confidence
                    error('not editted')                % remove after fix sigma_mc (see case 3)
                    d_new_sign = sign(d_new(:)+d0);                                                % -1 for respond new, +1 for respond old
                    newHisttemp = min(max(round(a.*log(abs(d_new(:)+d0))+b+randn.*sigma_mc),1),L) + L;        % confidence rating from 11 to 20
                    newHisttemp(d_new_sign < 0) = nConf+1 - newHisttemp(d_new_sign < 0);
                case 3 % log mapping on p(correct|evidence) instead of LPR
                    d_new_sign = sign(d_new(:)+d0);                                                % -1 for respond new, +1 for respond old
                    q = 1./(1+exp(-abs(d_new(:)+d0)));
                    conf = a.*log(q)+b;
                    newHisttemp = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                        (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
                    newHisttemp(d_new_sign < 0,:) = fliplr(newHisttemp(d_new_sign < 0,:));
                    newHist(iX,:) = (sum(newHisttemp)./nS)';                    % changing respond "new" words back to 1 to 10
                case 4 % log mapping on 1/(1-p(correct))
                    error('not editted')                % remove after fix sigma_mc (see case 3)
                    d_new_sign = sign(d_new(:)+d0);
                    q = 1./(1+exp(-abs(d_new(:)+d0)));
                    newHisttemp = min(max(round(a.*log(1-q)+b+randn.*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
                    newHisttemp(d_new_sign < 0) = nConf+1 - newHisttemp(d_new_sign < 0);    
                case 5 % generalized power law on p(correct)
                    d_new_sign = sign(d_new(:)+d0);                                                % -1 for respond new, +1 for respond old
                    q = 1./(1+exp(-abs(d_new(:)+d0)));
                    conf = a.*((q.^lambda - 1)./lambda)+b;
                    newHisttemp = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                        (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
                    newHisttemp(d_new_sign < 0,:) = fliplr(newHisttemp(d_new_sign < 0,:));
                    newHist(iX,:) = (sum(newHisttemp)./nS)';% changing respond "new" words back to 1 to 10
            end
        end
        
        % binning old words
        switch binningfn
            case 0 % linear binning
                oldHist = max(min(round(m.*d_old(:) + b),nConf),1);
                oldHist = histc(oldHist,1:nConf); % histogram
            case 1 % logistic binning
                oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
                oldHist = histc(oldHist,1:nConf); % histogram
            case 2
                error('not editted')                % remove after fix sigma_mc (see case 3)
                d_old_sign = sign(d_old(:)+d0);                                                % -1 for respond new, +1 for respond old
                oldHist = min(max(round(a.*log(abs(d_old(:)+d0))+b+randn.*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
                oldHist(d_old_sign < 0) = nConf+1 - oldHist(d_old_sign < 0);
            case 3 % log mapping on p(correct|evidence) instead of LPR
                d_old_sign = sign(d_old(:)+d0);                                                % -1 for respond new, +1 for respond old
                q = 1./(1+exp(-abs(d_old(:)+d0)));
                conf = a.*log(q)+b;
                oldHist = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                    (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
                oldHist(d_old_sign < 0,:) = fliplr(oldHist(d_old_sign < 0,:));
                oldHist = (sum(oldHist)./(nS*nX))';
            case 4 % log mapping on 1/(1-p(correct))
                error('not editted')                % remove after fix sigma_mc (see case 3)
                d_old_sign = sign(d_old(:)+d0);
                q = 1./(1+exp(-abs(d_old(:)+d0)));
                oldHist = min(max(round(a.*log(1-q)+b+randn.*sigma_mc),1),L)+L;        % confidence rating from 11 to 20
                oldHist(d_old_sign < 0) = nConf+1 - oldHist(d_old_sign < 0);                     % changing respond "new" words back to 1 to 10
            case 5 % power law mapping
                d_old_sign = sign(d_old(:)+d0);                                                % -1 for respond new, +1 for respond old
                q = 1./(1+exp(-abs(d_old(:)+d0)));
                conf = a.*((q.^lambda - 1)./lambda)+b;
                oldHist = 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
                    (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf-0.5)],conf)./(sigma_mc*sqrt(2))));
                oldHist(d_old_sign < 0,:) = fliplr(oldHist(d_old_sign < 0,:));
                oldHist = (sum(oldHist)./(nS*nX))';
        end
        
        switch binningfn
            case {0,1}
                nnew_part = round(mean(newHist,1)/(nS));
                nold_part = round(oldHist'/(nX*nS));
            case {2,3,4,5}
                nnew_part = mean(newHist);
                nold_part = oldHist';
        end
        
    case 'REM' % REM model
        M = theta(1);                   % number of features
        g = theta(2);                   % probability of success (for geometric distribution
        ustar = theta(3);               % probability of encoding something
        c = theta(4);                   % probability of encoding correct feature value
        m = theta(5);                   % number of storage attempts
        k = theta(6);                   % slope of logistic binning function
        d0 = theta(7);                  % shift of logistic binning function
        L = nConf/2;
        
        % figuring out X probabilities
        p0 = (1-ustar)^m;               % probability of x_ij= 0
        pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
        pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)
        
        % precalculated odds for mismatches with different X
        mismatchoddsVec = (c+(1-c).*(g.*(1-g).^(0:30)))./(g.*(1-g).^(0:30));
        
        newHist = nan(nX,nConf);
        d_old = nan(Nold*nS,nX);
        for iX = 1:nX;
            X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
            
            % generating new and old test words
            SNew = geornd(g,[1 M Nnew*nS])+1; % new words
            idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
            SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
            
            % decision variable values for new words
            idxmatch = bsxfun(@eq, SNew, X); % indices in which new words match X
            
            LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,(1-g).^-X) ; 
            LRmat(bsxfun(@eq,zeros(1,M,Nold*nS),X)) = 1;
            d_new = log(mean(prod(LRmat,2),1));   % log odds
            
            % decision variable values for old words
            idxmatch = bsxfun(@eq, permute(SOld, [3 2 1]), X); % indices in which new words match X
            
            LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,exp(-X.*log(1-g)));
            LRmat(bsxfun(@eq,zeros(1,M,Nold*nS),X)) = 1;
            d_old(:,iX) = log(mean(prod(LRmat,2),1));
            

            % binning new words.
            if (binningfn)
                newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
            else
                newHisttemp = min(max(round(slope.*d_new(:) + b),1),nConf);
            end
            newHist(iX,:) = histc(newHisttemp,1:nConf);
        end
        
        
        % binning old words
        if (binningfn) % log binning
            oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
        else % lin binning
            oldHist = max(min(round(m.*d_old(:) + b),nConf),1);
        end
        oldHist = histc(oldHist,1:nConf); % histogram
        
        nnew_part = round(mean(newHist,1)/(nS));
        nold_part = round(oldHist'/(nX*nS));
end




% % ========== PLOTTING STUFF ==========

if plotstuff(2) == 1;
    if plotstuff(1) ==1;
        subplot(2,1,1)
    end
    switch modelname
        case 'uneqVar'
            responses_uneqVar( theta, binningfn, 0, 1);
        case {'FP','FPheurs'}
            binningparameters = theta(3:end);
            memdistplot(d_old(:), d_new(:), binningparameters, binningfn);
        case 'REM'
            binningparameters = theta(3:end);
            pcorr_new = 1./(1+exp(-abs(d_new(:))));
            pcorr_old = 1./(1+exp(-abs(d_old(:))));
            memdistplot(pcorr_old(:), pcorr_new(:), binningparameters, binningfn);
    end
end

if plotstuff(1) == 1;
    if plotstuff(2) == 1;
        subplot(2,1,2)
    else
        figure;
    end
    confdistplot(nnew_part/N(1), nold_part/N(2))
end