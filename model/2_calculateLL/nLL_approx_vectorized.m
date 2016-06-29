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


switch modelname
    case 'uneqVar'
        [pnew, pold] = responses_uneqVar(theta, islogbinning);
        nLL = -sum(log(pnew).*nnew_part) - sum(log(pold).*nold_part);
    case {'FP','FPheurs'}
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
        
        LL_new = nan(1,nX);
        d_old = nan(Nold*nS,nX);
        for iX = 1:nX;
            
            X = randn(Nold, M)*sqrt(sigma^2+1);
            Xrepp = repmat(X,[nS 1]);
            
            SNew = randn(Nnew*nS,M);
            SOld = (Xrepp/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
            
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
        
    case 'REM'
        
        M = theta(1);                   % number of features
        g = theta(2);                   % probability of success (for geometric distribution
        ustar = theta(3);               % probability of encoding something
        c = theta(4);                   % probability of encoding correct feature value
        m = theta(5);                   % number of storage attempts
        L = nConf/2;
        if (islogbinning); 
            k = theta(6);  d0 = theta(7);
        else c1 = theta(6); c2 = theta(7);
        end
        Nold = sum(nold_part); Nnew = sum(nnew_part);
        lambda = 0.01;             % lapse rate
        
        if M ~= floor(M)
            M = floor(M);
        end
        
        if ~(islogbinning)
            slope =(nConf-2)/(c2-c1);       % slope
            b = 1.5-slope*c1;               % y-intercept
        end
        
        
        % figuring out X probabilities
        p0 = (1-ustar)^m;               % probability of x_ij= 0
        pM = (1-p0)*geocdf(m-1,c);      % probability of x_ij = s_ij (match)
        pQ = 1 - p0 - pM;               % probability drawn randomly (mismatch)
        
        % precalculated odds for mismatches with different X
        mismatchoddsVec = (c+(1-c).*(g.*(1-g).^(0:30)))./(g.*(1-g).^(0:30));

        d_old = nan(Nold*nS,nX);
        LL_new = nan(1,nX);
        
% %         % ====== NO FOR LOOP ======
% %         X = permute(binornd(1,1-p0,[Nold M nX]).*(geornd(g,[Nold M nX])+1),[1 2 4 3]);
% % %         Xrepp = repmat(X,[nS 1]);
% %         
% %         % generating new and old test words
% %         SNew = geornd(g,[1 M Nnew*nS nX])+1; % new words
% %         idx = logical(binornd(1,pQ,[Nold*nS M 1 nX]) + repmat((X == 0),[nS 1 1 1])); % indices of randomly drawn features
% %         SOld = permute((1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M 1 nX]) + 1),[3 2 1 4]);
% %         
% %         % new word match and mismatch
% %         idxmatch = bsxfun(@eq, SNew, X);
% %         mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0), ~idxmatch),2);
% %          
% %         tic;
% %         Xrep = repmat(X,[1 1 Nold*nS 1]);        % making X the proper size
% %         idxmatch = bsxfun(@eq, SNew, X);         % logicals of which elements match
% %         matmat = ones(Nold, M, Nold*nS, nX);     % setting size of matrix
% %         matmat(idxmatch) = mismatchoddsVec(Xrep(idxmatch)); % filling in matching elements with appropriate value
% %         toc
% %         
% %         % matmat without linear indexing
% %         tic;
% %         Xrep = repmat(X,[1 1 Nold*nS 1]);        % making X the proper size
% %         idxmatch = (Xrep.*bsxfun(@eq, SNew, X)) + 1;         % logicals of which elements match
% %         mismatchoddsVec = [1 mismatchoddsVec];
% %         matmat = mismatchoddsVec(idxmatch);
% %         toc
% %         
% %         d_new = squeeze(-log(Nnew) + log(sum((1-c).^mismatchSum.*prod(matmat,2),1)));   % log odds
% %         
% %         % binning new words (separately for each iX)
% %         if (islogbinning)
% %             newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new-d0)./k)) - 1)),nConf);
% %         else
% %             newHisttemp = min(max(round(slope.*d_new(:) + b),1),nConf);
% %         end
% %         for iX = 1:nX;
% %             newHist = histc(newHisttemp(:,iX),1:nConf);
% %             pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
% %             LL_new(iX) = nnew_part*log(pnew);
% %         end
% %         
% %         % old word match and mismatch
% %         idxmatch = bsxfun(@eq, SOld, X);
% %         mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0), ~idxmatch),2);
% %         
% %         Xrep = repmat(X,[1 1 Nold*nS 1]);
% %         idxmatch = bsxfun(@eq, SOld, X);
% %         matmat = ones(Nold, M, Nold*nS, nX);
% %         matmat(idxmatch) = mismatchoddsVec(Xrep(idxmatch));
% %         
% %         d_old = -log(Nnew) + log(sum((1-c).^mismatchSum.*prod(matmat,2),1));   % log odds
% %         
% %         % binning old words
% %         if (islogbinning) % log binning
% %             oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
% %         else % lin binning
% %             oldHist = max(min(round(slope.*d_old(:) + b),nConf),1);
% %         end
% %         oldHist = histc(oldHist,1:nConf); % histogram
% %         pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
% %         
% %         % calculating nLL
% %         LL_old = nold_part*log(pold);
% %         LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
% %         nLL = -LL_new-LL_old;
% %         
% % %          % indices in which new words match X
% % %         matmat = ones(Nold, M, Nold*nS, nX);
% % %         all_idx = mod(find(idxmatch),Nold*M*nX);
% % %         all_idx(all_idx==0) = Nold*M*nX;
% % %         matmat(idxmatch) = mismatchoddsVec(X(all_idx));
% % %         
% % %         % indices 2 (permuting the 3rd and 4th)
% % %         tic;
% % %         idxmatch = bsxfun(@eq, SNew, X);
% % %         matmat2 = ones(Nold, M, nX, Nold*nS); % 1 2 4 3
% % %         idxmatch = permute(idxmatch, [1 2 4 3]);
% % %         all_idx = mod(find(idxmatch),Nold*M*nX);
% % %         all_idx(all_idx==0) = Nold*M*nX;
% % %         matmat2(idxmatch) = mismatchoddsVec(X(all_idx));
% % %         matmat2 = permute(matmat2, [1 2 4 3]);
% % %         toc
% %         
        
        
        % ====== FOR LOOP =======
        for iX = 1:nX;
            X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
            Xrepp = repmat(X,[nS 1]);
            
            SNew = geornd(g,[Nnew*nS M])+1; % new words
            
            % old words
            idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
            SOld = (1-idx).*Xrepp + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
            
            % decision variable values for new words
            mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0),bsxfun(@ne,permute(SNew, [3 2 1]),X)),2);
            
            idxmatch = bsxfun(@eq, permute(SNew, [3 2 1]), X); % indices in which new words match X
            matmat = ones(Nold, M, Nold*nS);
            all_idx = mod(find(idxmatch),Nold*M);
            all_idx(all_idx==0) = Nold*M;
            matmat(idxmatch) = mismatchoddsVec(X(all_idx));
           
            d_new = -log(Nnew) + log(sum((1-c).^mismatchSum.*prod(matmat,2),1));   % log odds
            
            
            
            
            % decision variable values for new words
            mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0),bsxfun(@ne,permute(SOld, [3 2 1]),X)),2);
            
            idxmatch = bsxfun(@eq, permute(SOld, [3 2 1]), X); % indices in which new words match X
            matmat = ones(Nold, M, Nold*nS);
            all_idx = mod(find(idxmatch),Nold*M);
            all_idx(all_idx==0) = Nold*M;
            matmat(idxmatch) = mismatchoddsVec(X(all_idx));
           
            d_old(:,iX) = -log(Nold) + log(sum((1-c).^mismatchSum.*prod(matmat,2),1));   % log odds
            
            
%             % decision variable values for old words
%             Xrep = repmat(X,[1 1 Nold*nS]);
%             mismatchSum = sum(permute(repmat(SOld, [1 1 Nold]),[3 2 1]) ~= Xrep & repmat(X~= 0, [1 1 Nnew*nS]),2);
% %            mismatch = sum(bsxfun(@and,bsxfun(@ne,X,0),bsxfun(@ne,permute(SOld, [3 2 1]),X)),2);
%                         
%             idxmatch = permute(repmat(SOld,[1 1 Nold]),[3 2 1]) == Xrep;
%             matmat = ones(Nold, M, Nold*nS);
%             matmat(idxmatch) = mismatchoddsVec(Xrep(idxmatch));
%             d_old(:,iX) = -log(Nold) + log(sum((1-c).^mismatchSum.*prod(matmat,2),1));   % log odds
            
            
            % binning new words.
            if (islogbinning)
                newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
            else
                newHisttemp = min(max(round(slope.*d_new(:) + b),1),nConf);
            end
            newHist = histc(newHisttemp,1:nConf);
            pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
            LL_new(iX) = nnew_part*log(pnew);
        end
        
        
        % binning old words
        if (islogbinning) % log binning
            oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
        else % lin binning
            oldHist = max(min(round(slope.*d_old(:) + b),nConf),1);
        end
        oldHist = histc(oldHist,1:nConf); % histogram
        pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
        
        % calculating nLL
        LL_old = nold_part*log(pold);
        LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
        nLL = -LL_new-LL_old;
end

