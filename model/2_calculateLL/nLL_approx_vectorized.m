function [ nLL ] = nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part, fixparams, nX, nS, nConf )
% nLL_approx calculates the negative log likelihood using an approximation
% method
%
% [nLL] = nLL_approx_vectorized(MODELNAME, THETA, ISLOGBINNING, NNEW_PART,
% NOLD_PART) calculates the negative log p(nnew_part,nold_part|theta,model)
%
% ===== INPUT VARIABLES =====
% MODELNAME: 'FP','FPheurs','VP','VPheurs','uneqVar', 'REM'
% THETA: parameter values
% BINNINGFN: 0: linear, 1: logistic, 2: log
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
% Aspen Yoo - June 29, 2016

if nargin < 6; fixparams = []; end
if nargin < 7; nX = 30; end
if nargin < 8; if strcmp(modelname,'REM'); nS = 20; else nS = 20; end, end
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
        [pnew, pold] = responses_uneqVar(theta, binningfn);
        nLL = -sum(log(pnew).*nnew_part) - sum(log(pold).*nold_part);
    case {'FP','FPheurs'}
        M = theta(1);
        sigma = theta(2);
        switch binningfn
            case 0 % linear
                c1 = theta(3);
                c2 = theta(4); 
            case 1 % logistic
                 k = theta(3); 
                 d0 = theta(4); 
            case 2 % log
                a = theta(3);
                b = theta(4);
                sigma_mc = theta(5);
        end
        Nold = sum(nold_part); Nnew = sum(nnew_part);
        L = nConf/2;
        J = 1/sigma.^2+1;
        lambda = 0.01;             % lapse rate
        
        if M ~= floor(M)
            M = floor(M);
            %     assert('M must be a whole number')
        end
        
        if ~(binningfn) % if linear
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
            
            %             % plot so see how nS affects data
            %             for iS = 1:nS;
            %                 % binning new words.
            %                 if (islogbinning)
            %                     newHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(1:iS*Nold)-d0)./k)) - 1)),nConf); % bounds: [1 20]
            %                 else
            %                     newHist = min(max(round(m.*d_new(:) + b),1),20); % bounds: [1 20]
            %                 end
            %                 newHist = histc(newHist,1:nConf);
            %                 pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
            %                 LL_newtemp(iS) = nnew_part*log(pnew);
            %             end
            %             figure;
            %             plot(1:nS,LL_newtemp,'ok')
            %             defaultplot
            
            % binning new words.
            switch binningfn
                case 0 % linear
                    newHist = min(max(round(m.*d_new(:) + b),1),nConf);                            % bounds: [1 20]
                case 1 % logistic
                    newHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);   % bounds: [1 20]
                case 2 % log
                    d_new_sign = sign(d_new(:));                                                % -1 for respond new, +1 for respond old
                    newHist = min(max(round(a.*log(abs(d_new(:)+b))+randn.*sigma_mc),nConf/2+1),nConf);        % confidence rating from 11 to 20
                    newHist(d_new_sign < 0) = nConf+1 - newHist(d_new_sign < 0);                     % changing respond "new" words back to 1 to 10
            end
            newHist = histc(newHist,1:nConf);
            pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
            LL_new(iX) = nnew_part*log(pnew);
            
        end
        
        
%         % plotting for X samples
%         for iX = 1:nX;
%             % BINNING OLD WORDS. p(conf|X0,C)
%             oldHist = min(round(L.*(2./((1+exp(-(d_old(:,1:iX)-d0)./k))) - 1)+10.5),nConf);
%             oldHist = histc(oldHist,1:nConf); % histogram
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
%         defaultplot
%         xlabel('nX')
%         ylabel('LL for new words')
%         title('M = 12, \sigma = 1.02, k = 0.31, d0 = 0.1')
%         legend('mean LL','max of samples','LL for each sample')
        
        % BINNING OLD WORDS. p(conf|X0,C)
        switch binningfn
            case 0 % linear
                oldHist = min(max(round(m.*d_old(:) + b),1),nConf);                            % bounds: [1 20]
            case 1 % logistic
                oldHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_old(:)-d0)./k)) - 1)),nConf);   % bounds: [1 20]
            case 2 % log
                d_old_sign = sign(d_old(:));                                                % -1 for respond new, +1 for respond old
                oldHist = min(max(round(a.*abs(d_old(:)+b)+randn.*sigma_mc),nConf/2+1),nConf);        % confidence rating from 11 to 20
                oldHist(d_old_sign < 0) = nConf+1 - oldHist(d_old_sign < 0);                     % changing respond "new" words back to 1 to 10
        end
        oldHist = histc(oldHist,1:nConf); % histogram
        % oldHist(oldHist==0) = 1e-3; % changing any 0 freq to 1, (prevents LL from going to -Inf)
        pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
        
        % calculating nLL
        LL_old = nold_part*log(pold);
        LL_new = max(LL_new) + log(mean(exp(LL_new-max(LL_new)))); % average over X
        nLL = -LL_new-LL_old
        
    case 'REM'
        
        M = theta(1);                   % number of features
        g = theta(2);                   % probability of success (for geometric distribution
        ustar = theta(3);               % probability of encoding something
        c = theta(4);                   % probability of encoding correct feature value
        m = theta(5);                   % number of storage attempts
        L = nConf/2;
        if (binningfn); 
            k = theta(6);  d0 = theta(7);
        else c1 = theta(6); c2 = theta(7);
        end
        Nold = sum(nold_part); Nnew = sum(nnew_part);
        lambda = 0.01;             % lapse rate
        
        if M ~= floor(M)
            M = floor(M);
        end
        
        if ~(binningfn)
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
        
        % looping over X samples
        for iX = 1:nX;
            
            X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);

            % generating new and old test words
            SNew = geornd(g,[1 M Nnew*nS])+1; % new words
            idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
            SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories

            % decision variable for new words
            idxmatch = bsxfun(@eq, SNew, X); % indices in which new words match X

            LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,(1-g).^-X) ; 
            LRmat(bsxfun(@eq,zeros(1,M,Nnew*nS),X)) = 1;
            d_new = log(mean(prod(LRmat,2),1));   % log odds
            
            
            %             % plot so see how nS affects data
%             for iS = 1:nS;
%                 % binning new words.
%                 if (islogbinning)
%                     newHist= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(1:iS*Nold)-d0)./k)) - 1)),nConf); % bounds: [1 20]
%                 else
%                     newHist = min(max(round(m.*d_new(:) + b),1),20); % bounds: [1 20]
%                 end
%                 newHist = squeeze(histc(newHist,1:nConf));
%                 pnew = lambda/nConf + (1-lambda).*(newHist/sum(newHist));
%                 LL_newtemp(iS) = nnew_part*log(pnew);
%             end
%             subplot(3,3,iB);
%             plot(1:nS,LL_newtemp,'ok')
%             defaultplot
%             xlabel('nS')
%             ylabel('LLnew')
            
            % decision variable values for old words
            idxmatch = bsxfun(@eq, permute(SOld, [3 2 1]), X); % indices in which new words match X
            
            LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,exp(-X.*log(1-g))); 
            LRmat(bsxfun(@eq,zeros(1,M,Nnew*nS),X)) = 1;
            d_old(:,iX) = log(mean(prod(LRmat,2),1));
            

            % binning new words.
            if (binningfn)
                newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
            else
                newHisttemp = min(max(round(slope.*d_new(:) + b),1),nConf);
            end
            newHist = histc(newHisttemp,1:nConf);
            pnew = lambda/nConf + (1-lambda)*(newHist/sum(newHist));
            LL_new(iX) = nnew_part*log(pnew);
        end
        
%         % plotting for X samples
%         for iX = 1:nX;
%             % BINNING OLD WORDS. p(conf|X0,C)
%             oldHist = min(round(L.*(2./((1+exp(-(d_old(:,1:iX)-d0)./k))) - 1)+10.5),nConf);
%             oldHist = histc(oldHist,1:nConf); % histogram
%             pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
%             
%             % calculating nLL
%             LL_oldtemp(iX) = nold_part*log(pold);
%             LL_newtemp(iX) = max(LL_new(1:iX)) + log(mean(exp(LL_new(1:iX)-max(LL_new(1:iX))))); % average over X
%             maxnew(iX) = max(LL_new(1:iX));
%             %             blah(iX) = max(LL_new(iX))-LL_new(1:iX);
%         end
%         nLLtemp = -LL_newtemp-LL_oldtemp;
%                 figure; plot(1:nX,nLLtemp,'ko');defaultplot
%         %         hold on;
%         %         plot(1:nX,-LL_oldtemp,'bo');
%         
%         figure; hold on;
%         plot(1:nX,LL_newtemp,'o','Color',aspencolors('gold'));
%         plot(1:nX,LL_new,'*','Color',aspencolors('dustygold'));
%         plot(1:nX,maxnew,'.')
%         defaultplot
%         xlabel('nX')
%         ylabel('LL for new words')
% %         title('M = 12, \sigma = 1.02, k = 0.31, d0 = 0.1')
%         legend('mean LL','max of samples','LL for each sample')

        % binning old words
        if (binningfn) % log binning
            oldHist = min(round(L.*(2./((1+exp(-(d_old(:)-d0)./k))) - 1)+10.5),nConf);
        else % lin binning
            oldHist = max(min(round(slope.*d_old(:) + b),nConf),1);
        end
        oldHist = histc(oldHist,1:nConf); % histogram
        pold = lambda/nConf + (1-lambda)*(oldHist/sum(oldHist)); % normalizing
        
        % calculating nLL
        LL_old = nold_part*log(pold);
        LL_new = max(LL_new) + log(mean(exp(max(LL_new)-LL_new))); % average over X
        nLL = -LL_new-LL_old;
end

