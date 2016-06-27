function [nnew_part, nold_part] = simulate_data(modelname, theta, islogbinning, nX, nS, nConf, N, plotstuff)
% simulatefake simulates fake data

if nargin < 3; islogbinning = 1; end
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
        [pnew, pold] = responses_uneqVar(theta, islogbinning);
        nnew_part = round(Nnew*pnew);
        nold_part = round(Nold*pold);
        
    case {'FP','FPheurs'} % optimal or heuristic model
        M = theta(1);                   % number of features
        sigma = theta(2);               % memory noise
        k = theta(3);                   % slope of logistic binning function
        d0 = theta(4);                  % shift of logistic binning function
        L = nConf/2;
        
        sigs = 1;                       % width of word feature distribution
        J = 1/sigs^2 + 1/sigma^2;
        
        newHist = nan(nX,nConf);
        d_old = nan(Nold,nX,nS);
        d_new = nan(Nnew,nS);
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
                newHisttemp= min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf); % bounds: [1 20]
            else
                newHisttemp = min(max(round(m.*d_new(:) + b),1),20); % bounds: [1 20]
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
        mismatchoddsVec = (c+(1-c).*(g.*(1-g).^(0:20)))./(g.*(1-g).^(0:20));
        
        newHist = nan(nX,nConf);
        d_old = nan(Nold*nS,nX);
        for iX = 1:nX;
            X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
            Xrepp = repmat(X,[nS 1]);
            
            SNew = geornd(g,[Nnew*nS M])+1; % new words
            
            % old words
            idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
            SOld = (1-idx).*Xrepp + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
            
            % decision variable values for new words
            Xrep = repmat(X,[1 1 Nnew*nS]);            % expanded X
            mismatch = sum(permute(repmat(SNew, [1 1 Nold]),[3 2 1]) ~= Xrep & repmat(X~= 0, [1 1 Nnew*nS]),2); % number of mismatching features
            
            idxmatch = permute(repmat(SNew,[1 1 Nold]),[3 2 1]) == Xrep; % indices in which new words match X
            matmat = ones(Nold, M, Nold*nS);
            matmat(idxmatch) = mismatchoddsVec(Xrep(idxmatch));
            d_new = -log(Nnew) + log(sum((1-c).^mismatch.*prod(matmat,2),1));   % log odds
            
            % decision variable values for old words
            Xrep = repmat(X,[1 1 Nold*nS]);
            idxmatch = permute(repmat(SOld,[1 1 Nold]),[3 2 1]) == Xrep;
            mismatch = sum(permute(repmat(SOld, [1 1 Nold]),[3 2 1]) ~= Xrep & repmat(X~= 0, [1 1 Nnew*nS]),2);
            
            matmat = ones(Nold, M, Nold*nS);
            matmat(idxmatch) = mismatchoddsVec(Xrep(idxmatch));
            
            d_old(:,iX) = -log(Nold) + log(sum((1-c).^mismatch.*prod(matmat,2),1));   % log odds

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




% % ========== PLOTTING STUFF ==========

if plotstuff(2) == 1;
    if plotstuff(1) ==1;
        subplot(2,1,1)
    end
    if strcmp(modelname,'uneqVar')
        responses_uneqVar( theta, islogbinning, 0, 1);
    else
        binningparameters = theta(3:end);
        memdistplot(d_old(:), d_new(:), binningparameters, islogbinning);
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