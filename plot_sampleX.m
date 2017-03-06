
clear all

figure;
modelname = 'REM';
binningfn = 3;
isubj = 3;
nS = 50;
nX = 300;
nTimes = 15;
realtheta = 0;

load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'])
switch modelname
    case 'FP'
        %  M, sigma
        plb = [1 1e-3 ];
        pub = [50 3 ];
    case 'REM'
        %  M g ustar c m
        plb = [1 1e-3 1e-3 1e-3 1];
        pub = [50 1 1 1 15];
end
% setting binnfn parameters
switch binningfn
    case 3      % power law
        % a, b, gamma
        plb = [plb 0 -10 -5];
        pub = [pub 10 10 5];
    case 4 % weibull binning
        % scale, shift, a, b
        plb = [plb 0 0 0 -3];
        pub = [pub 10 10 3 3];
end
% d0, sigma_mc
plb = [plb -3 1e-6];
pub = [pub 3 3];

if realtheta
    theta = bestFitParam(isubj,:);
else
    theta = (pub-plb).*rand(1,length(plb))+plb;
    theta(1) = round(theta(1));
    % theta = [30 3 2 6 .4 0 3];
end

% subject data
load('subjdata.mat')
nnew_part = nNew_part(isubj,:);
nold_part = nOld_part(isubj,:);
nConf = 20;

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

for itimes = 1:nTimes;
    itimes
    
    if strcmp(modelname,'UVSD') % if uneqVar
        [pold, pold, confbounds] = responses_uneqVar(theta, binningfn);
        nLL = -sum(log(pold).*nnew_part) - sum(log(pold).*nold_part);
    else % if FP, FPheurs, or REM
        
        % stuff that doesn't depend on model, mapping, or memstrengthvar
        Nold = sum(nold_part); Nnew = sum(nnew_part);
        lapse = 0.01;             % lapse rate
        
        % parameter names
        switch modelname
            case {'FP','FPheurs'}
                M = theta(1);
                sigma = theta(2);
            case 'REM'
                M = theta(1);                   % number of features
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
            
            pold = lapse/nConf + (1-lapse)*(newHist/sum(newHist));
            LL_new(iX) = nnew_part*log(pold);
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
            %         oldHist = (sum(oldHist)./(nS*nX))';
        else
            conf = min(max(round(conf),1),nConf/2)+nConf/2;
            conf(d_old_sign<0) = nConf+1 - conf(d_old_sign<0);
            %         oldHist = histc(conf,1:nConf); % histogram
        end
        
        %     pold = lapse/nConf + (1-lapse)*(oldHist/sum(oldHist)); % normalizing
        %
        %     % calculating nLL
        %     LL_old = nold_part*log(pold);
        %
        %     nLL = -LL_new-LL_old;
        
    end
    
    for iX = 1:nX
        curroldHist = sum(oldHist(1:Nold*iX,:))';
        pold = lapse/nConf + (1-lapse)*(curroldHist/sum(curroldHist));
        
        LL_old(iX) = nold_part*log(pold);
        LL_new(iX) = max(LL_new(1:iX)) + log(mean(exp(LL_new(1:iX)-max(LL_new(1:iX))))); % average over X
    end
    finalLL = LL_new + LL_old;
    
    hold on
    % plot(LL_old,'ro')
    % plot(LL_new,'bo')
    plot(finalLL,'Color',0.7*ones(1,3))
    finalLLMat(itimes,:) = finalLL;
end
MEAN = mean(finalLLMat);
SEM = std(finalLLMat)/sqrt(nTimes);
plot_summaryfit(1:nX,[],[],MEAN,SEM,[],aspencolors('booger'));
hold on
plot(MEAN,'Color',aspencolors('booger'));

defaultplot
xlabel('number of X samples')
ylabel('LL')