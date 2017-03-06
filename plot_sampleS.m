% plot_sampleS

clear all

modelname = 'FP';
binningfn = 4;
nS = 300;

load('subjdata.mat')
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'])
isubj = 1;
realtheta = 0;

% getting model specific parameters
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

% real participant data
nnew_part = nNew_part(isubj,:);
nold_part = nOld_part(isubj,:);
nConf = 20;

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

switch binningfn
    case 3 % power law
        a = theta(end-4);
        b = theta(end-3);
        gamma = theta(end-2);
    case 4 % weibull
        scale = theta(end-5);
        shape = theta(end-4);
        a = theta(end-3);
        b = theta(end-2);
end
d0 = theta(end-1);
sigma_mc = theta(end);

[Nold, Nnew] = deal(150);

figure
nSamples = 10;
for isample = 1:nSamples
    isample
    switch modelname
        case 'FP'
            X = randn(Nold, M)*sqrt(sigma^2+1);
            SNew = randn(Nnew*nS,M);
            SOld = (repmat(X,[nS 1])/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
            
            [d_new, d_old] = calculate_d_FP(M, sigma, nS, Nnew, Nold, SNew, SOld, X);
        case 'REM'
            X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
            SNew = geornd(g,[Nnew*nS M])+1; % new words
            idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
            SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
            
            maxidx = max(X(:));
            matchoddsVec = (c+(1-c).*g.*(1-g).^(0:maxidx))./(g.*(1-g).^(0:maxidx));
            [d_new, d_old] = calculate_d_REM(M, g, c, nS, Nnew, Nold, SNew, SOld, X, matchoddsVec);
    end
    
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
    
    newHist = [zeros(length(conf),nConf/2) 0.5+0.5.*erf(bsxfun(@minus,[1.5:(nConf/2-0.5) Inf],conf)./(sigma_mc*sqrt(2))) - ...
        (0.5+0.5.*erf(bsxfun(@minus,[-Inf 1.5:(nConf/2-0.5)],conf)./(sigma_mc*sqrt(2))))];
    newHist(d_new_sign < 0,:) = fliplr(newHist(d_new_sign < 0,:));
    
    lapse = 0.01;             % lapse rate
    for iS = 1:nS
        
        currnewHist = sum(newHist(1:Nnew*iS,:))';
        
        pnew = lapse/nConf + (1-lapse)*(currnewHist/sum(currnewHist));
        LL_new(iS) = nnew_part*log(pnew);
        finalLL(iS) = max(LL_new) + log(mean(exp(LL_new-max(LL_new))));
        
    end
    
    % clf;
    % plot(LL_new,'.','Color',0.7*ones(1,3))
    hold on
    plot(finalLL,'k');%'Color',0.7*ones(1,3))
    
end
defaultplot
xlabel('number of S samples')
ylabel('LL for lures')