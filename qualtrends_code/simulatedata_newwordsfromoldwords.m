function [pnew, pold, d_newtotal, d_old] = simulatedata_newwordsfromoldwords(modelname, theta, binningfn, nX, nS, nConf, N, plotstuff)
% simulatedata_newwordsfromoldwords simulates fake data for a MODELNAME and
% THETA, in which the fifth element is the SD from which each new word is
% drawn around an old word. 

if nargin < 4; nX = 30; end
if nargin < 5; nS = 20; end
if nargin < 6; nConf = 20; end
if nargin < 7; N = [150 150]; end
if nargin < 8; plotstuff = zeros(1,2); end

if isempty(N); N = [150 150]; end
if length(N) == 1; N = [N N]; end
if isempty(nConf); nConf = 20; end
Nnew = N(1); Nold = N(2);       % number of words

if strcmp(modelname,'uneqVar') % if UVSD model
    [pnew, pold] = responses_uneqVar(theta, binningfn);
    
else % if a optimal or heuristic model
    
    switch modelname
        case {'FP','FPheurs','UVSD'}
            nParams = 2;
        case 'REM'
            nParams = 5;
    end
    switch binningfn
        case 0
            slope = theta(end-3);
            nParams = nParams + 4;
        case 1
            k = theta(end-3);
            nParams = nParams + 4;
        case 2 % logarithmic
            a = theta(end-4);
            b = theta(end-3);
            nParams = nParams + 5;
        case 3 % power law
            a = theta(end-5);
            b = theta(end-4);
            gamma = theta(end-3);
            nParams = nParams + 6;
        case 4 % weibull
            scale = theta(end-6);
            shift = theta(end-5);
            a = theta(end-4);
            b = theta(end-3);
            nParams = nParams + 7;
    end
    d0 = theta(end-2);
    sigma_mc = theta(end-1);
    sigma_new = theta(end);           % how wide of a distribution you draw new words from
    
    % check to make sure the length of theta is correct
    assert(nParams == length(theta),'length of theta is not correct')
    
    % stuff that doesn't depend on model, mapping, or memstrengthvar
    Nold = N(1); Nnew = N(2);
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
    
    d_old = nan(Nold*nS,nX);
    d_newtotal = nan(Nnew*nS,nX);
    newHisttotal = nan(nX,nConf);
    for iX = 1:nX
        
        % calculate d
        switch modelname
            case 'FP'
                X = randn(Nold, M)*sqrt(sigma^2+1);
                SOld = (repmat(X,[nS 1])/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
                SNew = SOld + randn(Nnew*nS,M).*sigma_new; % drawing new words
                
                [d_new, d_old(:,iX)] = calculate_d_FP(M, sigma, nS, Nnew, Nold, SNew, SOld, X);
                d_newtotal(:,iX) = d_new;
%             case 'FPheurs'
%                 X = randn(Nold, M)*sqrt(sigma^2+1);
%                 SOld = (repmat(X,[nS 1])/sigma^2)/(1/sigma^2 + 1) + randn(Nold*nS,M)*(1/sqrt(1/sigma^2 + 1));
%                 SNew = SOld + randn(Nnew*nS,M).*sigma_new; % drawing new words
%                 
%                 d_new = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew*nS])).^2,2),[],1));
%                 d_old(:,iX) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold*nS])).^2,2),[],1));
%                 d_newtotal(:,iX) = d_new;
%             case 'REM'
%                 X = binornd(1,1-p0,[Nold M]).*(geornd(g,[Nold M])+1);
%                 idx = logical(binornd(1,pQ,[Nold*nS M]) + repmat((X == 0),[nS 1])); % indices of randomly drawn features
%                 SOld = (1-idx).*repmat(X,[nS 1]) + idx.*(geornd(g,[Nold*nS M]) + 1); % old words from noisy memories
%                 
%                 maxidx = max(X(:));
%                 matchoddsVec = (c+(1-c).*g.*(1-g).^(0:maxidx))./(g.*(1-g).^(0:maxidx));
%                 [d_new, d_old(:,iX)] = calculate_d_REM(M, g, c, nS, Nnew, Nold, SNew, SOld, X, matchoddsVec);
%                 d_newtotal(:,iX) = d_new;
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
                conf = a.*(1-exp(-(q./scale).^shift)) + b;
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
            conf = a.*(1-exp(-(q./scale).^shift)) + b;
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
    

%     L = nConf/2;
%     
%     sigs = 1;                       % width of word feature distribution
%     J = 1/sigs^2 + 1/sigma^2;
%     
%     newHist = nan(nX,nConf);
%     dOld = nan(Nold,nX,nS);
%     d_new = nan(Nnew,nS);
%     dNew = nan(Nold,nX,nS);
%     for iX = 1:nX;
%         X = randn(Nold, M)*sqrt(sigma^2+1);
%         
%         for iS = 1:nS;
%             SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
%             SNew = SOld + randn(Nnew,M).*sigma_new; % drawing new words
%             
%             switch modelname
%                 case 'FP'
%                     d_new(:,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
%                         sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
%                     dOld(:,iX,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
%                         sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
%                 case 'FPheurs'
%                     d_new(:,iS) = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1));
%                     dOld(:,iX,iS) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2),[],1));
%             end
%         end
%         dNew(:,iX,:) = d_new;
%         
%         % binning new words.
%         switch (binningfn)
%             case 0
%                 newHisttemp = min(max(round(m.*d_new(:) + b),1),nConf);
%             case 1
%                 newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
%             case 3
%             case 4
%                 
%         end
%         newHist(iX,:) = histc(newHisttemp,1:nConf);
%     end
%     % figure;
%     % minbound = min([d_new(:)' d_old(:)']);
%     % maxbound = max([d_new(:)' d_new(:)']);
%     % interv = linspace(minbound,maxbound,50);
%     % memDist_new = hist(d_new(:),interv);
%     % memDist_old = hist(d_old(:),interv);
%     % plot(interv,memDist_new,interv,memDist_old./(nX));
%     
%     
%     % binning old words
%     if (binningfn) % log binning
%         oldHist = min(round(L.*(2./((1+exp(-(dOld(:)-d0)./k))) - 1)+10.5),nConf);
%     else % lin binning
%         oldHist = max(min(round(m.*dOld(:) + b),nConf),1);
%     end
%     oldHist = histc(oldHist,1:nConf); % histogram
%     
%     nnew_part = round(mean(newHist,1)/(nS));
%     nold_part = round(oldHist'/(nX*nS));
end




% % ========== PLOTTING STUFF ==========
% ===== PLOTS =====

if plotstuff(2) == 1;
    if plotstuff(1) ==1;
        subplot(2,1,1)
    end
    if strcmp(modelname,'uneqVar')
       responses_uneqVar( theta, binningfn, 0, 1);
    else
        binningparameters = theta(3:end);
        memdistplot(dOld(:), d_new(:), binningparameters, binningfn);
    end
end

if plotstuff(1) == 1;
    if plotstuff(2) == 1;
        subplot(2,1,2)
    else
        figure;
    end
    confdistplot(pnew/N(1), nold_part/N(2))
end

% nHist = 40;
% if plotstuff(1) == 1;
%     gold = aspencolors('gold');
%     greyblue = aspencolors('greyblue');
%     if cplot == 1;
%         figure;
%         subplot(2,1,1)
%     end
%     lb = min([-3.5 mu_old-3.5*sigma_old]);
%     ub = max([3.5 mu_old+3.5*sigma_old]);
%     binbounds2 = linspace(lb,ub,nHist);
%     hold on;
% 
%     centers = (binbounds2(1:end-1) + binbounds2(2:end))/2;
%     nelemnew = N* (normcdf(binbounds2(2:end)) - normcdf(binbounds2(1:end-1)));
%     nelemold = N*(normcdf(binbounds2(2:end),mu_old,sigma_old) - normcdf(binbounds2(1:end-1),mu_old,sigma_old));
%     histcount = [nelemold nelemnew];
% 
%     maxx = max(histcount) + mod(-max(histcount),30);
%     for j = 1:length(ratingBounds);
%         hLine = plot([confBounds(j) confBounds(j)], [0 maxx],'Color',[.85 .85 .85]);
%         set(get(get(hLine,'Annotation'),'LegendInformation'),...
%             'IconDisplayStyle','off'); % Exclude line from legend
%     end
% 
%     plot(centers,nelemold,'-','LineWidth',2,'Color',gold);
%     plot(centers,nelemnew,'-','LineWidth',2,'Color',greyblue);
% 
%     ylim([0 maxx])
%     defaultplot
% 
%         set(  gca                         ,...
%           'YTick'         , []        ,...
%           'YColor'        ,'w'        );
% end
% 
% 
% 
% if plotstuff(2) == 1;
%     if plotstuff(1) == 1;
%         subplot(2,1,2)
%     else
%         figure;
%     end
%     confdistplot(nnew_part/N(1), nold_part/N(2))
% end
% %
% %
% % % if plotstuff(1)  % plot confidence distribution
% % %
% % %
% % % end
% % %
% % %
% % % if plotstuff(2) % plot memory distribution
% % %     memdistplot(d_old(:), dNew(:), theta(3:end), islogbinning);
% % % end