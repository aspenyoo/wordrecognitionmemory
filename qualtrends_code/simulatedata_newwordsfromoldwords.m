function [nnew_part, nold_part, dNew, dOld] = simulatedata_newwordsfromoldwords(modelname, theta, islogbinning, nX, nS, nConf, N, plotstuff)
% simulatedata_newwordsfromoldwords simulates fake data for a MODELNAME and
% THETA, in which the fifth element is the SD from which each new word is
% drawn around an old word. 

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

if strcmp(modelname,'uneqVar') % if UVSD model
    [pnew, pold] = responses_uneqVar(theta, islogbinning);
    nnew_part = round(Nnew*pnew);
    nold_part = round(Nold*pold);
    
else % if a optimal or heuristic model
    M = theta(1);                   % number of features
    sigma = theta(2);               % memory noise
    k = theta(3);                   % slope of logistic binning function
    d0 = theta(4);                  % shift of logistic binning function
    sigma_new = theta(5);           % how wide of a distribution you draw new words from
    L = nConf/2;
    
    sigs = 1;                       % width of word feature distribution
    J = 1/sigs^2 + 1/sigma^2;
    
    newHist = nan(nX,nConf);
    dOld = nan(Nold,nX,nS);
    d_new = nan(Nnew,nS);
    dNew = nan(Nold,nX,nS);
    for iX = 1:nX;
        X = randn(Nold, M)*sqrt(sigma^2+1);
        
        for iS = 1:nS;
            SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
            SNew = SOld + randn(Nnew,M).*sigma_new; % drawing new words
            
            switch modelname
                case 'FP'
                    d_new(:,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
                    dOld(:,iX,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
                case 'FPheurs'
                    d_new(:,iS) = squeeze(-min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1));
                    dOld(:,iX,iS) = squeeze(-min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2),[],1));
            end
        end
        dNew(:,iX,:) = d_new;
        
        % binning new words.
        if (islogbinning)
            newHisttemp = min(round(L+0.5+ L.*(2./(1+exp(-(d_new(:)-d0)./k)) - 1)),nConf);
        else
            newHisttemp = min(max(round(m.*d_new(:) + b),1),nConf);
        end
        newHist(iX,:) = histc(newHisttemp,1:nConf);
    end
    % figure;
    % minbound = min([d_new(:)' d_old(:)']);
    % maxbound = max([d_new(:)' d_new(:)']);
    % interv = linspace(minbound,maxbound,50);
    % memDist_new = hist(d_new(:),interv);
    % memDist_old = hist(d_old(:),interv);
    % plot(interv,memDist_new,interv,memDist_old./(nX));
    
    
    % binning old words
    if (islogbinning) % log binning
        oldHist = min(round(L.*(2./((1+exp(-(dOld(:)-d0)./k))) - 1)+10.5),nConf);
    else % lin binning
        oldHist = max(min(round(m.*dOld(:) + b),nConf),1);
    end
    oldHist = histc(oldHist,1:nConf); % histogram
    
    nnew_part = round(mean(newHist,1)/(nS));
    nold_part = round(oldHist'/(nX*nS));
end




% % ========== PLOTTING STUFF ==========
% ===== PLOTS =====

if plotstuff(2) == 1;
    if plotstuff(1) ==1;
        subplot(2,1,1)
    end
    if strcmp(modelname,'uneqVar')
       responses_uneqVar( theta, islogbinning, 0, 1);
    else
        binningparameters = theta(3:end);
        memdistplot(dOld(:), d_new(:), binningparameters, islogbinning);
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