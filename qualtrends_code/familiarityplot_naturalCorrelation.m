% function [nnew_part, nold_part] = simulate_data(modelname, theta, islogbinning, nX, nS, nConf, N, plotstuff)
% simulate_data simulates fake data

modelname = 'FP';
theta = [20 1 4 0 nan];
islogbinning = 1; 
nX = 1; 
nS = 1;
nConf = 20;
N = [150 150];
Nold = N(1); Nnew = N(2);

if strcmp(modelname,'uneqVar') % if UVSD model
    [pnew, pold] = responses_uneqVar(theta, islogbinning);
    nnew_part = round(Nnew*pnew);
    nold_part = round(Nold*pold);
    
else % if a optimal or heuristic model
    M = theta(1);                   % number of features
    sigma = theta(2);               % memory noise
    k = theta(3);                   % slope of logistic binning function
    d0 = theta(4);                  % shift of logistic binning function
    L = nConf/2;
    
    sigs = 1;                       % width of word feature distribution
    J = 1/sigs^2 + 1/sigma^2;
    
%     newHist = nan(nX,nConf);
    dOld = nan(Nold,nX,nS);
    d_new = nan(Nnew,nS);
    dNew = nan(Nold,nX,nS);
    dist_new_X = nan(Nnew,nS); dist_new_old = dist_new_X;
    distNewX = nan(Nold,nX,nS); distNewOld = distNewX;
    distOld = nan(Nold,nS,nS);
    for iX = 1:nX;
        X = randn(Nold, M)*sqrt(sigma^2+1);
        
        for iS = 1:nS;
            SNew = randn(Nnew,M); % drawing new words
            SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(Nold,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
            
            switch modelname
                case 'FP'
                    d_new(:,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew])).^2,2)))));
                    dOld(:,iX,iS) = M/2*log(1+sigs^2/sigma^2) + 0.5*sum(SOld.^2,2) + log(squeeze(mean(exp(-0.5*J*...
                        sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nold])).^2,2)))));
                    dist_new_X(:,iS) = sqrt(squeeze(min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1))); % distance bewteen new words and the closest noisy memory
                    distOld(:,iX,iS) = sqrt(squeeze(min(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2),[],1))); % distance between old words and the closest noisy memory
                    dist_new_old(:,iS) = sqrt(squeeze(min(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(SOld,[1,1,Nnew])).^2,2),[],1))); % distance bewteen new words and the closest old word
                case 'FPheurs'
                    d_new(:,iS) = squeeze(-min(sqrt(sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nnew])).^2,2)),[],1));
                    dOld(:,iX,iS) = squeeze(-min(sqrt(sum((permute(repmat(SOld,[1,1,Nold]), [3,2,1])-repmat(X,[1,1,Nold])).^2,2)),[],1));
            end
        end
        dNew(:,iX,:) = d_new;
        distNewX(:,iX,:) = dist_new_X;
        distNewOld(:,iX,:) = dist_new_old;
    end

    
    % binning words
    if (islogbinning) % log binning
        confOld = min(round(L.*(2./(1+exp(-(dOld(:)-d0)./k)) - 1)+10.5),nConf);
        confNew = min(round(L.*(2./(1+exp(-(dNew(:)-d0)./k)) - 1)+L+0.5),nConf);
    else % lin binning
        confOld = max(min(round(m.*dOld(:) + b),nConf),1);
        confNew = min(max(round(m.*d_new(:) + b),1),nConf);
    end
    
    % statistics: correlation
    [rho_newX, pvalue_newX] = corr([distNewX confNew]); % new word: distance from noisy memory by confidence
    [rho_old, pvalue_old] = corr([distOld confOld]); % old word: distance from noisy memory by confidence
    [rho_newOld, pvalue_newOld] = corr([distNewOld confNew]); % new word: distance from old word by confidence
    [rho_newlpr, pvalue_newlpr] = corr([distNewX dNew]); % new word: distance by LPR
    [rho_oldlpr, pvalue_oldlpr] = corr([distOld dOld]); % old word: distance by LPR
    rho_newX = rho_newX(2); rho_old = rho_old(2); rho_newlpr = rho_newlpr(2); rho_oldlpr = rho_oldlpr(2);
    rho_newOld = rho_newOld(2); pvalue_newOld = pvalue_newOld(2);
    pvalue_newX = pvalue_newX(2); pvalue_old = pvalue_old(2);
    pvalue_newlpr = pvalue_newlpr(2); pvalue_oldlpr = pvalue_oldlpr(2);

    % statistics: regression
    [b_newX,~,~,~,stats_newX] = regress(confNew,[ones(length(confNew),1) distNewX]);
    [b_newlpr,~,~,~,stats_newlpr] = regress(dNew,[ones(length(dNew),1) distNewX]);
    [b_newOld,~,~,~,stats_newOld] = regress(confNew,[ones(length(confNew),1) distNewOld]);
    bestfitline_newX = [1 min(distNewX); 1 max(distNewX)]*b_newX; %[ones(length(confNew),1) distNew]*b_new;
    bestfitline_newlpr = [1 min(distNewX); 1 max(distNewX)]*b_newlpr; %[ones(length(confNew),1) distNew]*b_new;
    bestfitline_newOld = [1 min(distNewOld); 1 max(distNewOld)]*b_newOld; %[ones(length(confNew),1) distNew]*b_new;
    
    % ========== PLOTTING STUFF ==========
    % plotting distance as a function of distance from closest noisy memory
    
    grey = 0.7*ones(1,3);
    
    % new word: confidence rating by distance to closest noisy memory
    figure; 
    plot(distNewX,confNew,'.','Color',aspencolors('greyblue'),'MarkerSize',12); hold on
    plot([min(distNewX) max(distNewX)],bestfitline_newX,'Color',grey);
    xlabel('distance from closest noisy memory')
    ylabel('confidence rating')
    if pvalue_newX < 0.001;
        signn = '<';
        pvalue_newX = 0.001;
    else
        signn = '=';
    end
    tt = annotation('textbox',[.6 .8 .3 .1],...
               'String', sprintf('\\rho = %2.3f, p %c %2.3f, \nR^2 = %2.3f', rho_newX, ...
               signn, pvalue_newX, stats_newX(1)));
           set(tt,'LineStyle','none')
    defaultplot
    axis([3 7 1 20])
    ax = gca;
    ax.YTick = [ 1 10 20];
    ax.XTick = [3 4 5 6 7];
    
    % new word: confidence rating by distance to closest old word
    figure; 
    plot(distNewOld,confNew,'.','Color',aspencolors('greyblue'),'MarkerSize',12); hold on
    plot([min(distNewOld) max(distNewOld)],bestfitline_newOld,'Color',grey);
    xlabel('distance from closest old word')
    ylabel('confidence rating')
    if pvalue_newOld < 0.001;
        signn = '<';
        pvalue_newOld = 0.001;
    else
        signn = '=';
    end
    tt = annotation('textbox',[.6 .8 .3 .1],...
               'String', sprintf('\\rho = %2.3f, p %c %2.3f, \nR^2 = %2.3f', rho_newOld, ...
               signn, pvalue_newOld, stats_newOld(1)));
           set(tt,'LineStyle','none')
    defaultplot
    axis([3 7 1 20])
    ax = gca;
    ax.YTick = [ 1 10 20];
    ax.XTick = [3 4 5 6 7];
    
%     % old word: confidence rating by distance to closest old word
%     figure; 
%     plot(distOld,confOld,'.','Color',aspencolors('dustygold')); hold on
%     plot([min(distOld) max(distOld)],bestfitline_old,'Color',grey);
%     ylim([1 20]);
%     xlabel('distance from closest old word')
%     ylabel('confidence rating')
%     tt = annotation('textbox',[.6 .8 .3 .1],...
%                'String', sprintf('\\rho = %2.3f, p = %2.3f, \nR^2 = %2.3f', rho_old, ...
%                pvalue_old, stats_old(1)));
%            set(tt,'LineStyle','none')
%     defaultplot
    
%     % new word: LPR by distance to closest old word
%     figure; 
%     plot(distNew,dNew,'.','Color',aspencolors('greyblue')); hold on
%     plot([min(distNew) max(distNew)],bestfitline_newlpr,'Color',grey);
%     xlabel('distance from mean old word')
%     ylabel('LPR')
%     tt = annotation('textbox',[.6 .8 .3 .1],...
%                'String', sprintf('\\rho = %2.3f, p = %2.3f, \nR^2 = %2.3f', rho_newlpr, ...
%                pvalue_newlpr, stats_newlpr(1)));
%            set(tt,'LineStyle','none')
%     defaultplot

end

