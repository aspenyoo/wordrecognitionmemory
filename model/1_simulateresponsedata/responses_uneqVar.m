function [pnew, pold, confBounds] = responses_uneqVar( theta, nConf, cplot, dplot)
% [pnew, pold] = function responses_uneqVar( THETA, PPLOT, RPLOT) gives
% probability of "old" and "new" response confidences.
%
% [pnew, pold] = function responses_uneqVar(THETA) outputs the
% probability of responding at difference confidences given a group of
% parameters (THETA).
% CPLOT = 1 plots pnew and poldoepn 
% DPLOT = 1 plots the memory distribution (of d)
%
% theta = [mu_old, sigma_old, k, d0, gamMax] for log binning  and
% theta = [mu_old, sigma_old, c1, c2] for linear binning
%
% Aspen Yoo
% march 6, 2015

% theta = [2 1.25 1 1.5 10];
% pplot = 1;
% rplot = 1;

% -----------------------------------------------------------------------
if nargin < 2; nConf = 20; end
if nargin < 3; cplot = 0;end
if nargin < 4; dplot = 0; end

% getting parameter names
mu_old = theta(1);
sigma_old = theta(2);

switch binningfn
    case {0,1}              % linear (0) and logistic (1) mapping
        k = theta(end-1);
        d0 = theta(end);
        nParams = 4;
    case 2        % logarithmic mapping on d (not other x axis)
        a = theta(end-3);
        b = theta(end-2);
        d0 = theta(end-1);
        sigma_mc = theta(end);
        nParams = 6;
    case 3              % power law mapping from p(correct)
        a = theta(end-4);
        b = theta(end-3);
        d0 = theta(end-2);
        lambda = theta(end-1);
        sigma_mc = theta(end);
        nParams = 6;
    case 4              % weibull mapping
        scale = theta(end-3);
        shape = theta(end-2);
        d0 = theta(end-1);
        sigma_mc = theta(end);
        nParams = 6;
end

assert(nParams == length(theta),'check number of parameters');

if ~(sigma_mc)
    ratingBounds = [ -Inf 1.5:(nConf - 0.5) Inf];
    decisionboundary = fminsearch(@(x) abs(normpdf(x,mu_old,sigma_old) - normpdf(x)),rand); % finding decisionboundary and shifting the distributions over by it so that 0 is decision boundary
    switch binningfn
        case 0 % linear
            confBounds = (ratingBounds-d0-nConf/2-0.5)./k +decisionboundary;
        case 1
            confBounds = -k.*log(nConf./(ratingBounds - 0.5)-1)+d0+decisionboundary;
%         case 2 % logarithmic
%             ratingBounds = [ 0 1.5:(nConf/2 - 0.5) Inf];
%             confBounds = exp((ratingBounds-b)./a)-d0;
%             confBounds(confBounds < 0) = 0;
%             confBounds = [-confBounds(11:-1:2) confBounds]+decisionboundary;
%         case 3 % power law
% % %         case 4 % weibull %% I don't think this is correct
%             ratingBounds = [ 0 1.5:(nConf/2 - 0.5) Inf];
%             confBounds = -scale/shape.*log(ratingBounds - 1);
%             confBounds(confBounds < 0) = 0;
%             confBounds = [-confBounds(11:-1:2) confBounds]+decisionboundary;
    end
    
    try
        switch binningfn
            case {0,1}
                pnew = normcdf(confBounds(2:end)) - normcdf(confBounds(1:end-1));
                pold = normcdf(confBounds(2:end),mu_old,sigma_old) - normcdf(confBounds(1:end-1),mu_old,sigma_old);
            case {2,3}
                pnew = normcdf(confBounds(2:end)) - normcdf(confBounds(1:end-1));
                pold = normcdf(confBounds(2:end),mu_old,sigma_old) - normcdf(confBounds(1:end-1),mu_old,sigma_old);
        end
    catch
        pnew = zeros(1,nConf);
        pold = zeros(1,nConf);
    end
else % logarithmic binning on |d|
    
    nSamples = 100;
    nSDs = 4; % how far out you wanna bin
    
    % finding a good region to numerically integrate over
    x = linspace(min([0 mu_old]-nSDs.*[1 sigma_old]),max([0 mu_old]+nSDs.*[1 sigma_old]),nSamples);
    y = linspace(-nSDs.*sigma_mc,nSDs.*sigma_mc,nSamples);
    
    % calculate joint probability distribution
    pold_x = normpdf(x,mu_old,sigma_old);
    pnew_x = normpdf(x,0,1);
    p_y = normpdf(y,0,sigma_mc);
    jointdist_old = pold_x'*p_y; % first dimension p(old), second p(y)
    jointdist_new = pnew_x'*p_y;
    jointdist_old = jointdist_old/sum(jointdist_old(:)); % normalize
    jointdist_new = jointdist_new/sum(jointdist_new(:));
    
    % put through nonlinearity;
    d = -log(sigma_old) - 1/2.*( ((x-mu_old).^2)./sigma_old.^2 - x.^2);
    [yy,dd] = meshgrid(y,d); % get 2D values of x and y
    switch binningfn
        case 2
            conf = round(a.*log(abs(dd)) + b + yy);
        case 3 % power law
            conf = round(a.*((abs(dd).^lambda - 1)./lambda) + b + yy); 
        case 4
            conf = round(1 - exp(-(abs(dd)./scale).^shape));
    end
    conf(conf < 1) = 1;
    conf(conf > nConf/2) = nConf/2;
    
    % get proportion of responses 
    idx_new = dd < 0; % indices corresponding to "new" responses
    idx_old = dd >= 0; % indices corresponding to "old" responses
    pold = nan(1,nConf);
    pnew = nan(1,nConf);
    for iconf = 1:nConf/2
        
        % indices that correspond to the current confidence value and resp
        idx_new_conf = idx_new & (conf == iconf);
        idx_old_conf = idx_old & (conf == iconf);
        pnew(iconf+nConf/2) = sum(sum(jointdist_new.*idx_old_conf)); 
        pnew(nConf/2 + 1 - iconf) = sum(sum(jointdist_new.*idx_new_conf)); 
        pold(iconf+nConf/2) = sum(sum(jointdist_old.*idx_old_conf)); 
        pold(nConf/2 + 1 - iconf) = sum(sum(jointdist_old.*idx_new_conf)); 
    end
    
    confBounds = nan(1,nConf);
end
pnew = pnew./sum(pnew);
pold = pold./sum(pold);

% making 0 probabilities very small (prevents LL from going to -Inf)
% picked this value because FP and VP model have 150*6*6 total simulations
pnew(pnew==0) = 1.8519e-04;
pold(pold==0) = 1.8519e-04;

% ===== PLOTS =====
nHist = 40;
if dplot == 1;
    
    N = 150;
    gold = aspencolors('dustygold');
    greyblue = aspencolors('greyblue');
    if cplot == 1;
        figure;
        subplot(2,1,1)
    end
    lb = min([-3.5 mu_old-3.5*sigma_old]);
    ub = max([3.5 mu_old+3.5*sigma_old]);
    binbounds2 = linspace(lb,ub,nHist);
    hold on;
    
    centers = (binbounds2(1:end-1) + binbounds2(2:end))/2;
    nelemnew = N* (normcdf(binbounds2(2:end)) - normcdf(binbounds2(1:end-1)));
    nelemold = N*(normcdf(binbounds2(2:end),mu_old,sigma_old) - normcdf(binbounds2(1:end-1),mu_old,sigma_old));
    histcount = [nelemold nelemnew];
    
    maxx = max(histcount) + mod(-max(histcount),30);
    for j = 1:length(ratingBounds);
        hLine = plot([confBounds(j) confBounds(j)], [0 maxx],'Color',[.85 .85 .85]);
        set(get(get(hLine,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend
    end
    
    plot(centers,nelemold,'-','LineWidth',2,'Color',gold);
    plot(centers,nelemnew,'-','LineWidth',2,'Color',greyblue);
    
    ylim([0 maxx])
    defaultplot
    
    set(  gca                         ,...
        'YTick'         , []        ,...
        'YColor'        ,'w'        );
end



if cplot == 1;
    if dplot == 1;
        subplot(2,1,2)
    else
        figure;
    end
    confdistplot(pnew, pold)
end


