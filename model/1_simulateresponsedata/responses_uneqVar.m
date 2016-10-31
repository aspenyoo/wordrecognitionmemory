function [pnew, pold, confBounds] = responses_uneqVar( theta, nConf, cplot, dplot)
% [pnew, pold] = function responses_uneqVar( THETA, PPLOT, RPLOT) gives
% probability of "old" and "new" response confidences.
% 
% [pnew, pold] = function responses_uneqVar(THETA) outputs the
% probability of responding at difference confidences given a group of
% parameters (THETA).
% CPLOT = 1 plots pnew and pold
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

d0 = theta(end-4);
a = theta(end-3);
b = theta(end-2);
gamma = theta(end-1);
k = theta(end);
nParams = 7;


assert(nParams == length(theta),'check number of parameters');

% if any(binningfn == [2 3])
%     d_new_sign = sign(d_new(:)+d0); % -1 for respond new, +1 for respond old
%     q = abs(d_new(:)+d0);
% else
%     q = d_new(:) + d0;
% end

% calculate confBounds
% ratingBounds = [ -Inf 1.5:19.5 Inf];
binvalues = 1.5:(nConf/2 - 0.5);
decisionboundary = fminsearch(@(x) abs(normpdf(x,mu_old,sigma_old) - normpdf(x)),rand); % finding decisionboundary and shifting the distributions over by it so that 0 is decision boundary
tempp = 1-((binvalues-b)./a);
tempp(tempp < 0) = nan;
tempp(tempp > 1) = nan;
confBounds = gamma.*(-log(tempp)).^(1/k) - d0;
confBounds
confBounds = [-Inf decisionboundary-confBounds decisionboundary decisionboundary+confBounds Inf];
confBounds

% % PNEW AND POLD

% getting "new" response confidence
confbounds = confBounds;
confbounds(confbounds < 0) = 0;
confbounds(isnan(confbounds)) = 0;
confbounds = [-Inf decisionboundary-confbounds decisionboundary];
pnew = normcdf(confbounds(2:end)) - normcdf(confbounds(1:end-1));
pold = normcdf(confbounds(2:end),mu_old,sigma_old) - normcdf(confbounds(1:end-1),mu_old,sigma_old);

% getting "old" response confidence



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


