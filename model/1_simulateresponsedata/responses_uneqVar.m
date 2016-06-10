function [pnew, pold] = responses_uneqVar( theta, islogbinning, cplot, dplot)
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

if nargin < 4;
    cplot = 0;
    dplot = 0;
end

N = 150;

% getting parameter names
mu_old = theta(1);
sigma_old = theta(2);
if (islogbinning)
    k = theta(3);
    if length(theta) < 5;
        gamMax = 10;
        if length(theta) < 4;
            d0 = 0;
        else
            d0 = theta(4);
        end
    else
        d0 = theta(4);
        gamMax = theta(5);
    end
else
    c1 = theta(3);
    c2 = theta(4);
end

% calculating pnew and pold
if (islogbinning)
    ratingBounds = 0.5:20.5;
    confBounds = -k.*log(((2*gamMax)./(ratingBounds - 10.5 + gamMax)) -1) + d0;
else
    c0 = mean([c1 c2]);
    
    ratingBounds = [-Inf 1.5:19.5 Inf];

    confBounds = ([-Inf linspace(c1,c2,19) Inf]);
end
pnew = normcdf(confBounds(2:end)) - normcdf(confBounds(1:end-1));
pold = normcdf(confBounds(2:end),mu_old,sigma_old) - normcdf(confBounds(1:end-1),mu_old,sigma_old);

% making 0 probabilities very small (prevents LL from going to -Inf)
% picked this value because FP and VP model have 150*6*6 total simulations
pnew(pnew==0) = 1.8519e-04;
pold(pold==0) = 1.8519e-04;

% ===== PLOTS =====
nHist = 40;
if dplot == 1;
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


