function plotparamfits(modelname, bestFitParam, binningfn, nConf, fakedata, subjnum, selectiveplot)
% function that plots confdist (data and model overlaid)

nSubj = size(bestFitParam,1);
if nargin < 4; nConf = 20; end
if nargin < 5; fakedata = 0; end % default: real subject data
if nargin < 6; subjnum = 1:nSubj; end
if isempty(subjnum); subjnum = 1:nSubj; end
if nargin < 7; selectiveplot = [1 1 1 0]; end % automaticlyy plots everything


subplotsize = ceil(sqrt(nSubj));

if isstruct(fakedata)
    pNew_part = fakedata.new;
    pOld_part = fakedata.old;
else
    % load real data
    load('subjdata.mat')
    pNew_part = bsxfun(@rdivide,nNew_part,sum(nNew_part,2));
    pOld_part = bsxfun(@rdivide,nOld_part,sum(nOld_part,2));
end

% ====================================================================
% PLOTS! (indvl, ave (average), d (memDist))
% ====================================================================

gold = aspencolors('dustygold');
greyblue = aspencolors('greyblue');

switch modelname
    case 'UVSD'
        nParams = 2;
    case {'FP','FPheurs'}
        nParams = 2;
    case 'REM'
        nParams = 5;
end

switch binningfn
    case {0,1} % linear, logistic
        % slope, y-int, sigma_mc
        nParams = nParams + 3;
    case 2     % logarithmic
        % a, b, d0, sigma_mc
        nParams = nParams + 4;
    case 3      % power law
        % a, b, gamma, d0, sigma_mc
        nParams = nParams + 5;
    case 4 % weibull
        % a, b, shape, scale, d0, sigma_mc
        nParams = nParams + 6;
end

% getting simulated data
if any(selectiveplot(1:2))
    pNew_est = nan(nSubj,20); pOld_est = pNew_est;
    for isubjnum = 1:nSubj
        
        subjnum(isubjnum)
        [pNew_est(isubjnum,:), pOld_est(isubjnum,:)] = nLL_approx_vectorized(modelname, bestFitParam(isubjnum,:), binningfn, nNew_part(isubjnum,:), nOld_part(isubjnum,:));
    end
end

% % change numbers to probabilities
% N = 150;
% pNew_est = pNew_est/N;
% pOld_est = pOld_est/N;

if nConf ~= 20 && (size(pNew_est,2)==20)
    pNew_esttemp = nan(nSubj,nConf);
    pOld_esttemp = pNew_esttemp;
    num = round(20/nConf);
    for i = 1:nConf-1
        pNew_esttemp(:,i) = sum(pNew_est(:,num*(i-1)+1:num*i),2);
        pOld_esttemp(:,i) = sum(pOld_est(:,num*(i-1)+1:num*i),2);
    end
    pNew_esttemp(:,nConf) = sum(pNew_est(:,(num*(nConf-1)+1):end),2);
    pOld_esttemp(:,nConf) = sum(pOld_est(:,(num*(nConf-1)+1):end),2);
    pNew_est = pNew_esttemp;
    pOld_est = pOld_esttemp;
end


% ===== INDVL PLOT =====
if (selectiveplot(1))
    if sum(selectiveplot)>1; figure('units','normalized','outerposition',[0 0 1 1]);end
    
    for isubjnum = 1:nSubj;
        
        isubj = subjnum(isubjnum);
        
        % plots
        if nSubj > 1;  subplot(subplotsize,subplotsize,isubjnum); end
        hold on;
        fNew = plot(pNew_est(isubjnum,:)); % pnewMat(:,MInd,sigInd,c1Ind,c2Ind) are best fit thing
        fOld = plot(pOld_est(isubjnum,:)); % poldMat(:,MInd,sigInd,c1Ind,c2Ind)
        sNew = plot(pNew_part(isubj,:));
        sOld = plot(pOld_part(isubj,:));
        
        defaultplot;
        % aesthetics
        axis([0 nConf+1 0 .55]);
        set(fNew                                    ,...
            'Color'             , greyblue          ,...
            'LineWidth'         , 2                 );
        set(fOld                                    ,...
            'Color'             , gold               ,...
            'LineWidth'         , 2                 );
        set(sNew                            ,...
            'Marker'            ,'.'        ,...
            'MarkerFaceColor'   , greyblue        ,...
            'MarkerSize'        ,14        ,...
            'Color'             , greyblue        ,...
            'LineWidth'         , 2         ,...
            'LineStyle'         ,'none'     );
        set(sOld                            ,...
            'Marker'            ,'.'        ,...
            'MarkerFaceColor'   , gold       ,...
            'MarkerSize'        , 14        ,...
            'Color'             , gold       ,...
            'LineWidth'         , 2         ,...
            'LineStyle'         ,'none'     );
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.02 .02] , ...
            'XTick'       , 0:5:nConf  , ...
            'YMinorTick'  , 'off'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 1         );
        %  'YTick'       , 0:.25:.5    ,...%                 'YTick'       , 0:.25:.5    ,...
        % labels
%         labeltext = ['participant ' num2str(isubj)];
%         text(.39,.5,labeltext);
        
        
    end
    
    %big title
    h = axes('Position',[0 0 1 1],'Visible','off'); %add an axes on the left side of your subplots
%     hTitle = text(0.45,0.95,[modelname]);
    %big ylabel
    set(gcf,'CurrentAxes',h)
    hYlabel = text(.1,.45,'Proportion');
    % big xlabel
    set(gcf,'CurrentAxes',h)
    hXlabel = text(.45,0.05,'Bin Number');
    % legend
    hLegend  = legend(...
        [sNew, sOld, fNew, fOld],...
        'New Words (participant)'       ,...
        'Old Words (participant)'       ,...
        'New Words (estimated)'         ,...
        'Old Words (estimated)'         );
    
    
    % label aesthetics
%     set(hTitle,...
%         'FontSize'      , 18                ,...
%         'FontWeight'    ,'bold'             );
    set(hYlabel,...
        'VerticalAlignment'     ,'bottom'   ,...
        'HorizontalAlignment'   ,'left'     ,...
        'Rotation'              , 90        ,...
        'FontSize'              , 14         );
    set(hXlabel,...
        'FontSize'              , 14        );
    
    %     subplotlegend(hLegend,[4 4 nSubj+1]);
    if nSubj < subplotsize^2; subplotlegend(hLegend,[subplotsize subplotsize subplotsize^2]); end
    
end

% ===== AVERAGE CONFIDENCE PLOT =====

if (selectiveplot(2))
    % averaged participant (summary)
    pNew_partSum = mean(pNew_part);
    pOld_partSum = mean(pOld_part);
    sem_pNew_part = std(pNew_part)/sqrt(nSubj);
    sem_pOld_part = std(pOld_part)/sqrt(nSubj);
    pNew_estSum = mean(pNew_est);
    pOld_estSum = mean(pOld_est);
    sem_pNew_est = std(pNew_est)/sqrt(nSubj);
    sem_pOld_est = std(pOld_est)/sqrt(nSubj);
    
    % plots
    if sum(selectiveplot)>1; figure('units','normalized','outerposition',[0 0 2/3 2/3]); end
    
    hold on;
    errfNew = fill([1:nConf, nConf:-1:1],[pNew_estSum-sem_pNew_est ...
        fliplr(pNew_estSum+sem_pNew_est) ], greyblue);
    errfOld = fill([1:nConf, nConf:-1:1],[pOld_estSum-sem_pOld_est ...
        fliplr(pOld_estSum+sem_pOld_est) ], gold);
    errsNew = errorbar(pNew_partSum, sem_pNew_part);
    errsOld = errorbar(pOld_partSum, sem_pOld_part);
    
    defaultplot;
    
    % aesthetics
    axis([0 nConf+1 0 .2])
    set(errsNew                     ,...
        'Color'         , greyblue  ,...
        'LineStyle'     ,'none'     ,...
        'LineWidth'     , 1.5       );
    set(errsOld                     ,...
        'Color'         ,gold        ,...
        'LineStyle'     ,'none'     ,...
        'LineWidth'     , 1.5       );
    set(errfNew                     ,...
        'FaceColor'     , greyblue        ,...
        'LineStyle'     ,'none'     ,...
        'LineWidth'     , 1.5         ,...
        'FaceAlpha'     , 0.3       );
    set(errfOld                     ,...
        'FaceColor'     , gold        ,...
        'LineStyle'     ,'none'     ,...
        'LineWidth'     , 1.5         ,...
        'FaceAlpha'     , 0.3       );
    set(gca, ...
        'Box'         , 'off'    , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XTick'       , 0:5:nConf    , ...
        'YMinorTick'  , 'on'      , ...
        'YTick'       , [0 0.1 0.2] ,...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1         ,...
        'FontName'    ,'Helvetica' );
    
    % LABELS
    
    % label
    text(.39,.39,'average participant');
    title([modelname]);
    ylabel('Proportion');
    xlabel('Bin Number');
    
    hold off;
    
    % legend
    % shows up behind figure...too cluttered anyways
    %     hLegend  = legend(...
    %         [sNew, sOld, errfNew, errfOld],...
    %         'New Words (participant)'       ,...
    %         'Old Words (participant)'       ,...
    %         'New Words (estimated)'         ,...
    %         'Old Words (estimated)'         );
    
   
end

% ===== INDVL DPLOT =====

if (selectiveplot(3))
    if sum(selectiveplot)>1; figure('units','normalized','outerposition',[0 0 1 1]);end
    for isubjnum = 1:nSubj;
        
        isubj = subjnum(isubjnum);
        subplot(subplotsize,subplotsize,isubjnum);
        
        [centers_new, counts_new, centers_old, counts_old, confBounds] = nLL_approx_vectorized( modelname, bestFitParam(isubjnum,:), binningfn, nNew_part(isubj,:), nOld_part(isubj,:));
        
        % plots
        subplot(subplotsize,subplotsize,isubjnum); hold on;
        fNew = plot(centers_new,counts_new);
        fOld = plot(centers_old,counts_old);
        fConf = plot(repmat(confBounds(:),1,2)', repmat([0 max([counts_new(:); counts_old(:)])],length(confBounds),1)');
        
        % aesthetics
        defaultplot;
        set(fNew                                    ,...
            'Color'             , greyblue          ,...
            'LineWidth'         , 2                 );
        set(fOld                                    ,...
            'Color'             , gold               ,...
            'LineWidth'         , 2                 );
        set(fConf                                   ,...
            'Color'             ,0.7*ones(1,3)      );
        title(['participant ' num2str(isubj)])
    end
    
    
    hold off;
    legend('old','new');
    hLeg = legend('old','new');
    %     subplotlegend(hLeg,[4 4 nSubj+1])
    if nSubj > subplotsize^2; subplotlegend(hLeg,[subplotsize subplotsize subplotsize^2]); end
    
end


% ===== BOXPLOT AVERAGE =====

if (selectiveplot(4))
    
    % plots
    if sum(selectiveplot)>1; figure('units','normalized','outerposition',[0 0 2/3 2/3]); end
    
    hold on;
    modelNew = boxplot(pNew_est       ,...
        'plotstyle'    , 'compact'  ,...
        'boxstyle'     ,'outline'   ,...
        'symbol'       , 'o'        ,...
        'colors'       , greyblue+0.2   );
    modelOld = boxplot(pOld_est       ,...
        'plotstyle'    , 'compact'  ,...
        'boxstyle'     ,'outline'   ,...
        'symbol'       , 'o'        ,...
        'colors'       , gold+0.2   );
    partNew = boxplot(pNew_part       ,...
        'plotstyle'    , 'compact'  ,...
        'symbol'       , '+'        ,...
        'colors'       , greyblue   );
    partOld = boxplot(pOld_part       ,...
        'plotstyle'    , 'compact'  ,...
        'symbol'       , '+'        ,...
        'colors'       , gold   );
%     errsNew = errorbar(pNew_partSum, sem_pNew_part);
%     errsOld = errorbar(pOld_partSum, sem_pOld_part);
    
    defaultplot;
    
    % aesthetics
    axis([0 nConf+1 0 0.6])
    
%     set(errsNew                     ,...
%         'Color'         , greyblue  ,...
%         'LineStyle'     ,'none'     ,...
%         'LineWidth'     , 1.5       );
%     set(errsOld                     ,...
%         'Color'         ,tan        ,...
%         'LineStyle'     ,'none'     ,...
%         'LineWidth'     , 1.5       );
    set(gca, ...
        'Box'         , 'off'    , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XTick'       , 0:5:nConf    , ...
        'YMinorTick'  , 'on'      , ...
        'YTick'       , [0 0.1 0.2] ,...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1         ,...
        'FontName'    ,'Helvetica' );
    
    % LABELS
    
    % label
%     text(.39,.39,'average participant');
    title([modelname]);
    ylabel('Proportion');
    xlabel('Bin Number');
    
    hold off;

end
