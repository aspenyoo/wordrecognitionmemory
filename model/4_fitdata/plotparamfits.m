function plotparamfits(modelname, binningfn, memstrengthvar, optimizationMethod, bestFitParam, nConf, isdatasave, fakedata, subjnum, selectiveplot)
% function that plots confdist (data and model overlaid)

nSubj = size(bestFitParam,1);
nParams = size(bestFitParam,2);
if nargin < 6; nConf = 20; end
if nargin < 7; isdatasave = 0; end
if nargin < 8; fakedata = 0; end % default: real subject data
if nargin < 9; subjnum = 1:nSubj; end
if isempty(subjnum); subjnum = 1:nSubj; end
if nargin < 9; selectiveplot = [1 1 1 0]; end % automaticlyy plots everything
date = clock;

subplotsize = ceil(sqrt(nSubj));

if isstruct(fakedata)
    pNew_part = fakedata.new;
    pOld_part = fakedata.old;
else
    % load real data
    load('subjdata.mat')
    pNew_part = nNew_part/150;
    pOld_part = nOld_part/150;
end

% ====================================================================
% PLOTS! (indvl, ave (average), d (memDist))
% ====================================================================

gold = aspencolors('dustygold');
greyblue = aspencolors('greyblue');

% switch modelname
%     case {'FP','FPheurs','uneqVar'}
%         nParams = 4;
%     case 'REM'
%         nParams = 7;
% end
% if any(binningfn == [2 3 4]); nParams = nParams + 2; end 
% if binningfn == 5; nParams = nParams + 3; end

% getting simulated data
pNew_est = nan(nSubj,20); pOld_est = pNew_est;
nX = 30; nS = 50;
for isubjnum = 1:nSubj;
    
    subjnum(isubjnum)
    [pNew_est(isubjnum,:), pOld_est(isubjnum,:)] = nLL_approx_vectorized(modelname, bestFitParam(isubjnum,:), binningfn, memstrengthvar, nNew_part(isubjnum,:), nOld_part(isubjnum,:), [], nX, nS, nConf );
end

% % change numbers to probabilities
% N = 150;
% pNew_est = pNew_est/N;
% pOld_est = pOld_est/N;

if nConf ~= 20 && (size(pNew_est,2)==20);
    pNew_esttemp = nan(nSubj,nConf);
    pOld_esttemp = pNew_esttemp;
    num = round(20/nConf);
    for i = 1:nConf-1;
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
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    for isubjnum = 1:nSubj;
        
        isubj = subjnum(isubjnum);
        
        % plots
        subplot(subplotsize,subplotsize,isubj); hold on;
        fNew = plot(pNew_est(isubjnum,:)); % pnewMat(:,MInd,sigInd,c1Ind,c2Ind) are best fit thing
        fOld = plot(pOld_est(isubjnum,:)); % poldMat(:,MInd,sigInd,c1Ind,c2Ind)
        sNew = plot(pNew_part(isubjnum,:));
        sOld = plot(pOld_part(isubjnum,:));
        
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
        labeltext = ['participant ' num2str(isubj)];
        text(.39,.5,labeltext);
        
        
    end
    
    %big title
    h = axes('Position',[0 0 1 1],'Visible','off'); %add an axes on the left side of your subplots
    hTitle = text(0.45,0.95,[num2str(nParams) ' Parameter Model']);
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
    set(hTitle,...
        'FontSize'      , 18                ,...
        'FontWeight'    ,'bold'             );
    set(hYlabel,...
        'VerticalAlignment'     ,'bottom'   ,...
        'HorizontalAlignment'   ,'left'     ,...
        'Rotation'              , 90        ,...
        'FontSize'              , 14         );
    set(hXlabel,...
        'FontSize'              , 14        );
    
    %     subplotlegend(hLegend,[4 4 nSubj+1]);
    subplotlegend(hLegend,[subplotsize subplotsize subplotsize^2]);
    
    if (isdatasave)
        %         saveas(gcf,['paramfit_' modelname linlog '_' optimizationMethod ...
        %              '_confdist' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.emf'])
        img = getframe(gcf);
        imwrite(img.cdata, ['paramfit_' modelname num2str(binningfn) '_' optimizationMethod ...
            '_confdist_' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.png']);
    end
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
    figure('units','normalized','outerposition',[0 0 2/3 2/3]); hold on;
    
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
    title([modelname num2str(binningfn) ' Parameter Group Fit']);
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
    
    
    if (isdatasave)
        %         saveas(gcf,['paramfit_' modelname linlog '_' optimizationMethod ...
        %              '_confdistAverage' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.emf'])
        img = getframe(gcf);
        imwrite(img.cdata, ['paramfit_' modelname num2str(binningfn) '_' optimizationMethod ...
            '_confdistAverage_' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.png']);
    end
end

% ===== INDVL DPLOT =====

if (selectiveplot(3))
    figure('units','normalized','outerposition',[0 0 1 1]);
    for isubjnum = 1:nSubj;
        
        isubj = subjnum(isubjnum);
        subplot(subplotsize,subplotsize,isubj);
        
%         switch modelname
%             case 'FP'
                simulate_data(modelname,bestFitParam(isubjnum,:),binningfn, 30, 20,20,150,[0 1]);
%             case 'VP'
%                 responses_VP( bestFitParam(isubjnum,:),islogbinning, 0,1);
%             case 'FPheurs'
%                 responses_FPheurs( bestFitParam(isubjnum,:),islogbinning, 0,1);
%             case 'VPheurs'
%                 responses_VPheurs( bestFitParam(isubjnum,:),islogbinning, 0,1);
%             case 'uneqVar'
%                 responses_uneqVar( bestFitParam(isubjnum,:), islogbinning,0,1);
%         end
        
        title(['participant ' num2str(isubj)])
    end
    hold off;
    legend('old','new');
    hLeg = legend('old','new');
    %     subplotlegend(hLeg,[4 4 nSubj+1])
    subplotlegend(hLeg,[subplotsize subplotsize subplotsize^2])
    
    if (isdatasave)
        %         saveas(gcf,['paramfit_' modelname linlog '_' optimizationMethod ...
        %              '_memdist' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.emf'])
        img = getframe(gcf);
        imwrite(img.cdata, ['paramfit_' modelname num2str(binningfn) '_' optimizationMethod ...
            '_memdist_' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.png']);
    end
end


% ===== BOXPLOT AVERAGE =====

if (selectiveplot(4))
    
    % plots
    figure('units','normalized','outerposition',[0 0 2/3 2/3]); hold on;
    
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
    text(.39,.39,'average participant');
    title([num2str(nParams) ' Parameter Group Fit']);
    ylabel('Proportion');
    xlabel('Bin Number');
    
    hold off;
    
    
    if (isdatasave)
        img = getframe(gcf);
        imwrite(img.cdata, ['paramfit_' modelname num2str(binningfn) '_' optimizationMethod ...
            '_confdistBoxplot_' num2str(date(2)) num2str(date(3)) num2str(date(1)) '.png']);
    end
end
