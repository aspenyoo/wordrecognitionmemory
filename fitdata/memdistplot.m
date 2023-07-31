function memdistplot(dold, dnew, binningparameters, binningfn)
% plots memory distribution
%
% ========== INPUT VARIABLES ==========
% DNEW: vector of DNEW values. if DNEW is a matrix, then MEMDISTPLOT will
% plot the mean and SEM for each row.

if length(size(dold))>2; error('matrix dold too large'); end
nX = min(size(dold));
if nX ~= 1;
    summarystat = 1;
else
    summarystat = 0;
end

L = 10;
nHist = 100;

gold = aspencolors('dustygold');
greyblue = aspencolors('greyblue');
%     figure; subplot(2,1,1)
%     figure;
confidence = 0.5:19.5;
hold on;    

[~,centers] = hist([dold(:)', dnew(:)'],nHist);
[nelemold] = hist(dold,centers);
[nelemnew] = hist(dnew,centers);

% eliminating zeros and normalizing (normalizing for graphing)
nelemold = (1+nelemold)/sum(nelemold+1);
nelemnew = (1+nelemnew)/sum(nelemnew+1);
histcount = [nelemold(:)' nelemnew(:)'];

if(summarystat)
    mean_dold = mean(nelemold,2);
    mean_dnew = mean(nelemnew,2);
    sem_dold = std(nelemold,[],2)/sqrt(nX);
    sem_dnew = std(nelemnew,[],2)/sqrt(nX);
    
    % plotting
    plot(centers,mean_dold,'-','LineWidth',2,'Color',gold);
    plot(centers,mean_dnew,'-','LineWidth',2,'Color',greyblue);
    errfNew = fill([centers fliplr(centers)],[mean_dnew-sem_dnew; ...
        flipud(mean_dnew+sem_dnew) ]', greyblue);
    errfOld = fill([centers fliplr(centers)],[mean_dold-sem_dold; ...
        flipud(mean_dold+sem_dold) ]', gold);
    
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
else
    plot(centers,nelemold,'-','LineWidth',2,'Color',gold);
    plot(centers,nelemnew,'-','LineWidth',2,'Color',greyblue);
    set(gca,'YTick',[])
end
% 
% switch binningfn
%     case 3
%         a = binningparameters(1);
%         b = binningparameters(2);
%         d0 = binningparameters(3);
%         d = [-log(1./exp((confidence(1:length(confidence)/2)-b)./a)-1)-d0 log(1./exp((confidence(1:length(confidence)/2)-b)./a)-1)-d0];
%     case 2
%         a = binningparameters(1);
%         b = binningparameters(2);
%         d0 = binningparameters(3);
%         d = [exp((confidence(1:length(confidence)/2) - b)./a)-d0 -exp((confidence(1:length(confidence)/2) - b)./a)-d0];
%     case 1
%         k = binningparameters(1);
%         d0 = binningparameters(2);
%         d = -k.*log(((2*L)./(confidence - 10.5 + L)) -1) + d0;
%     case 0
%         c1 = binningparameters(1);
%         c2 = binningparameters(2);
%         d = [-Inf linspace(c1,c2,19) Inf];
% end
% 
% for j = 1:length(confidence);
%     maxx = max(histcount);
%     hLine = plot([d(j) d(j)], [0 maxx],'Color',[.85 .85 .85]);
%     set(get(get(hLine,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','off'); % Exclude line from legend
% end

ylim([0 max(histcount)])
defaultplot

% set(  gca                         ,...
%     'YTick'         , []        ,...
%     'YColor'        ,'w'        );
%     hold off