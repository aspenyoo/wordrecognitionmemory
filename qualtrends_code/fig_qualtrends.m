function fig_qualtrends(isubj)

%% % % % % % % % % % % % % % % % % % % % % 
% plot qualitative figure plot
% % % % % % % % % % % % % % % % % % % % % 

% get theta
modelname = 'FP';
binningfn = 4;
load(['paramfit_patternbayes_' modelname num2str(binningfn) '.mat'],'bestFitParam')
theta = bestFitParam(isubj,:);

% subplot dimensions
dim1 = 2; 
dim2 = 4; 

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

%% list length effect

listlengthVec = [20 40 100 150 300];
nListlengths = length(listlengthVec);

[pnewVec, poldVec] = deal(nan(nListlengths,20));
dprimeVec = nan(1,nListlengths);
for ilistlength = 1:nListlengths
    ilistlength
    listlength = listlengthVec(ilistlength);
    [nnew_part, nold_part] = deal([listlength zeros(1,19)]);
    
    [ pnewVec(ilistlength,:), poldVec(ilistlength,:) ] = nLL_approx_vectorized( modelname, theta, binningfn, nnew_part, nold_part );
    mean_pnew = pnewVec(ilistlength,:)*[1:20]';
    dprimeVec(ilistlength) = (poldVec(ilistlength,:)*[1:20]' - mean_pnew)/sqrt(listlength/(listlength-1).*pnewVec(ilistlength,:)*(([1:20]'- mean_pnew).^2));
end

% plot
if 0
colorss = [0.7*ones(1,3); 0.4*ones(1,3); zeros(1,3)];
figure; hold on
for ilistlength = 1:nListlengths
    plot(1:20,pnew(ilistlength,:),1:20,pold(ilistlength,:),'Color',colorss(ilistlength,:));
end
defaultplot
end

if 1
    subplot(dim1,dim2,2);
    plot(listlengthVec,dprimeVec,'k.-','MarkerSize',16);
    axis([0 300 0 3])
    set(gca,'YTick',0:3,'XTick',0:50:300,'XTickLabel',{'0',' ','100',' ','200',' ','300'})
    xlabel('list length')
    ylabel('dprime')
    defaultplot
end

%% ROC of one subject

pnew = pnewVec(listlengthVec == 150,:);
pold = poldVec(listlengthVec == 150,:);
cumHit = cumsum(fliplr(pold));
cumFA = cumsum(fliplr(pnew));

subplot(dim1,dim2,1)
plot(cumFA,cumHit,'k.-','MarkerSize',16)
defaultplot
axis([0 1 0 1])
xlabel('FA')
ylabel('Hit')


%% plot mirror and variance effect

bigM = round(theta(1)*1.2);
smallSigma = theta(2)*0.85;

[nnew_part, nold_part] = deal([150 zeros(1,19)]);
[ pnewM, poldM ] = nLL_approx_vectorized( modelname, [bigM theta(2:end)], binningfn, nnew_part, nold_part );
[ pnewSigma, poldSigma ] = nLL_approx_vectorized( modelname, [theta(1) smallSigma theta(3:end)], binningfn, nnew_part, nold_part );

colorss = [0.7*ones(1,3); zeros(1,3)];
% newsimcolors = aspencolors(nCond+1,'blue');
% oldsimcolors = [0.9 0.9 0.7; aspencolors('dustygold')];

% increasing M plot
subplot(dim1,dim2,3)
plot(1:20,pnew,1:20,pold,'Color',colorss(1,:)); hold on;
plot(1:20,pnewM,1:20,poldM,'Color',colorss(2,:)); hold off
hold off
defaultplot
xlabel('confidence report')
ylabel('probability density')
ax = gca;
axis([1 20 0 0.2])
ax.XTick = [1 10 20];
ax.YTick = [0 0.1 0.2];

% decreasing sigma plot
subplot(dim1,dim2,7)
plot(1:20,pnew,1:20,pold,'Color',colorss(1,:)); hold on;
plot(1:20,pnewSigma,1:20,poldSigma,'Color',colorss(2,:)); hold off
hold off
defaultplot
xlabel('confidence report')
ylabel('probability density')
ax = gca;
axis([1 20 0 0.2])
ax.XTick = [1 10 20];
ax.YTick = [0 0.1 0.2];

%% plot zROC plots

zFA = norminv(cumFA);
zHit = norminv(cumHit);
zFAM = norminv(cumsum(fliplr(pnewM)));
zHitM = norminv(cumsum(fliplr(poldM)));
zFASigma = norminv(cumsum(fliplr(pnewSigma)));
zHitSigma = norminv(cumsum(fliplr(poldSigma)));

% increasing M plot
subplot(dim1,dim2,4)
plot(zFA,zHit,'Color',colorss(1,:)); hold on;
plot(zFAM,zHitM,'Color',colorss(2,:)); hold off;
defaultplot;
axis([-4 4 -4 4])
set(gca,'XTick',-4:2:4,'YTick',-4:2:4)
xlabel('zFA')
ylabel('zHit')

% decreasing \sigma plot
subplot(dim1,dim2,8)
plot(zFA,zHit,'Color',colorss(1,:)); hold on;
plot(zFASigma,zHitSigma,'Color',colorss(2,:)); hold off;
defaultplot;
axis([-4 4 -4 4])
set(gca,'XTick',-4:2:4,'YTick',-4:2:4)
xlabel('zFA')
ylabel('zHit')

%% similarity plot: correlation

operationalization = 'mean';
[confidence, distancee, rho, pvalue, b, statss] = qualtrend_correlation(modelname, theta, binningfn, operationalization);

subplot(dim1,dim2,5)
grey = 0.7*ones(1,3);
plot(distancee,confidence,'.','Color',aspencolors('greyblue'),'MarkerSize',12); hold on
plot([min(distancee) max(distancee)],[1 min(distancee); 1 max(distancee)]*b,'Color',grey); hold off
ylabel('confidence rating')
xlabel('dissimilarity')
if pvalue < 0.001;
    signn = '<';
    pvalue = 0.001;
else
    signn = '=';
end
% tt = annotation('textbox',[.6 .8 .3 .1],...
%     'String', ...
    sprintf('\\rho = %2.3f, p %c %2.3f, \nR^2 = %2.3f', rho, ...
    signn, pvalue, statss(1))
% set(tt,'LineStyle','none')
defaultplot
ax = gca;
ax.YTick = [ 1 10 20];

%% similarity plot: changing \sigma_{new}

thetaa = [theta nan]; % last parameter is \sigma_{new}
nSimilarities = 15;
similarityVec = linspace(0,2,nSimilarities);

[pnew, pold] = deal(nan(nSimilarities,20));
for isim = 1:nSimilarities;
    thetaa(end) = similarityVec(isim);
    [pnew(isim,:), pold(isim,:)] = simulatedata_newwordsfromoldwords(modelname, thetaa, binningfn);
end

normpnewFA = bsxfun(@rdivide, pnew(:,11:20),sum(pnew(:,11:20),2));
meanConfFA = normpnewFA*[11:20]';

subplot(dim1,dim2,6);
meanConf = pnew*[1:20]';
plot(similarityVec,meanConf,'k-','LineWidth',1)
defaultplot;
xlabel('\sigma_{new}')
ylabel('mean confidence')
axis([0 2 1 15])
ax = gca;
ax.XTick = [0 1 2];
ax.YTick = [1 5 10 15];