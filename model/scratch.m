

for i = 1:5;
    i
    filename = 'blahhhh.txt';
    permission = 'a+'; % open or create new file for reading and writing. append data to the end of the file
    fileID = fopen(filename,permission);
    
    A1 = nan(5,4);
    formatSpec = '%2.0f \t %2.4f \t %2.4f \t %4.4f \r\n';
    fprintf(fileID, formatSpec, A1);
%     pause;
    fclose(fileID);
end

%% load the data in

nSubj = 25;
nM = 50;
blah = nan(nSubj,nM);
for isubj = 1:nSubj;
filename = ['paramfit_patternbayes_FPheurs_subj' num2str(isubj) '.txt'];
bleh = dlmread(filename);
blah(isubj,:) = histc(bleh(:,1),1:nM);
end

%% see how many params fit for each M for each subject

modelname = 'FPheurs';
nSubj = 14;
counts = nan(50,nSubj);
for isubj = 1:nSubj;
    removetxtspaces(modelname,isubj);
    blah = dlmread(['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt']);
    counts(:,isubj) = histc(blah(:,1),1:50);
end

%%
fileid = fopen(['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt']);
C = fread(fileid);


%% mean of chi distribution
MVec = 1:50;
mu = sqrt(2)*(gamma((MVec+1)/2)./gamma(MVec/2));
plot(MVec,mu)




%% MARCH 31, 2016

% heatmap of correlation of LPR and distanceNew (distance between new words
% and closest noisy memory)

modelname = 'FP';
MVec = 1:50;
sigmaVec = exp(linspace(.01,1.5,50));

rho_new = nan(length(MVec),length(sigmaVec));
pvalue_new = nan(length(MVec),length(sigmaVec));
for iM = 1:length(MVec);
    M = MVec(iM);
    
    for isigma = 1:length(sigmaVec);
        sigma = sigmaVec(isigma);
        
        [rho_new(iM,isigma), pvalue_new(iM,isigma)] = correlation_LPR_distanceNew(modelname, M, sigma);
        
    end
end

[Mmesh,sigmamesh] = meshgrid(MVec,sigmaVec);

% figure;
% imagesc(Mmesh(:),sigmamesh(:),rho_new(:))

figure;
Mmesh = Mmesh(:); sigmamesh = sigmamesh(:); rho = rho_new'; rho = rho(:);
Mmesh = Mmesh(~isnan(rho)); sigmamesh = sigmamesh(~isnan(rho));
rho = rho(~isnan(rho));
scatter(Mmesh(:),log(sigmamesh(:)),[],rho(:),'filled');
defaultplot
xlabel('M'); ylabel('log(\sigma)'); title('Spearman \rho')
% HeatMap(rho_new)

blah = [.3 .8 .8; 0.5*ones(1,3); .8 .4 .3];
bleh = [linspace(blah(1,1),blah(2,1),35) linspace(blah(2,1),blah(3,1),20);...
linspace(blah(1,2),blah(2,2),35) linspace(blah(2,2),blah(3,2),20);...
linspace(blah(1,3),blah(2,3),35) linspace(blah(2,3),blah(3,3),20)];
colormap(bleh')

%% APRIL 1, 2016

% heatplot of correlation of mean confidence as a function of how far new words were drawn
% from the old word, for different M and sigmas

MVec = 1:10;
sigmaVec = linspace(.01,5,50);
newsigmaVec = 0:0.5:2;

modelname = 'FP';
theta = [nan nan 1 0 nan];
islogbinning = 1;
nX = 1; 
nS = 1;

r = nan(length(MVec), length(sigmaVec));
pvalue = nan(length(MVec), length(sigmaVec));
for iM = 1:length(MVec);
    iM
    theta(1) = MVec(iM);
    
    for isigma = 1:length(sigmaVec);
        theta(2) = sigmaVec(isigma);
        
        mean_nNew = nan(1,length(newsigmaVec));
        for inewsigma = 1:length(newsigmaVec);
            theta(end) = newsigmaVec(inewsigma);
            [nNew,nOld,dNew, dOld] = simulatedata_newwordsfromoldwords(modelname, theta, islogbinning, nX, nS);
        
            mean_nNew(inewsigma) = wmean(1:20,nNew/sum(nNew));
%             mean_dNew(inewsigma) = nanmean(dNew);
        end
        [tempr,tempp] = corr(newsigmaVec',mean_nNew','type','Spearman');
        r(iM,isigma) = tempr;
        pvalue(iM,isigma) = tempp;
    end
end

% plot stuff
figure;
imagesc(MVec,sigmaVec,r')
set(gca,'YDir','normal');
xlabel('M'); ylabel('\sigma')
title('correlation with \sigma_{new} and mean confidence rating')
defaultplot

% get colormap
posN = 20; negN = posN;
blah = [aspencolors('seablue'); 0.5*ones(1,3); aspencolors('coral')];
% blah = [.3 .8 .8; 0.5*ones(1,3); .8 .4 .3];
bleh = [linspace(blah(1,1),blah(2,1),posN) linspace(blah(2,1),blah(3,1),negN);...
linspace(blah(1,2),blah(2,2),posN) linspace(blah(2,2),blah(3,2),negN);...
linspace(blah(1,3),blah(2,3),posN) linspace(blah(2,3),blah(3,3),negN)];
colormap(bleh')