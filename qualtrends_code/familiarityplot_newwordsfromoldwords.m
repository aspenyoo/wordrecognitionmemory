% each new word is drawn from a distribution centered around an old word
% "similarity" is manipulated by varying the variance of the distribution
% that these new words are drawn from. 
% larger variance: less similar. 
% 
% plots the mean confidence rating for increasing levels of dissimilarity
% (i.e., increasing variance of the distribution new words are drawn from
% around old words)

modelname = 'FP';
theta = [20 1 4 0 0 nan];
islogbinning = 0; 
nX = 20; 
nS = 20;
nConf = 20;
N = [150 150];

nSimilarities = 50;
similarityVec = linspace(0,2,nSimilarities);%[0 0.5 1 1.5];
% nSimilarities = length(similarityVec);


pnew = nan(nSimilarities,nConf);
pold = nan(nSimilarities,nConf);
for isim = 1:nSimilarities;
    isim
    theta(end) = similarityVec(isim);
    [pnew(isim,:), pold(isim,:)] = simulatedata_newwordsfromoldwords(modelname, theta, islogbinning, nX, nS, nConf, N);
end

normpnewFA = bsxfun(@rdivide, pnew(:,11:20),sum(pnew(:,11:20),2));
meanConfFA = normpnewFA*[11:20]';

% figure;
% clf
% simcolors = flipud(aspencolors(nSimilarities,'blue'));
% % simcolors = bsxfun(@times,linspace(0,0.8,nSimilarities)',ones(1,3));
% for isim = 1:nSimilarities;
%     plot(1:20,nNew(isim,:),'Color',simcolors(isim,:),'LineWidth',2); hold on;
% end
% hold off
% defaultplot
% legend(cellfun(@(x)['\sigma_{new} = ' num2str(x)],num2cell(similarityVec),'UniformOutput',false));
% % title('response proportions with similarity')
% xlabel('confidence reports')
% ylabel('probability density')
% axis([1 20 0 0.5])
% ax = gca;
% ax.XTick = [ 1 10 20];
% ax.YTick = [0 0.5];

figure;
meanConf = pnew*[1:20]';
plot(similarityVec,meanConf,'Color',0.7*ones(1,3),'LineWidth',2)
defaultplot;
xlabel('\sigma_{new}')
ylabel('mean confidence response')
axis([0 2 1 15])
ax = gca;
ax.XTick = [0 1 2];
ax.YTick = [1 5 10 15];
