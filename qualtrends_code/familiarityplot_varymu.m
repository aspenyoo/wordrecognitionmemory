
% "similarity" is manipulated by varying the mean of the distributions new
% words are drawn from
% greater distance from the true mean, less similar

modelname = 'FP';
theta = [20 1 4 0 nan 1];
islogbinning = 1; 
nX = 10; 
nS = 10;
nConf = 20;
N = [150 150];

% nSimilarities = 5;
% similarityVec = linspace(1e-3,2,nSimilarities); 
similarityVec = [0 0.5 1 1.5];
nSimilarities = length(similarityVec);

nIter = 1;
nNew = nan(nIter,nSimilarities,nConf);
nOld = nan(nIter,nSimilarities,nConf);
for isim = 1:nSimilarities;
    theta(5) = similarityVec(isim);
    
    for iiter = 1:nIter;
    [nNew(iiter,isim,:), nOld(iiter,isim,:)] = simulatedata_varySnew(modelname, theta, islogbinning, nX, nS, nConf, N);
    end
    
end

nNew = squeeze(mean(nNew,1))./N(1);     % averaging across iterations and normalizing
nOld = squeeze(mean(nOld,1))./N(2);     % averaging across iterations and normalizing

meanConfFA = sum(bsxfun(@times,11:20,bsxfun(@rdivide,nNew(:,11:20),sum(nNew(:,11:20),2))),2)

% figure;
clf
simcolors = flipud(aspencolors(nSimilarities,'blue'));
for isim = 1:nSimilarities;
    plot(1:20,nNew(isim,:),'Color',simcolors(isim,:)); hold on;
end
hold off
defaultplot
legend('high',' ',' ', 'low')
title('response proportions with similarity')