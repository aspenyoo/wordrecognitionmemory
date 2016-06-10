
% "similarity" is manipulated by varying the SD of the distributions new
% words are drawn from
% a little confusing to think about how smaller variance would affect
% decisions. is it more similar or less similar than the distribution S
% words are drawn from?



modelname = 'FP';
theta = [20 1 3 0 0 nan];
islogbinning = 1; 
nX = 10; 
nS = 10;
nConf = 20;
N = [150 150];

% nSimilarities = 5;
% similarityVec = linspace(1e-3,2,nSimilarities); 
similarityVec = [0.01 0.5 1 1.5 2]; % 1 is the same sigma as S distribution
nSimilarities = length(similarityVec);

nIter = 1;
nNew = nan(nIter,nSimilarities,nConf);
nOld = nan(nIter,nSimilarities,nConf);
for isim = 1:nSimilarities;
    theta(6) = similarityVec(isim);
    
    for iiter = 1:nIter;
    [nNew(iiter,isim,:), nOld(iiter,isim,:)] = simulatedata_varySnew(modelname, theta, islogbinning, nX, nS, nConf, N);
    end
    
end

nNew = squeeze(mean(nNew,1))./N(1);     % averaging across iterations and normalizing
nOld = squeeze(mean(nOld,1))./N(2);     % averaging across iterations and normalizing

figure;
simcolors = flipud(aspencolors(nSimilarities,'blue'));
for isim = 1:nSimilarities;
    plot(1:20,nNew(isim,:),'Color',simcolors(isim,:)); hold on;
end
defaultplot
legend('lower variance',' ','same variance',' ', 'higher variance')
title('response proportions with similarity')