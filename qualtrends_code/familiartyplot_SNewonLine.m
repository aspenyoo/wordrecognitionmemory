
modelname = 'FP';
theta = [20 1 3 0];
islogbinning = 1; 
nX = 10; 
nS = 10;
nConf = 20;
N = [150 150];

nIter = 10;

% plot mean histograms as a function of similarity
nSimilarities = 3;
dissimilarityMax = 3;
similarities = linspace(0,dissimilarityMax,nSimilarities);
% for iiter = 1:nIter;
    [nnew, nold] = simulate_data(modelname,theta,islogbinning,nX,nS,nConf,N);

% end

% plot d(memory strength) as a function of dissimilarity (larger number,
% less similar)
figure
dMeanMat = nan(nIter,N(2));
for iiter = 1:nIter;
[d_mean, lengths] = simulatedata_snewonline(modelname,theta,islogbinning,nX,nS,nConf,N);
hold on
dMeanMat(iiter,:) = d_mean;

end
hold off

figure; 
meann = mean(dMeanMat);
stdd = std(dMeanMat)/sqrt(nIter);

errorbar(lengths,meann,stdd);
defaultplot