% estimate relationship between M and sigma in word
% recognition memory model by looking at the best fit parameter for each
% fixed value of M and plot that for each subject. (and average)

clear all
modelname = 'FP3';
MVec = [1:65 70:5:90];
nMs = length(MVec);
Mcol = 1;       % column with the M listed
sigmacol = 2;   % column with the sigmas listed
nLLcol = 7;     % column with the nLLs listed
nSubj = 14;
subplotsize = 4;

mu = sqrt(2)*(gamma((MVec + 1)/2))./(gamma(MVec/2));
sigmaVec = nan(nSubj,nMs);
nLLMat = nan(nSubj,nMs);
for isubj = 1:nSubj;
    isubj
    % get best fit parameter for each M
    filename = ['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt'];
    alldata = dlmread(filename);

    % plot data
    subplot(subplotsize,subplotsize,isubj)
    colormap('parula')
    scatter(alldata(:,1),alldata(:,2),[],-alldata(:,size(alldata,2)/2),'filled'); defaultplot
    colorbar;
    title(sprintf('Subject %d',isubj))
    if mod(isubj,4) == 1; ylabel('\sigma'); end
    if isubj > (subplotsize*subplotsize-subplotsize); xlabel('M'); end
    
end