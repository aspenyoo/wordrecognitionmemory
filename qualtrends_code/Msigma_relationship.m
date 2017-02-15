% estimate relationship between M and sigma in word
% recognition memory model by looking at the best fit parameter for each
% fixed value of M and plot that for each subject. (and average)

clear all
modelname = 'FP3';
MMax = 65;
MMin = 1;
MVec = MMin:MMax;
Mcol = 1;       % column with the M listed
sigmacol = 2;   % column with the sigmas listed
nLLcol = 7;     % column with the nLLs listed
nSubj = 14;
subplotsize = 4;

mu = sqrt(2)*(gamma((MVec + 1)/2))./(gamma(MVec/2));
sigmaVec = nan(nSubj,MMax-MMin+1);
nLLMat = nan(nSubj,MMax-MMin+1);
for isubj = 1:nSubj;
    isubj
    % get best fit parameter for each M
    filename = ['paramfit_patternbayes_' modelname '_subj' num2str(isubj) '.txt'];
    alldata = dlmread(filename);
    
    clear bestdata
    for iM = MMin:MMax;
        try
            [bestdata(iM,:), nLLMat(isubj,iM)] = getbestfitparams(modelname(1:end-1),str2double(modelname(end)),isubj,[],[1; iM; iM]);
        catch
            bestdata(iM,:) = nan(1,size(bestdata,2));
            nLLMat(isubj,iM) = nan;
        end
    end
    sigmaVec(isubj,:) = bestdata(MMin:MMax,sigmacol)';
    
    % find bestfitline
    E_chi = @(x) x*sqrt(2)*(gamma((MVec+1)/2)./gamma(MVec/2));
    costfunc = @(x) abs(sum(sigmaVec(isubj,:) - E_chi(x)));
    [sigmatilde(isubj)] = fminsearch(costfunc,1);
    
    % plot data
    subplot(subplotsize,subplotsize,isubj)
    colormap('parula')
    scatter(MVec,sigmaVec(isubj,:),[],-nLLMat(isubj,MMin:MMax),'filled'); defaultplot
    colorbar;
    title(sprintf('Subject %d',isubj))
    if mod(isubj,4) == 1; ylabel('\sigma'); end
    if isubj > (subplotsize*subplotsize-subplotsize); xlabel('M'); end
    
    % plot best fit line
    hold on;
    bestFitLine(isubj,:) = E_chi(sigmatilde(isubj));
    plot(MVec,bestFitLine(isubj,:),'Color',0.6*ones(1,3));
end
biglabelplot('Estimated \sigma based on M')

mean_sigma = mean(sigmaVec);
sem_sigma = std(sigmaVec)/sqrt(nSubj);
mean_model = mean(bestFitLine);
sem_model = std(bestFitLine)/sqrt(nSubj);
figure;
plot_summaryfit(MVec,mean_sigma,sem_sigma,mean_model,sem_model,0.5*ones(1,3),aspencolors('greyblue'));
% plot(MVec,mu,'LineWidth',2); hold on
% herr = errorbar(MMin:MMax,mean_sigma,sem_sigma);
% set(herr,'Color',0.7*ones(1,3));
xlabel('M')
ylabel('\sigma')
defaultplot;
title('average \sigma across M')