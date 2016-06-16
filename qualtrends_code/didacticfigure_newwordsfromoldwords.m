function didacticfigure_newwordsfromoldwords

M = 2;
sigma = 0.1;
N = 10;
signewVec = [0.05 0.5];

figure
for isig = 1:2;
sigma_new = signewVec(isig);

% getting new and old words
X = randn(N,M)*sqrt(sigma^2+1);
SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(N,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X
SNew = SOld + randn(N,M).*sigma_new;

% idx of minimum distance
distances = squeeze(sum((permute(repmat(SNew,[1,1,N]), [3,2,1])-repmat(X,[1,1,N])).^2,2));
mindist = squeeze(min(distances,[],1));
[row,~] = find(distances == repmat(mindist,[N 1]));

% plot the new and old words with the distance lines between them
subplot(1,2,isig); hold on
plot(SNew(:,1),SNew(:,2),'k.','MarkerSize',12); % plot new words
plot(SOld(:,1),SOld(:,2),'ko'); % plot old words
plot([SNew(:,1) SOld(row,1)]',[SNew(:,2) SOld(row,2)]','--k')
axis equal
axismax = 3;
axis([-axismax axismax -axismax axismax])
set(gca,'XTick',[-3 0 3],'YTick',[-3 0 3]);
% ax = gca;
% ax.XTick = [-3 0 3];
% ax.YTick = [-3 0 3];
defaultplot
end