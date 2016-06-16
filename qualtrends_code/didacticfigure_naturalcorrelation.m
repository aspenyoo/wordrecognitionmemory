function didacticfigure_naturalcorrelation

M = 2;
% sigma = 0.1;
N = 10;

% getting new and old words
SNew = randn(N,M);
SOld = randn(N,M);

% idx of minimum distance
distances = squeeze(sum((permute(repmat(SNew,[1,1,N]), [3,2,1])-repmat(SOld,[1,1,N])).^2,2));
mindist = squeeze(min(distances,[],1));
[row,~] = find(distances == repmat(mindist,[N 1]));

% plot the new and old words with the distance lines between themhold on
figure; hold on
plot(SNew(:,1),SNew(:,2),'k.','MarkerSize',12); % plot new words
plot(SOld(:,1),SOld(:,2),'ko','MarkerSize',4); % plot old words
plot([SNew(:,1) SOld(row,1)]',[SNew(:,2) SOld(row,2)]',':k')

axis equal
axismax = 3;
axis([-axismax axismax -axismax axismax])
set(gca,'XTick',[-3 0 3],'YTick',[-3 0 3]);
% ax = gca;
% ax.XTick = [-3 0 3];
% ax.YTick = [-3 0 3];
defaultplot