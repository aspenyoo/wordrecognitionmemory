function didacticfigure_newwordsfromoldwords

M = 2;
sigma = 0.1;
N = 10;
signewVec = [0.25 0.5];
xx = linspace(0,2*pi,50);

% same old words and noisy memories
X = randn(N,M)*sqrt(sigma^2+1);
SOld = (X/sigma^2)/(1/sigma^2 + 1) + randn(N,M)*(1/sqrt(1/sigma^2 + 1)); % drawing old words from X

figure
for isig = 1:2;
sigma_new = signewVec(isig);

% getting new and old words
SNew = SOld + randn(N,M).*sigma_new;

% idx of minimum distance
distances = squeeze(sum((permute(repmat(SNew,[1,1,N]), [3,2,1])-repmat(SOld,[1,1,N])).^2,2));
mindist = squeeze(min(distances,[],1));
[row,~] = find(distances == repmat(mindist,[N 1]));

% plot the new and old words with the distance lines between them
subplot(1,2,isig); hold on
plot(SNew(:,1),SNew(:,2),'k.','MarkerSize',12); % plot new words
plot(SOld(:,1),SOld(:,2),'ko','MarkerSize',4); % plot old words
% plot([SNew(:,1) SOld(row,1)]',[SNew(:,2) SOld(row,2)]',':k')

% plot circle around each old word with sigma_new
x = cos(xx).*sigma_new;
y = sin(xx).*sigma_new;
plot(bsxfun(@plus,x',SOld(:,1)'),bsxfun(@plus,y',SOld(:,2)'),':','Color',0.6*ones(1,3));

axis equal
axismax = 3;
axis([-axismax axismax -axismax axismax])
set(gca,'XTick',[-3 0 3],'YTick',[-3 0 3]);
% ax = gca;
% ax.XTick = [-3 0 3];
% ax.YTick = [-3 0 3];
defaultplot
end