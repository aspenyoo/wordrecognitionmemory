function confdistplot(pnew, pold)

gold = aspencolors('gold');
greyblue = aspencolors('greyblue');

plot(1:20,pnew,'.','MarkerSize', 12,'Color',greyblue);
hold on;
plot(1:20,pold,'.','MarkerSize', 12,'Color',gold);
hold off;

defaultplot
set(gca     ,...
    'XTick'     ,0:5:20     );