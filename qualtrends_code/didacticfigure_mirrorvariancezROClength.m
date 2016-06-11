% didacticfigure_mirrorvariancezROClength
% plots a didactic figure of the three qualitative trends that arise from
% two-condition word recognitio memory studies


% plotting actual distributions
xx = linspace(-5,5,100);
strongM = 2;
strongSD = 1;
weakSD = 0.7;
weakM = 0.5;

SO = normpdf(xx,strongM,strongSD);
WO = normpdf(xx,weakM,weakSD);
WN = normpdf(xx,-weakM,weakSD);
SN = normpdf(xx,-strongM,strongSD);

plot(xx,[SO' WO' WN' SN'],'Color','k')
defaultplot


% plotting zROC length effect
criteria = linspace(-4,4,6);
zFA_SN = norminv(1-normcdf(criteria,-strongM,strongSD));
zFA_WN = norminv(1-normcdf(criteria,-weakM,weakSD));
zH_WO = norminv(1-normcdf(criteria,weakM,weakSD));
zH_SO = norminv(1-normcdf(criteria,strongM,strongSD));

figure;
plot(zFA_SN,zH_SO,'k-'); hold on;
plot(zFA_WN,zH_WO,'Color',0.3*ones(1,3))
defaultplot
