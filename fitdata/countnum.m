function count = countnum(testmodelname,binningfn, isubj,M, truemodelname, optimMethod)
if nargin < 5; truemodelname = testmodelname; end
if nargin < 6; optimMethod = 'patternbayes'; end

filename = ['paramfit_' optimMethod '_' testmodelname num2str(binningfn) '_subj' num2str(isubj) '.txt'];
if isubj > 14;
    filename = ['modelrecovery_patternbayes_' testmodelname num2str(binningfn) '_' truemodelname 'subj' num2str(isubj) '.txt'];
end

% ensuring that there are no space in the txt to fuck up the reading of the file
% removetxtspaces(testmodelname,isubj,optimMethod)

filepath = 'model/4_fitdata/BPSfits/';
try
    data = dlmread([filepath filename]);
    count = sum(data(:,1) == M);
catch
    count = 0;
end

