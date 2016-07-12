function count = countnum(testmodelname,isubj,M, truemodelname, optimMethod)
if nargin < 4; truemodelname = testmodelname; end
if nargin < 5; optimMethod = 'patternbayes'; end

filename = ['paramfit_' optimMethod '_' testmodelname '_subj' num2str(isubj) '.txt'];
if isubj > 14;
    filename = ['modelrecovery_patternbayes_' testmodelname '_' truemodelname 'subj' num2str(isubj) '.txt'];
end

% ensuring that there are no space in the txt to fuck up the reading of the file
% removetxtspaces(testmodelname,isubj,optimMethod)

filepath = 'model/4_fitdata/BPSfits/';
data = dlmread([filepath filename]);
count = sum(data(:,1) == M);

