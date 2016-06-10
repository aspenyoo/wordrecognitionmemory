function count = countnum(modelname,isubj,M, optimMethod)
if nargin < 4; optimMethod = 'patternbayes'; end

% modelname = 'FP';
% isubj = 21;
% optimMethod = 'patternbayes';

filename = ['paramfit_' optimMethod '_' modelname '_subj' num2str(isubj) '.txt'];

% ensuring that there are no space in the txt to fuck up the reading of the file
% removetxtspaces(modelname,isubj,optimMethod)

data = dlmread(filename);
count = sum(data(:,1) == M);

