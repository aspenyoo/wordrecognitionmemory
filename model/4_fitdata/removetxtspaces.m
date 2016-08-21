function removetxtspaces(modelname,binningfn,memstrengthvar,isubj,optimMethod,filepath)
if nargin < 5; optimMethod = 'patternbayes'; end
if nargin < 6; filepath = 'model/4_fitdata/BPSfits/'; end

filename = [filepath 'paramfit_' optimMethod '_' modelname num2str(binningfn) ...
    num2str(memstrengthvar) '_subj' num2str(isubj) '.txt'];

% Read the file as cell string line by line
fid = fopen(filename,'r');
% FileName = fopen(fid);
if fid < 0, error('Cannot open file: %s', filename); end
Data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

% remove empty lines
idxsMat = [];
for i = 1:9;
    tempidx = cellfun(@(x) min([find(x == num2str(i)) Inf]),Data{1},'UniformOutput',false);
    idxsMat = [idxsMat tempidx];
end
idxMat = min(cell2mat(idxsMat),[],2);
data = Data{1};

% get rid of any empty lines
data(idxMat == Inf,:) = [];
idxMat(idxMat == Inf) = [];

% get rid of m
idxMat = mat2cell(idxMat,ones(1,length(idxMat)));
data = cellfun(@(x,y) x(y:end),data,idxMat,'UniformOutput',false);

% Write the cell string
fid = fopen(filename, 'w');
fprintf(fid, '%s\r\n', data{:});
fclose(fid);

