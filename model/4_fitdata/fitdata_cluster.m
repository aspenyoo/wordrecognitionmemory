function [bestFitParam, nLL_est, startTheta, outputt] = fitdata_cluster(isubj, testmodelname, binningfn, fixparams, truemodelname, nConf, nStartVals)

if nargin < 4; fixparams = []; end
if nargin < 5; truemodelname = []; end
if isempty(truemodelname); truemodelname = testmodelname; end
if nargin < 6; nConf = 20; end
if isempty(nConf); nConf = 20; end 
if (size(fixparams,2) > 1) && (size(fixparams,1) < 2) % if it is a vector of Ms, instead of a 2 x fixed parameter things
    nMs = length(fixparams);
    if nargin < 7; nStartVals = 1; end
else
    nMs = 1;
    if nargin < 7; nStartVals = 10; end
end


filepath = 'wordrecognitionmemory/model/4_fitdata/BPSfits/';

% ===== INPUT VARIABLES =====
% ISUBJ: number of subject fitting. can also enter isubj and fixedM as a
% single number (for cluster) to fix M
% OPTIMMETHOD: 'patternbayes','patternsearch' or 'GS'
% FIXPARAMS: 2 x nFixparams matrix if 'patternsearch'.
%            nParams x 3 nGridsVec matrix if 'GS'
%            1 x (number of different Ms) if 'patternbayes'

switch testmodelname
    case {'UVSD','FP','FPheurs'}
        nParams = 2;
    case {'VP','VPheurs'}
        nParams = 3;
    case 'REM'
        nParams = 5;
end

switch binningfn
    case {0,1} % linear, logistic
        % slope, y-int, sigma_mc
        nParams = nParams + 3;
    case 2     % logarithmic
        % a, b, d0, sigma_mc
        nParams = nParams + 4;
    case 3      % power law
        % a, b, gamma, d0, sigma_mc
        nParams = nParams + 5;
    case 4 % weibull
        % a, b, shape, scale, d0, sigma_mc
        nParams = nParams + 6;
end


% loading real subject data
[nnew_part, nold_part] = loadsubjdata(isubj,truemodelname,nConf);

% open txt file
permission = 'a+'; % open or create new file for reading and writing. append data to the end of the file

% parameter fitting function
% iterMs = [];
formatSpec = repmat('%4.4f \t ',1,2*nParams+2);
formatSpec = [formatSpec(1:end-3) '\r\n'];

for iM = 1:nMs
    if (size(fixparams,2) > 1) && (size(fixparams,1) < 2) % if it is a vector of Ms, instead of a 2 x fixed parameter things
        fixparam = [1; fixparams(iM)];
    else
        fixparam = fixparams;
    end
    
    for istartval = 1:nStartVals
        istartval
        
        filename = [filepath 'paramfit_patternbayes_' testmodelname num2str(binningfn) '_subj' num2str(isubj) '.txt'];
        if isubj > 14
            filename = [filepath 'modelrecovery_patternbayes_' testmodelname num2str(binningfn) '_' truemodelname 'subj' num2str(isubj) '.txt'];
        end
        
        [bestFitParam, nLL_est, startTheta, outputt] = paramfit_patternbayes(testmodelname, binningfn, nnew_part, nold_part, fixparam ,1);
        fileID = fopen(filename,permission);
        A1 = [bestFitParam, nLL_est, startTheta, outputt.fsd];
        fprintf(fileID, formatSpec, A1); % save stuff in txt file
        disp('saved')
        fclose(fileID);
                
    end
end

end
