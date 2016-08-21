function [bestFitParam, nLL_est, startTheta, outputt] = fitdata_cluster(isubj, testmodelname, binningfn, memstrengthvar, optimMethod, fixparams, truemodelname, nConf, nStartVals)

if nargin < 5; optimMethod = 'patternbayes'; end
if nargin < 6; fixparams = []; end
if nargin < 7; truemodelname = []; end
if isempty(truemodelname); truemodelname = testmodelname; end
if nargin < 8; nConf = 20; end
if isempty(nConf); nConf = 20; end 
if (size(fixparams,2) > 1) && (size(fixparams,1) < 2); % if it is a vector of Ms, instead of a 2 x fixed parameter things
    nMs = length(fixparams);
    if nargin < 9; nStartVals = 1; end
else
    nMs = 1;
    if nargin < 9; nStartVals = 10; end
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
    case {'FP','FPheurs','UVSD'}
        nParams = 2;
    case {'VP','VPheurs'}
        nParams = 3;
    case 'REM'
        nParams = 5;
end
switch binningfn
    case {0,1}; % linear, logistic
        nParams = nParams + 3;
    case 2      % logarithmic
        nParams = nParams + 4;
    case 3      % power law
        nParams = nParams + 5; 
end

if strcmp(optimMethod,'GS'); nGridsVec = fixparams; clear fixparams; end

% % used in cluster to indicate both subject and fixed M
% if isubj > 100;
%     blah = num2str(isubj);
%     isubj = str2double(blah(1:end-2));
%     fixparams(2,1) = str2double(blah(end-1:end)); % fixing M value
% end

% loading real subject data
[nnew_part, nold_part] = loadsubjdata(isubj,truemodelname,nConf);

% open txt file
permission = 'a+'; % open or create new file for reading and writing. append data to the end of the file

% parameter fitting function
iterMs = [];
formatSpec = repmat('%4.4f \t ',1,2*nParams+2);
formatSpec = [formatSpec(1:end-3) '\r\n'];

for iM = 1:nMs;
    if (size(fixparams,2) > 1) && (size(fixparams,1) < 2); % if it is a vector of Ms, instead of a 2 x fixed parameter things
        fixparam = [1; fixparams(iM)];
    else
        fixparam = fixparams;
    end
    
    for istartval = 1:nStartVals;
        switch optimMethod
            case 'patternbayes'
                filename = [filepath 'paramfit_patternbayes_' testmodelname num2str(binningfn) num2str(memstrengthvar) '_subj' num2str(isubj) '.txt'];
                if isubj > 14;
                    filename = [filepath 'modelrecovery_patternbayes_' testmodelname num2str(binningfn) num2str(memstrengthvar) '_' truemodelname 'subj' num2str(isubj) '.txt'];
                end
                
                [bestFitParam, nLL_est, startTheta, outputt] = paramfit_patternbayes(testmodelname, binningfn, memstrengthvar, nnew_part, nold_part, fixparam ,1);
                fileID = fopen(filename,permission);
                A1 = [bestFitParam, nLL_est, startTheta, outputt.fsd];
                fprintf(fileID, formatSpec, A1); % save stuff in txt file
                fclose(fileID);
                
                %     case 'patternsearch'
                %         nStartVals = 1;
                %         [bestFitParam, nLL_est, startTheta, outputt] = paramfit_debug(modelname, nnew_part, nold_part, fixparams, nStartVals, nConf);
                %         filename = ['paramfit_patternsearch_' modelname '_subj' num2str(isubj) '_fixparam' num2str(fixparams(1,:)) '_value' num2str(fixparams(2,:)) '_'...
                %             num2str(c(4)) num2str(c(5)) num2str(c(6)) '.mat'];
                %     case 'GS'
                %         nX = 50;
                %         nrep = 50;
                %         [bestFitParam, nLL_est] = paramfit_GS_correct(modelname,nnew_part, nold_part, nGridsVec, nX, nrep, nConf);
                %         startTheta = nan;
                %         filename = ['fitdata_correct_subj' num2str(isubj) '_M' num2str(nGridsVec(1,1)) '_09132015.mat'];
        end
    end
end

end
