function makejoblist(testmodelname,maxtime, optimMethod,subjids,truemodelname)

if nargin < 3; optimMethod = 'patternbayes'; end
if nargin < 4; subjids = 1:14; end
if nargin < 5; truemodelname = []; end

Mmax = 50;
filepath = 'model/4_fitdata/';
jobtime = 12;

for isubj = 1:length(subjids);
    subjid = subjids(isubj)
    
    jobnumVec = []; esttimeVec = [];
    for iM = 1:Mmax;
        numM = countnum(testmodelname,subjid,iM,truemodelname);

        jobnumVec = [jobnumVec repmat(iM, [1 10-numM])];
        esttimeVec = [esttimeVec repmat(jobtime/numM,[1 10-numM])];
    end
    esttimeVec(isinf(esttimeVec)) = jobtime;
    
    jobfilename = [filepath 'joblist_' optimMethod '_subj' num2str(subjid) '_maxtime' num2str(maxtime) '.txt'];
    create_joblist(jobfilename, jobnumVec, esttimeVec, maxtime)
    
%     % put everything in the new joblist
%     fileID = fopen(jobfilename,permission);
%     A1 = [bestFitParam, nLL_est, startTheta, outputt.fsd];
%     fprintf(fileID, formatSpec, A1); % save stuff in txt file
%     fclose(fileID);
end