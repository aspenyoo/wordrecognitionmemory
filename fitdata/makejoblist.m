function makejoblist(testmodelname,binningfn, maxtime, optimMethod,subjids,Mrange,nJobs, truemodelname)

if nargin < 4; optimMethod = 'patternbayes'; end
if nargin < 5; subjids = 1:14; end
if nargin < 6; Mrange = [1 50]; end
if nargin < 7; nJobs = []; end
if nargin < 8; truemodelname = []; end

filepath = 'model/4_fitdata/';
jobtime = 12;

for isubj = 1:length(subjids);
    subjid = subjids(isubj)
    
    jobnumVec = []; esttimeVec = [];
    for iM = Mrange(1):Mrange(2);
        numM = countnum(testmodelname,binningfn,subjid,iM,truemodelname);

        jobnumVec = [jobnumVec repmat(iM, [1 10-numM])];
        esttimeVec = [esttimeVec repmat(jobtime/numM,[1 10-numM])];
    end
    esttimeVec(isinf(esttimeVec)) = jobtime;

    jobfilename = [filepath 'joblist_' testmodelname num2str(binningfn) '_' optimMethod '_subj' num2str(subjid) '_maxtime' num2str(maxtime) '.txt'];
    create_joblist(jobfilename, jobnumVec, esttimeVec, maxtime, nJobs)
    
%     % put everything in the new joblist
%     fileID = fopen(jobfilename,permission);
%     A1 = [bestFitParam, nLL_est, startTheta, outputt.fsd];
%     fprintf(fileID, formatSpec, A1); % save stuff in txt file
%     fclose(fileID);
end