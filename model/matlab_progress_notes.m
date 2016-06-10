%% THIS IS A SCRIPT TO DOCUMENT LITTLE NOTES/CODE I RUN THROUGHOUT THE PROJECT
% useful to see what I've done and when I've done it, and to write scratch
% code

%% MARCH 30, 2016
% setting up the cluster to be able to run less-but-longer code

% the below is used to create a joblist given the parameters I want to fit
nStartVals = 10;
jobnumVec = repmat(1:50, 1,nStartVals);
estTimeVec = repmat(linspace(1/6,7,50),1,nStartVals); % in hours
maxTime = 24;
nJobs = [];
create_joblist(jobnumVec, estTimeVec, maxTime, nJobs)
% the joblist created from this has 75 jobs. 

% this will be plugged into my new function cluster_wrapper.m

%% MARCH 31, 2016

% running cluster_wrapper.m on the cluster with paramfit_wrapper.sh

%% APRIL 5, 2015


% check number of completed runs for each subject and make sure that it's
% nStartVals

nStartVals = 10;
nM = 50;
nSubj = 14;
modelname = 'FPheurs';
optimMethod = 'patternbayes';

count = nan(nSubj,nM);
for isubj = 1:nSubj;
    for iM = 1:nM;
        count(isubj,iM) = countnum(modelname,isubj,iM, optimMethod);
    end
end

count

% since counts are equal across Ms
countsleft = 10 - count(:,1);

% =======================================================
% create joblist for all Ms
nStartVals = 1;
jobnumVec = repmat(1:50, 1,nStartVals);
estTimeVec = repmat(linspace(1/6,7,50),1,nStartVals); % in hours
maxTime = 24;
nJobs = [];
create_joblist(jobnumVec, estTimeVec, maxTime, nJobs)

%% APRIL 12, 2016

% increasing Mmax for all subjects
nStartVals = 10;
jobnumVec = repmat(51:75, 1, nStartVals);
estTimeVec = jobnumVec.*1/6; % in hours
maxTime = 48;
nJobs = [];
create_joblist(jobnumVec, estTimeVec, maxTime, nJobs)
