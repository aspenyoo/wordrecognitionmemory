function create_joblist(jobfilename, jobnumVec, estTimeVec, maxTime, nJobs)
% create_joblist creates a txt file with a list of NJOBS that run for a 
% maximum of MAXTIME
% 
% ========= INPUT VARIABLES ==========
% JOBNUMVEC: vector with length = number of jobs (currently M values in word recogntion memory)
% ESTIMEVEC: vector with length = number of jobs, each element is the 
%   estimated time of corresponding job in JOBNUMVEC
% NJOBS: number of jobs you want to be created (will be overridden in
% MAXTIME * NJOBS < sum(ESTTIMEVEC))
% MAXTIME: maximum time for each job

if isempty(jobfilename); c = clock; jobfilename = ['joblist_' sprintf('%02d%02d%04d',c(2),c(3),c(1)) '.txt']; end
if nargin < 4; maxTime = 48; end % 48 hours per job
if nargin < 5; nJobs = []; end
if ~isempty(maxTime) && isempty(nJobs); % nJob is determed such that each job takes 48 hours
    nJobs = ceil(sum(estTimeVec)/maxTime); 
elseif isempty(maxTime) && ~isempty(nJobs)
    maxTime = ceil(sum(estTimeVec)/nJobs);
end 

combos = sortrows([jobnumVec(:) estTimeVec(:)],2);

% check to make sure that they aren't less than minTime. if so, change
% nJobs
% THIS PART IS NOT DONE YET

% make a paramcombination x nParam+1 matrix, where last column is estTime
jobs = cell(1,nJobs);

% allocating jobs so that they are evenly spaced in terms of time taken
for i = 1:nJobs; % start by giving each joblist one job
    jobs{i} = [jobs{i}; combos(end,:)];
    combos(end,:) = [];
end
while size(combos,1);
    % index of current job taking the least amount of time
    currMin = cellfun(@(x) sum(x(:,2)),jobs,'UniformOutput',false);
    currMin = cell2mat(currMin);
    currMin = find(currMin == min(currMin),1,'first');
    
    % add longest taking param value to current job with min time
    jobs{currMin} = [jobs{currMin}; combos(end,:)];
    
    % delete that param value from combos matrix.
    combos(end,:) = [];
end

% % time it will take on average
% mean(cell2mat(cellfun(@(x) sum(x(:,2)),jobs,'UniformOutput',false)))

% making file
fid = fopen(jobfilename,'w');
for ijob = 1:nJobs;
    currjoblist = jobs{ijob}(:,1);
    textFormat = cell2mat(repmat({'%d \t '},1,length(currjoblist)));
    textFormat = [textFormat(1:end-4) '\r\n'];
    fprintf(fid,textFormat,currjoblist);
end
fclose(fid);

